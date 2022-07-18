#pragma once
/*
We do several optimization and correctness patches for CCEH, including:
(1) remove fence between storing value and storing key during insert() because
these two stores are in the same cacheline and will mot be reordered. (2) remove
bucket-level lock described in their original paper since frequent
lock/unlocking will severly degrade its performance (actually their original
open-sourced code also does not have bucket-level lock). (3) add epoch manager
in the application level (mini-benchmark) to gurantee correct memory
reclamation. (4) avoid the perssitent memory leak during the segment split by
storing the newly allocated segment in a small preallocated area (organized as a
hash table). (5) add uniqnuess check during the insert opeartion to avoid
inserting duplicate keys. (6) add support for variable-length key by storing the
pointer to the key object. (7) use persistent lock in PMDK library to aovid
deadlock caused by sudden system failure.
*/
#include <bitset>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <shared_mutex>
#include <thread>
#include <unordered_map>
#include <vector>

#include "../util/hash.h"
#include "../util/pair.h"
#include "../Hash.h"
#include "../PMAllocator.h"
#ifdef PMEM
#include <libpmemobj.h>
#endif

#define PERSISTENT_LOCK 1
#define INPLACE 1
#define EPOCH 1
#define LOG_NUM 1024

namespace cceh {
#define CCEH_HYBRID
// #define SEARCH_CLOCK
// #define THROUGHPUT


/* search clock */
double hash_time;
double read_segment_time;
double lock_time;
double search_time;
double search_in_target_time;
double search_in_neighbor_time;
double init_time;
double double_time;
double get_time;

/* search count */
int target_cnt;
int neighbor_cnt;
int insert_count;

template <class T>
struct _Pair {
  T key;
  Value_t value;
};

struct log_entry {
  uint64_t lock;
  PMEMoid temp;

  void Lock_log() {
    uint64_t temp = 0;
    while (!CAS(&lock, &temp, 1)) {
      temp = 0;
    }
  }

  void Unlock_log() { lock = 0; }
};

// const size_t kCacheLineSize = 64;
constexpr size_t kSegmentBits = 8;
constexpr size_t kMask = (1 << kSegmentBits) - 1;
constexpr size_t kShift = kSegmentBits;
constexpr size_t kSegmentSize = (1 << kSegmentBits) * 16 * 4;
constexpr size_t kNumPairPerCacheLine = kCacheLineSize / 16;
constexpr size_t kNumCacheLine = 4;

uint64_t clflushCount;

inline bool var_compare(char *str1, char *str2, int len1, int len2) {
  if (len1 != len2) return false;
  return !memcmp(str1, str2, len1);
}

template <class T>
struct Segment {
  static const size_t kNumSlot = kSegmentSize / sizeof(_Pair<T>);

  Segment(void)
      : local_depth{0}, sema{0}, count{0}, seg_lock{0}, mutex(), rwlock() {
    memset((void *)&_[0], 255, sizeof(_Pair<T>) * kNumSlot);
  }

  Segment(size_t depth)
      : local_depth{depth}, sema{0}, count{0}, seg_lock{0}, mutex(), rwlock() {
    memset((void *)&_[0], 255, sizeof(_Pair<T>) * kNumSlot);
  }

  static void New(Segment<T> **seg, size_t depth) {
#ifdef PMEM
    auto callback = [](PMEMobjpool *pool, void *ptr, void *arg) {
      auto value_ptr = reinterpret_cast<size_t *>(arg);
      auto seg_ptr = reinterpret_cast<Segment *>(ptr);
      seg_ptr->local_depth = *value_ptr;
      seg_ptr->sema = 0;
      seg_ptr->count = 0;
      seg_ptr->seg_lock = 0;
      memset((void *)&seg_ptr->mutex, 0, sizeof(std::shared_mutex));
      memset((void *)&seg_ptr->rwlock, 0, sizeof(PMEMrwlock));
      memset((void *)&seg_ptr->_[0], 255, sizeof(_Pair<T>) * kNumSlot);
      return 0;
    };
    PMAllocator::Allocate((void**)seg, kCacheLineSize, sizeof(Segment), callback,
                        reinterpret_cast<void *>(&depth));
#else
    PMAllocator::ZAllocate((void **)seg, kCacheLineSize, sizeof(Segment));
    new (*seg) Segment(depth);
#endif
  }

  static void New(PMEMoid *seg, size_t depth) {
#ifdef PMEM
    auto callback = [](PMEMobjpool *pool, void *ptr, void *arg) {
      auto value_ptr = reinterpret_cast<size_t *>(arg);
      auto seg_ptr = reinterpret_cast<Segment *>(ptr);
      seg_ptr->local_depth = *value_ptr;
      seg_ptr->sema = 0;
      seg_ptr->count = 0;
      seg_ptr->seg_lock = 0;
      memset((void *)&seg_ptr->mutex, 0, sizeof(std::shared_mutex));
      memset((void *)&seg_ptr->rwlock, 0, sizeof(PMEMrwlock));
      memset((void *)&seg_ptr->_[0], 255, sizeof(_Pair<T>) * kNumSlot);
      pmemobj_persist(pool, seg_ptr, sizeof(Segment<T>));
      return 0;
    };
    PMAllocator::Allocate(seg, kCacheLineSize, sizeof(Segment), callback,
                        reinterpret_cast<void *>(&depth));
#endif
  }

  ~Segment(void) {}

  int Insert(PMEMobjpool *, T, Value_t, size_t, size_t);
  int Insert4split(T, Value_t, size_t);
  bool Put(T, Value_t, size_t);
  PMEMoid *Split(PMEMobjpool *, size_t, log_entry *);

  void get_lock(PMEMobjpool *pop) {
#ifdef PERSISTENT_LOCK
    pmemobj_rwlock_wrlock(pop, &rwlock);
#else
    mutex.lock();
#endif
  }

  void release_lock(PMEMobjpool *pop) {
#ifdef PERSISTENT_LOCK
    pmemobj_rwlock_unlock(pop, &rwlock);
#else
    mutex.unlock();
#endif
  }

  void get_rd_lock(PMEMobjpool *pop) {
#ifdef PERSISTENT_LOCK
    pmemobj_rwlock_rdlock(pop, &rwlock);
#else
    mutex.lock_shared();
#endif
  }

  void release_rd_lock(PMEMobjpool *pop) {
#ifdef PERSISTENT_LOCK
    pmemobj_rwlock_unlock(pop, &rwlock);
#else
    mutex.unlock_shared();
#endif
  }

  bool try_get_lock(PMEMobjpool *pop) {
#ifdef PERSISTENT_LOCK
    if (pmemobj_rwlock_trywrlock(pop, &rwlock) == 0) {
      return true;
    }
    return false;
#else
    return mutex.try_lock();
#endif
  }

  bool try_get_rd_lock(PMEMobjpool *pop) {
#ifdef PERSISTENT_LOCK
    if (pmemobj_rwlock_tryrdlock(pop, &rwlock) == 0) {
      return true;
    }
    return false;
#else
    return mutex.try_lock_shared();
#endif
  }

  _Pair<T> _[kNumSlot];
  size_t local_depth;
  int64_t sema = 0;
  size_t pattern = 0;
  int count = 0;
  std::shared_mutex mutex;
  uint64_t seg_lock;
  PMEMrwlock rwlock;
};

template <class T>
struct Seg_array {
  typedef Segment<T> *seg_p;
  size_t global_depth;
  seg_p _[0];

  static void New(PMEMoid *sa, size_t capacity) {
#ifdef PMEM
    auto callback = [](PMEMobjpool *pool, void *ptr, void *arg) {
      auto value_ptr = reinterpret_cast<size_t *>(arg);
      auto sa_ptr = reinterpret_cast<Seg_array *>(ptr);
      sa_ptr->global_depth = static_cast<size_t>(log2(*value_ptr));
      memset(sa_ptr->_, 0, (*value_ptr) * sizeof(uint64_t));
      return 0;
    };
    PMAllocator::Allocate(sa, kCacheLineSize,
                        sizeof(Seg_array) + sizeof(uint64_t) * capacity,
                        callback, reinterpret_cast<void *>(&capacity));
#else
    PMAllocator::ZAllocate((void **)sa, kCacheLineSize, sizeof(Seg_array));
    new (*sa) Seg_array(capacity);
#endif
  }

  static void NewHybrid(Seg_array **sa, size_t capacity)
  {
    *sa = (Seg_array *)malloc(sizeof(Seg_array) + sizeof(uint64_t) * capacity);
    Seg_array *sa_ptr = *sa;
    sa_ptr->global_depth = log2(capacity);
    memset(sa_ptr->_, 0, (capacity) * sizeof(uint64_t));
  }
};

template <class T>
struct Directory {
  static const size_t kDefaultDirectorySize = 1024;
  Seg_array<T> *sa;
  PMEMoid new_sa;
  size_t capacity;
  bool lock;
  int sema = 0;

  Directory(Seg_array<T> *_sa) {
    capacity = kDefaultDirectorySize;
    sa = _sa;
    new_sa = OID_NULL;
    lock = false;
    sema = 0;
  }

  Directory(size_t size, Seg_array<T> *_sa) {
    capacity = size;
    sa = _sa;
    new_sa = OID_NULL;
    lock = false;
    sema = 0;
  }

  static void New(Directory **dir, size_t capacity) {
#ifdef PMEM
    auto callback = [](PMEMobjpool *pool, void *ptr, void *arg) {
      auto value_ptr =
          reinterpret_cast<std::pair<size_t, Seg_array<T> *> *>(arg);
      auto dir_ptr = reinterpret_cast<Directory *>(ptr);
      dir_ptr->capacity = value_ptr->first;
      dir_ptr->sa = value_ptr->second;
      dir_ptr->new_sa = OID_NULL;
      dir_ptr->lock = false;
      dir_ptr->sema = 0;
      dir_ptr = nullptr;
      return 0;
    };
    auto call_args = std::make_pair(capacity, nullptr);
    PMAllocator::Allocate((void **)dir, kCacheLineSize, sizeof(Directory),
                        callback, reinterpret_cast<void *>(&call_args));
#else
    PMAllocator::ZAllocate((void **)dir, kCacheLineSize, sizeof(Directory));
    new (*dir) Directory(capacity, temp_sa);
#endif
  }

#ifdef CCEH_HYBRID
  static void NewHybrid(Directory **dir, size_t capacity)
  {
    *dir = (Directory *)malloc(sizeof(Directory));
    Directory *dir_ptr = *dir;
    dir_ptr->capacity = capacity;
    dir_ptr->sa = nullptr;
    dir_ptr->new_sa = OID_NULL;
    dir_ptr->lock = false;
    dir_ptr->sema = 0;
  }
#endif

  ~Directory(void) {}

  void cal_load_factory()
  {
    size_t count = 0;
    size_t seg_num = 0;
    Seg_array<T> *seg = sa;
    Segment<T> **dir_entry = seg->_;
    Segment<T> *ss;
    auto global_depth = seg->global_depth;
    size_t depth_diff;
    for (int i = 0; i < capacity;) {
      ss = dir_entry[i];
      depth_diff = global_depth - ss->local_depth;

      for (unsigned i = 0; i < Segment<T>::kNumSlot; ++i) {
        if constexpr (std::is_pointer_v<T>) {
          if ((ss->_[i].key != (T)INVALID) &&
              ((h(ss->_[i].key->key, ss->_[i].key->length) >>
                (64 - ss->local_depth)) == ss->pattern)) {
            ++count;
          }
        } else {
          if ((ss->_[i].key != (T)INVALID) &&
              ((h(&ss->_[i].key, sizeof(Key_t)) >> (64 - ss->local_depth)) ==
               ss->pattern)) {
            ++count;
          }
        }
      }

      seg_num++;
      i += pow(2, depth_diff);
    }
    std::cout << "load_factor: " << (double)count / (seg_num * 256 * 4) << std::endl;
  }

  void get_item_num() {
    size_t count = 0;
    size_t seg_num = 0;
    Seg_array<T> *seg = sa;
    Segment<T> **dir_entry = seg->_;
    Segment<T> *ss;
    auto global_depth = seg->global_depth;
    size_t depth_diff;
    for (int i = 0; i < capacity;) {
      ss = dir_entry[i];
      depth_diff = global_depth - ss->local_depth;

      for (unsigned i = 0; i < Segment<T>::kNumSlot; ++i) {
        if constexpr (std::is_pointer_v<T>) {
          if ((ss->_[i].key != (T)INVALID) &&
              ((h(ss->_[i].key->key, ss->_[i].key->length) >>
                (64 - ss->local_depth)) == ss->pattern)) {
            ++count;
          }
        } else {
          if ((ss->_[i].key != (T)INVALID) &&
              ((h(&ss->_[i].key, sizeof(Key_t)) >> (64 - ss->local_depth)) ==
               ss->pattern)) {
            ++count;
          }
        }
      }

      seg_num++;
      i += pow(2, depth_diff);
    }
    std::cout << "#items: " << count << std::endl;
    std::cout << "load_factor: " << (double)count / (seg_num * 256 * 4) << std::endl;
    std::cout << "========================== time ==========================" << std::endl;
    std::cout << "hash time: " << hash_time / CLOCKS_PER_SEC << std::endl;
    std::cout << "segment time: " << read_segment_time / CLOCKS_PER_SEC << std::endl;
    std::cout << "lock time: " << lock_time / CLOCKS_PER_SEC << std::endl;
    std::cout << "search in target time: " << search_in_target_time / CLOCKS_PER_SEC << std::endl;
    std::cout << "search in target count: " << target_cnt << std::endl;
    std::cout << "search in target time: " << search_in_neighbor_time / CLOCKS_PER_SEC << std::endl;
    std::cout << "search in neighbor time: " << neighbor_cnt << std::endl;
    std::cout << "init time: " << init_time / CLOCKS_PER_SEC << std::endl;
    std::cout << "double time: " << double_time / CLOCKS_PER_SEC << std::endl;
    std::cout << "========================== time ==========================" << std::endl;
  }

  bool Acquire(void) {
    bool unlocked = false;
    return CAS(&lock, &unlocked, true);
  }

  bool Release(void) {
    bool locked = true;
    return CAS(&lock, &locked, false);
  }

  void SanityCheck(void *);
};
#ifdef CCEH_HYBRID
  Directory<uint64_t> *dir;
#endif
template <class T>
class CCEH : public Hash<T> {
 public:
  CCEH(void);
  CCEH(int, PMEMobjpool *_pool);
  ~CCEH(void);
  int Insert(T key, Value_t value);
  int Insert(T key, Value_t value, bool);
  bool Delete(T);
  bool Delete(T, bool);
  Value_t Get(T);
  Value_t Get(T, bool is_in_epoch);
  Value_t FindAnyway(T);
  double Utilization(void);
  size_t Capacity(void);
  void Recovery(void);
  void Directory_Doubling(int x, Segment<T> *s0, 
                      #ifdef CCEH_HYBRID
                        Segment<T> *s1
                      #else
                        PMEMoid *s1
                      #endif
                      );
                      
  void Directory_Update(int x, Segment<T> *s0, 
                      #ifdef CCEH_HYBRID
                        Segment<T> *s1
                      #else
                        PMEMoid *s1
                      #endif
                      );
  void Lock_Directory();
  void Unlock_Directory();
  void TX_Swap(void **entry, PMEMoid *new_seg);
  void getNumber() { dir->get_item_num(); }
#ifndef CCEH_HYBRID
  Directory<T> *dir;
#endif
  log_entry log[LOG_NUM];
  int seg_num;
  int restart;
#ifdef PMEM
  PMEMobjpool *pool_addr;
#endif
};
//#endif  // EXTENDIBLE_PTR_H_

template <class T>
int Segment<T>::Insert(PMEMobjpool *pool_addr, T key, Value_t value, size_t loc,
                       size_t key_hash) {
  if (sema == -1) {
    return 2;
  };
  get_lock(pool_addr);
  if ((key_hash >> (8 * sizeof(key_hash) - local_depth)) != pattern ||
      sema == -1) {
    release_lock(pool_addr);
    return 2;
  }
  int ret = 1;
  T LOCK = (T)INVALID;

  /*uniqueness check*/
  auto slot = loc;
  for (unsigned i = 0; i < kNumCacheLine * kNumPairPerCacheLine; ++i) {
    slot = (loc + i) % kNumSlot;
    if constexpr (std::is_pointer_v<T>) {
      if (_[slot].key != (T)INVALID &&
          (var_compare(key->key, _[slot].key->key, key->length,
                       _[slot].key->length))) {
        release_lock(pool_addr);
        return -3;
      }
    } else {
      if (_[slot].key == key) {
        release_lock(pool_addr);
        return -3;
      }
    }
  }

  for (unsigned i = 0; i < kNumPairPerCacheLine * kNumCacheLine; ++i) {
    slot = (loc + i) % kNumSlot;
    if constexpr (std::is_pointer_v<T>) {
      if ((_[slot].key != (T)INVALID) &&
          ((h(_[slot].key->key, _[slot].key->length) >>
            (8 * sizeof(key_hash) - local_depth)) != pattern)) {
        _[slot].key = (T)INVALID;
      }
      if (CAS(&_[slot].key, &LOCK, SENTINEL)) {
        _[slot].value = value;
        _[slot].key = key;
        PMAllocator::Persist(&_[slot], sizeof(_Pair<T>));
        ret = 0;
        break;
      } else {
        LOCK = (T)INVALID;
      }
    } else {
      if ((h(&_[slot].key, sizeof(Key_t)) >>
           (8 * sizeof(key_hash) - local_depth)) != pattern) {
        _[slot].key = INVALID;
      }
      if (CAS(&_[slot].key, &LOCK, SENTINEL)) {
        _[slot].value = value;
        _[slot].key = key;
        PMAllocator::Persist(&_[slot], sizeof(_Pair<T>));
        ret = 0;
        break;
      } else {
        LOCK = INVALID;
      }
    }
  }
  release_lock(pool_addr);
  return ret;
}

template <class T>
int Segment<T>::Insert4split(T key, Value_t value, size_t loc) {
  for (unsigned i = 0; i < kNumPairPerCacheLine * kNumCacheLine; ++i) {
    auto slot = (loc + i) % kNumSlot;
    if (_[slot].key == (T)INVALID) {
      _[slot].key = key;
      _[slot].value = value;
      return 0;
    }
  }
  return -1;
}

template <class T>
PMEMoid *Segment<T>::Split(PMEMobjpool *pool_addr, size_t key_hash,
                           log_entry *log) {
  using namespace std;
  if (!try_get_lock(pool_addr)) {
    return nullptr;
  }
  sema = -1;

  size_t new_pattern = (pattern << 1) + 1;
  size_t old_pattern = pattern << 1;

  uint64_t log_pos = key_hash % LOG_NUM;
  log[log_pos].Lock_log();

  Segment::New(&log[log_pos].temp, local_depth + 1);
  Segment<T> *split =
      reinterpret_cast<Segment<T> *>(pmemobj_direct(log[log_pos].temp));

  for (unsigned i = 0; i < kNumSlot; ++i) {
    uint64_t key_hash;
    if constexpr (std::is_pointer_v<T>) {
      if (_[i].key != (T)INVALID) {
        key_hash = h(_[i].key->key, _[i].key->length);
      }
    } else {
      key_hash = h(&_[i].key, sizeof(Key_t));
    }
    if ((_[i].key != (T)INVALID) &&
        (key_hash >> (8 * 8 - local_depth - 1) == new_pattern)) {
      split->Insert4split(_[i].key, _[i].value,
                          (key_hash & kMask) * kNumPairPerCacheLine);
      if constexpr (std::is_pointer_v<T>) {
        _[i].key = (T)INVALID;
      }
    }
  }

#ifdef PMEM
  PMAllocator::Persist(split, sizeof(Segment<T>));
#endif
  if constexpr (std::is_pointer_v<T>) {
#ifdef PMEM
    PMAllocator::Persist(this, sizeof(Segment<T>));
#endif
  }
  return &log[log_pos].temp;
}

template <class T>
CCEH<T>::CCEH(int initCap, PMEMobjpool *_pool) {
  timeval start, stop;
  gettimeofday(&start, NULL);
#ifdef CCEH_HYBRID
  Directory<T>::NewHybrid(&dir, initCap);
  Seg_array<T>::NewHybrid(&dir->sa, initCap);
#else
  Directory<T>::New(&dir, initCap);
  Seg_array<T>::New(&dir->new_sa, initCap);
  dir->sa = reinterpret_cast<Seg_array<T> *>(pmemobj_direct(dir->new_sa));
#endif
  dir->new_sa = OID_NULL;
  auto dir_entry = dir->sa->_;
  for (int i = 0; i < dir->capacity; ++i) {
    Segment<T>::New(&dir_entry[i], dir->sa->global_depth);
    dir_entry[i]->pattern = i;
  }
  /*clear the log area*/
  for (int i = 0; i < LOG_NUM; ++i) {
    log[i].lock = 0;
    log[i].temp = OID_NULL;
  }

  seg_num = 0;
  restart = 0;
  pool_addr = _pool;
  gettimeofday(&stop, NULL);
  init_time += (double)(stop.tv_usec - start.tv_usec) + 
              (double)(stop.tv_sec - start.tv_sec) * 1000000;
}

template <class T>
CCEH<T>::CCEH(void) {
  std::cout << "Reintialize Up for CCEH" << std::endl;
}

template <class T>
CCEH<T>::~CCEH(void) {}

template <class T>
void CCEH<T>::Recovery(void) {
  PMAllocator::EpochRecovery();
  for (int i = 0; i < LOG_NUM; ++i) {
    if (!OID_IS_NULL(log[i].temp)) {
      pmemobj_free(&log[i].temp);
    }
  }

  if (dir != nullptr) {
    dir->lock = 0;
    if (!OID_IS_NULL(dir->new_sa)) {
      pmemobj_free(&dir->new_sa);
    }

    if (dir->sa == nullptr) return;
    auto dir_entry = dir->sa->_;
    size_t global_depth = dir->sa->global_depth;
    size_t depth_cur, buddy, stride, i = 0;
    /*Recover the Directory*/
    size_t seg_count = 0;
    while (i < dir->capacity) {
      auto target = dir_entry[i];
      depth_cur = target->local_depth;
      target->sema = 0;
      stride = pow(2, global_depth - depth_cur);
      buddy = i + stride;
      for (int j = buddy - 1; j > i; j--) {
        target = dir_entry[j];
        if (dir_entry[j] != dir_entry[i]) {
          dir_entry[j] = dir_entry[i];
          target->pattern = i >> (global_depth - depth_cur);
        }
      }
      seg_count++;
      i = i + stride;
    }
  }
}

template <class T>
void CCEH<T>::TX_Swap(void **entry, PMEMoid *new_seg) {
  TX_BEGIN(pool_addr) {
    pmemobj_tx_add_range_direct(entry, sizeof(void *));
    pmemobj_tx_add_range_direct(new_seg, sizeof(PMEMoid));
    *entry = pmemobj_direct(*new_seg);
    *new_seg = OID_NULL;
  }
  TX_ONABORT { std::cout << "Error in TXN Swap!" << std::endl; }
  TX_END
}

template <class T>
void CCEH<T>::Directory_Doubling(int x, Segment<T> *s0, 
                              #ifdef CCEH_HYBRID
                                Segment<T> *s1
                              #else  
                                PMEMoid *s1
                              #endif
                              ) {
  timeval stop, start;
  gettimeofday(&start, NULL);
  Seg_array<T> *sa = dir->sa;
  Segment<T> **d = sa->_;
  auto global_depth = sa->global_depth;

  /* new segment array*/
#ifdef CCEH_HYBRID
  Seg_array<T>::NewHybrid(&dir->sa, 2 * dir->capacity);
  auto dd = dir->sa->_;
#else
  Seg_array<T>::New(&dir->new_sa, 2 * dir->capacity);
  auto new_seg_array =
      reinterpret_cast<Seg_array<T> *>(pmemobj_direct(dir->new_sa));
  auto dd = new_seg_array->_;
#endif
  for (unsigned i = 0; i < dir->capacity; ++i) {
    dd[2 * i] = d[i];
    dd[2 * i + 1] = d[i];
  }
#ifdef CCEH_HYBRID
  dd[2 * x + 1] = s1;
#else
  TX_Swap((void **)&dd[2 * x + 1], s1);
#endif

#if (defined PMEM) && (! defined CCEH_HYBRID) 
  PMAllocator::Persist(
      new_seg_array,
      sizeof(Seg_array<T>) + sizeof(Segment<T> *) * 2 * dir->capacity);
#endif

#ifdef CCEH_HYBRID
  free(sa);
  dir->new_sa = OID_NULL;
  dir->capacity *= 2;
#else
  auto reserve_item = PMAllocator::ReserveItem();
  TX_BEGIN(pool_addr) {
    pmemobj_tx_add_range_direct(reserve_item, sizeof(*reserve_item));
    // pmemobj_tx_add_range_direct(&dir->sa, sizeof(dir->sa));
    // pmemobj_tx_add_range_direct(&dir->new_sa, sizeof(dir->new_sa));
    // pmemobj_tx_add_range_direct(&dir->capacity, sizeof(dir->capacity));
    PMAllocator::Free(reserve_item, sa);
    dir->sa = reinterpret_cast<Seg_array<T> *>(pmemobj_direct(dir->new_sa));
    dir->new_sa = OID_NULL;
    dir->capacity *= 2;
  }
  TX_ONABORT {
    std::cout << "TXN fails during doubling directory" << std::endl;
  }
  TX_END
#endif
  gettimeofday(&stop, NULL);
  double_time += (double)(stop.tv_usec - start.tv_usec) +
                (double)(stop.tv_sec - start.tv_sec) * 1000000;
}

template <class T>
void CCEH<T>::Lock_Directory() {
  while (!dir->Acquire()) {
    asm("nop");
  }
}

template <class T>
void CCEH<T>::Unlock_Directory() {
  while (!dir->Release()) {
    asm("nop");
  }
}

template <class T>
void CCEH<T>::Directory_Update(int x, Segment<T> *s0, 
                      #ifdef CCEH_HYBRID
                        Segment<T> *s1
                      #else
                        PMEMoid *s1
                      #endif
                      ) {
  Segment<T> **dir_entry = dir->sa->_;
  auto global_depth = dir->sa->global_depth;
  unsigned depth_diff = global_depth - s0->local_depth;
  if (depth_diff == 1) {
    if (x % 2 == 0) {
#ifdef CCEH_HYBRID
      dir_entry[x + 1] = s1;
#else
      TX_Swap((void **)&dir_entry[x + 1], s1);
#endif
#ifdef PMEM
      PMAllocator::Persist(&dir_entry[x + 1], sizeof(Segment<T> *));
#endif
    } else {
#ifdef CCEH_HYBRID
      dir_entry[x] = s1;
#else
      TX_Swap((void **)&dir_entry[x], s1);
#endif
#ifdef PMEM
      PMAllocator::Persist(&dir_entry[x], sizeof(Segment<T> *));
#endif
    }
  } else {
    int chunk_size = pow(2, global_depth - (s0->local_depth));
    x = x - (x % chunk_size);
    int base = chunk_size / 2;
#ifdef CCEH_HYBRID
    dir_entry[x + base + base - 1] = s1;
#else
    TX_Swap((void **)&dir_entry[x + base + base - 1], s1);
#endif
    auto seg_ptr = dir_entry[x + base + base - 1];
    for (int i = base - 2; i >= 0; --i) {
      dir_entry[x + base + i] = seg_ptr;
      PMAllocator::Persist(&dir_entry[x + base + i], sizeof(uint64_t));
    }
  }
}

template <class T>
int CCEH<T>::Insert(T key, Value_t value, bool is_in_epoch) {
  if (!is_in_epoch) {
    // auto epoch_guard = PMAllocator::AquireEpochGuard();
    return Insert(key, value);
  }
  // if(insert_count % 1000000 == 0)
  // {
  //   dir->cal_load_factory();    
  // }
  // insert_count ++;
  return Insert(key, value);
}

template <class T>
int CCEH<T>::Insert(T key, Value_t value) {
STARTOVER:
  uint64_t key_hash;
  if constexpr (std::is_pointer_v<T>) {
    key_hash = h(key->key, key->length);
  } else {
    key_hash = h(&key, sizeof(key));
  }
  auto y = (key_hash & kMask) * kNumPairPerCacheLine;

RETRY:
  auto old_sa = dir->sa;
  auto x = (key_hash >> (8 * sizeof(key_hash) - old_sa->global_depth));
  auto dir_entry = old_sa->_;
  Segment<T> *target = dir_entry[x];
  if (old_sa != dir->sa) {
    goto RETRY;
  }
  auto ret = target->Insert(pool_addr, key, value, y, key_hash);

  if(ret == -3) return -1;

  if (ret == 1) {
    auto s = target->Split(pool_addr, key_hash, log);
    if (s == nullptr) {
      goto RETRY;
    }

    auto ss = reinterpret_cast<Segment<T> *>(pmemobj_direct(*s));
    ss->pattern =
        ((key_hash >> (8 * sizeof(key_hash) - ss->local_depth + 1)) << 1) + 1;
    PMAllocator::Persist(&ss->pattern, sizeof(ss->pattern));

    // Directory management
    Lock_Directory();
    {  // CRITICAL SECTION - directory update
      auto sa = dir->sa;
      dir_entry = sa->_;

      x = (key_hash >> (8 * sizeof(key_hash) - sa->global_depth));
      target = dir_entry[x];
      if (target->local_depth < sa->global_depth) {
        // std::cout << "local double" << std::endl;
        #ifdef CCEH_HYBRID
        Directory_Update(x, target, ss);
        #else
        Directory_Update(x, target, s);
        #endif
      } else {  // directory doubling
        std::cout << "global double" << std::endl;
        dir->cal_load_factory();   
      #ifdef CCEH_HYBRID
        Directory_Doubling(x, target, ss);
      #else
        Directory_Doubling(x, target, s);
      #endif
      }
      target->pattern =
          (key_hash >> (8 * sizeof(key_hash) - target->local_depth)) << 1;
      PMAllocator::Persist(&target->pattern, sizeof(target->pattern));
      target->local_depth += 1;
      PMAllocator::Persist(&target->local_depth, sizeof(target->local_depth));
#ifdef INPLACE
      target->sema = 0;
      target->release_lock(pool_addr);
#endif
    }  // End of critical section
    Unlock_Directory();
    uint64_t log_pos = key_hash % LOG_NUM;
    log[log_pos].Unlock_log();
    goto RETRY;
  } else if (ret == 2) {
    goto STARTOVER;
  }

  return 0;
}

template <class T>
bool CCEH<T>::Delete(T key, bool is_in_epoch) {
  if (!is_in_epoch) {
    // auto epoch_guard = PMAllocator::AquireEpochGuard();
    return Delete(key);
  }
  return Delete(key);
}

template <class T>
bool CCEH<T>::Delete(T key) {
  uint64_t key_hash;
  if constexpr (std::is_pointer_v<T>) {
    key_hash = h(key->key, key->length);
  } else {
    key_hash = h(&key, sizeof(key));
  }
  auto y = (key_hash & kMask) * kNumPairPerCacheLine;

RETRY:
  auto old_sa = dir->sa;
  auto x = (key_hash >> (8 * sizeof(key_hash) - old_sa->global_depth));
  auto dir_entry = old_sa->_;
  Segment<T> *dir_ = dir_entry[x];

  auto sema = dir_->sema;
  if (sema == -1) {
    goto RETRY;
  }
  dir_->get_lock(pool_addr);

  if ((key_hash >> (8 * sizeof(key_hash) - dir_->local_depth)) !=
          dir_->pattern ||
      dir_->sema == -1) {
    dir_->release_lock(pool_addr);
    goto RETRY;
  }

  for (unsigned i = 0; i < kNumPairPerCacheLine * kNumCacheLine; ++i) {
    auto slot = (y + i) % Segment<T>::kNumSlot;
    if constexpr (std::is_pointer_v<T>) {
      if ((dir_->_[slot].key != (T)INVALID) &&
          (var_compare(key->key, dir_->_[slot].key->key, key->length,
                       dir_->_[slot].key->length))) {
        dir_->_[slot].key = (T)INVALID;
        PMAllocator::Persist(&dir_->_[slot], sizeof(_Pair<T>));
        dir_->release_lock(pool_addr);
        return true;
      }
    } else {
      if (dir_->_[slot].key == key) {
        dir_->_[slot].key = (T)INVALID;
        PMAllocator::Persist(&dir_->_[slot], sizeof(_Pair<T>));
        dir_->release_lock(pool_addr);
        return true;
      }
    }
  }
  dir_->release_lock(pool_addr);
  return false;
}

template <class T>
Value_t CCEH<T>::Get(T key, bool is_in_epoch) {
  if (is_in_epoch) {
#ifdef EPOCH
    // auto epoch_guard = PMAllocator::AquireEpochGuard();
#endif
#ifdef THROUGHPUT
    timeval start, stop;
    gettimeofday(&start, NULL);
#endif
    auto value = Get(key);
#ifdef THROUGHPUT
    gettimeofday(&stop, NULL);
    get_time += (double)(stop.tv_usec - start.tv_usec) + 
                (stop.tv_sec - start.tv_sec) * 1000000;
#endif
    return value;
  }
  return Get(key);
}

template <class T>
Value_t CCEH<T>::Get(T key) {
#ifdef SEARCH_CLOCK
  timeval start, stop;
  gettimeofday(&start, NULL);
#endif
  uint64_t key_hash;
  if constexpr (std::is_pointer_v<T>) {
    key_hash = h(key->key, key->length);
  } else {
    key_hash = h(&key, sizeof(key));
  }
  auto y = (key_hash & kMask) * kNumPairPerCacheLine;
#ifdef SEARCH_CLOCK
  gettimeofday(&stop, NULL);
  hash_time += (double)(stop.tv_usec - start.tv_usec) +
                (double) (stop.tv_sec - start.tv_sec) * 1000000;
  gettimeofday(&start, NULL);
#endif
RETRY:
  /*
    x: segment index
    y: bucket index
  */

  auto old_sa = dir->sa;
  auto x = (key_hash >> (8 * sizeof(key_hash) - old_sa->global_depth));
  auto dir_entry = old_sa->_;
  Segment<T> *dir_ = dir_entry[x];
#ifdef SEARCH_CLOCK
  gettimeofday(&stop, NULL);
  read_segment_time += (double)(stop.tv_usec - start.tv_usec) +
                (double) (stop.tv_sec - start.tv_sec) * 1000000;
  gettimeofday(&start, NULL);
#endif
  if (!dir_->try_get_rd_lock(pool_addr)) {
    goto RETRY;
  }

  if ((key_hash >> (8 * sizeof(key_hash) - dir_->local_depth)) !=
          dir_->pattern ||
      dir_->sema == -1) {
    dir_->release_rd_lock(pool_addr);
    goto RETRY;
  }
#ifdef SEARCH_CLOCK
  gettimeofday(&stop, NULL);
  lock_time += (double)(stop.tv_usec - start.tv_usec) +
                (double) (stop.tv_sec - start.tv_sec) * 1000000;
  gettimeofday(&start, NULL);
#endif
  unsigned i;
  Value_t value = nullptr;
  for (i = 0; i < kNumPairPerCacheLine * kNumCacheLine; ++i) {
    auto slot = (y + i) % Segment<T>::kNumSlot;

    _Pair<T> kvs = dir_->_[slot];
// #ifdef SEARCH_CLOCK
//     gettimeofday(&stop, NULL);
//     search_time += (double)(stop.tv_usec - start.tv_usec) +
//                 (double) (stop.tv_sec - start.tv_sec) * 1000000;
// #endif
    if constexpr (std::is_pointer_v<T>) {
      if ((kvs.key != (T)INVALID) &&
          (var_compare(key->key, kvs.key->key, key->length,
                       kvs.key->length))) {
        value = kvs.value;
        dir_->release_rd_lock(pool_addr);
        return value;
      }
    } else {
      if (kvs.key == key) {
        value = kvs.value;
        dir_->release_rd_lock(pool_addr);
        break;
      }
    }
  }

  #ifdef SEARCH_CLOCK
  if(i <= 4)
  {
    target_cnt ++;
    gettimeofday(&stop, NULL);
    search_in_target_time += (double)(stop.tv_usec - start.tv_usec) +
            (double) (stop.tv_sec - start.tv_sec) * 1000000;

  }
  else
  {
    neighbor_cnt ++;
    gettimeofday(&stop, NULL);
    search_in_neighbor_time += (double)(stop.tv_usec - start.tv_usec) +
            (double) (stop.tv_sec - start.tv_sec) * 1000000;
  }
#endif
  if(value != nullptr)
    return value;
  dir_->release_rd_lock(pool_addr);
  return NONE;
}
}  // namespace CCEH