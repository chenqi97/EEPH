#ifndef EEPH_H
#define EEPH_H

#include <cstring>
#include <random>
#include <set>
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifdef PMEM
#include <libpmemobj.h>
#endif
#include "immintrin.h"
#include "murmur3.h"
#include "PMAllocator.h"
#include "Hash.h"

// xxhash
#define XXH_INLINE_ALL
#include "xxhash.h"

// begin kv pair
#define KEY_LEN          8
#define VAL_LEN          8
#define KV_NUM    30000000

#define CACHELINESIZE   64
#define BUCKETS_IN_LEVEL 4

#define MAX_LAYER       12
#define M               32  // max number of hash functions in each layer
#define BUCKET_CAPACITY 11  // item number in a bucket
#define MULTI_REHASH    5   // the number of muti-rehash
#define INIT_LEVEL      32  // initial bucket levels
#define SEED            0   // seed index

const uint32_t lock_pos = ((uint32_t)1 << 31);
const uint32_t lock_mask = ((uint32_t)1 << 31) - 1;
const uint32_t fp_mask = (1 << 8) - 1;
const uint64_t key_mask = ((uint64_t)1 << 8)  - 1;
namespace eeph 
{
// switches
// #define INSERT_CLOCK
// #define SEARCH_CLOCK
// #define DELETE_CLOCK
// #define DIRECTORY
// #define USING_SIMD
// #define XXHASH
#define RECOVERY

// optimes
// #define THRESHOLD
#define USING_SSE
#define LOOP_UNROLLING

/* search clock */
double cell_time;
double init_guide_time;
double bucket_time;
double search_time;
double hash_time;
double level_time;
double level_1_time;
double guide_1_time;
double read_bucket_time;

/* insert clock */
double ihash_time;
double iread_cell_time;
double iread_cal_time;
double icalculate_bucket_time;
double insert_time;
double insert_failed_time;
double collect_time;
double threshold_time;
double move_time;
double extend_time;
    double double_bucket_time;
    double copy_bucket_time;
    double double_cell_time;
int flag = 0;
int bucket_flag = 1;

/* delete clock */
double dcell_time;
double dinit_guide_time;
double dbucket_time;
double dsearch_time;
double dohter_time;

//debug
int query_cnt = 0;
int cnt = 0; // debug use
int insert_dire_cnt = 0;
int move_cnt = 0;
int ins_cnt = 0;

int extend_cnt;
int rehash_cnt;
int rehash_muti_cnt;
int threshold_cnt;
int insert_count;
// exp
int seq_cnt = 0;
int seq_bucket = 0;


int bucket_items;               // the number of items in sys
int capacity;                   // the capacity of the sys
int bucket_items_per_level[256]; // the number of items in per level


template <class T>
struct _Pair {
  T key;
  Value_t value;
};
template <class T>
struct Bucket{
    
    uint32_t version_lock;                      // lock & version lock
    uint32_t k[BUCKET_CAPACITY];                // item's cell index
    _Pair<T> kv[BUCKET_CAPACITY];               // kv
    uint16_t bitmap;                            // valid bitmap
    uint8_t used;                               // used slots
    uint8_t init_guide[BUCKET_CAPACITY];
    uint8_t fp[BUCKET_CAPACITY];                 // fingerprint
    // uint8_t offset;
    uint8_t padding[7];                         // padding

    static void New(Bucket<T> **bkt, uint32_t bucket_number, size_t arg) {
        auto init_callback = [](PMEMobjpool *pool, void *ptr, void *arg) {
            auto bucket_p = reinterpret_cast<Bucket *>(ptr);
            uint32_t bucket_number__ = *(reinterpret_cast<uint32_t *>(arg));
            for(int i = 0; i < bucket_number__; i ++)
            {
                auto tmp_bucket = bucket_p + i;
                memset(tmp_bucket, 0, 256);
            }
            pmemobj_persist(pool, bucket_p, sizeof(Bucket<T>) * bucket_number__);
            return 0;
        };
        PMAllocator::Allocate((void **)bkt, CACHELINESIZE, sizeof(Bucket<T>) * bucket_number, init_callback,
                            reinterpret_cast<void *>(&bucket_number));
    }

    static void New(PMEMoid *bkt, uint32_t bucket_number_, size_t arg) {
        auto init_callback = [](PMEMobjpool *pool, void *ptr, void *arg) {
            auto bucket_p = reinterpret_cast<Bucket *>(ptr);
            uint32_t bucket_number__ = *(reinterpret_cast<uint32_t *>(arg));
            for(int i = 0; i < bucket_number__; i ++)
            {
                auto tmp_bucket = bucket_p + i;
                memset(tmp_bucket, 0, 256);
            }
            pmemobj_persist(pool, bucket_p, sizeof(Bucket<T>) * bucket_number__);
            return 0;
        };
        PMAllocator::Allocate(bkt, CACHELINESIZE, sizeof(Bucket<T>) * bucket_number_, init_callback,
                            reinterpret_cast<void *>(&bucket_number_));
    }
    /* lock operations */
    // get the lock: try get the lock util locked
    inline void get_lock()
    {
        uint32_t new_version = 0;
        uint32_t old_version = 0;
        do
        {
            while(true)
            {
                old_version = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
                if(!(old_version & lock_pos))
                {
                    old_version &= lock_mask;
                    break;
                }
            }
            new_version = old_version | lock_pos;
        } while(!CAS(&version_lock, &old_version, new_version));
    }

    // try get the lock: try once, return the result
    inline bool try_get_lock()
    {
        uint32_t version = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
        if(version & lock_pos) // if lock pos is set, return false.
            return false;
        auto old_version = version & lock_mask;
        auto new_version = version | lock_pos;
        return CAS(&version_lock, &old_version, new_version);
    }

    // release the lock: set lock pos and version num ++
    inline void release_lock()
    {
        uint32_t v = version_lock;
        __atomic_store_n(&version_lock, v + 1 - lock_pos, __ATOMIC_RELEASE);
    }

    // test lock set: test if the lock pos is set, return lock pos
    inline bool test_lock_set(uint32_t &version)
    {
        version = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
        return (version & lock_pos) != 0;
    }

    inline bool test_lock_version_change(uint32_t old_version)
    {
        auto value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
        return (old_version != value);
    }

    inline int find_empty_pos()
    {
        auto mask = ~(bitmap);
        return __builtin_ctz(mask);
    }

    inline void set_bitmap_and_fp(T key, int empty_pos)
    {
        bitmap |= 1 << empty_pos;
        fp[empty_pos] = key & fp_mask;
    }

    bool verify_unique(T key_)
    {
        int mask = 0;
        uint8_t meta_hash = key_ & fp_mask;
        int i;
        #ifdef USING_SSE
            SSE_CMP8(fp, key_ & fp_mask);
            mask &= (1 << BUCKET_CAPACITY) - 1;
        #else
            for(i = 0; i < BUCKET_CAPACITY;  i ++)
            {
                if(fp[i] == meta_hash)
                {
                    mask |= (1 << i);
                }
            }
        #endif
        mask = mask & bitmap;
        if(mask == 0)
            return false;
        // OPTIMSE: loop unrolling
        #ifndef LOOP_UNROLLING
        for(i = 0; i < BUCKET_CAPACITY; i += 1){
            // int valid = bitmap & (1 << i);
            if(CHECK_BIT(mask, i) && kv[i].key == key_)
            {
                return true;
            }
        }
        #else
        for(i = 0; i < BUCKET_CAPACITY; i += 4){
            if(CHECK_BIT(mask, i) && kv[i].key == key_)
            {
                return true;
            }

            if(CHECK_BIT(mask, i + 1) && kv[i + 1].key == key_)
            {
                return true;
            }

            if(CHECK_BIT(mask, i + 2) && kv[i + 2].key == key_)
            {
                return true;
            }

            if(CHECK_BIT(mask, i + 3) && kv[i + 3].key == key_)
            {
                return true;
            }
        }
        
        if(CHECK_BIT(mask, 8) && kv[8].key == key_)
        {
            return true;
        }

        if(CHECK_BIT(mask, 9) && kv[9].key == key_)
        {
            return true;
        }

        if(CHECK_BIT(mask, 10) && kv[10].key == key_)
        {
            return true;
        }
        #endif
        return false;
    }

    /*
        return: 
            -1: insert failed
            0: redundant insertion
            1: redundant insertion
    */
    int insert_kv_to_bucket(T key_, Value_t value_, int cell_index, int iguide){
        // verify first
        bool ret = verify_unique(key_);
        if(ret)
            return 0;

        if(used >= BUCKET_CAPACITY)
            return -1;
        // append kv pair
        int empty_pos = find_empty_pos();
        kv[empty_pos].key   = key_;
        kv[empty_pos].value = value_;
        PMAllocator::Persist(&kv[empty_pos], sizeof(kv[empty_pos]));
        // append metadata
        k[empty_pos] = cell_index;
        init_guide[empty_pos] = iguide;
        set_bitmap_and_fp(key_, empty_pos);
        used++;
        PMAllocator::Persist(&bitmap, sizeof(bitmap) + sizeof(used) + 
                        sizeof(k) + sizeof(init_guide));
        // debug
        bucket_items++;
        ins_cnt++;
        return 1;
    }

    // collect kvs from bucket indexed by guide
    inline void collect_kv(int cell_k, _Pair<T> *kvs, uint8_t local_depth, 
                        int *init_guides, int iguide, int &size, uint32_t initial_cell){   
        for(int i = 0; i < BUCKET_CAPACITY; i++)
        {
            int valid = bitmap & (1 << i);
            if(valid == 0)
                continue;
            if(k[i] % initial_cell == cell_k % initial_cell
                && iguide == init_guide[i])
            {
                init_guides[size] = iguide;
                kvs[size].key = kv[i].key;
                kvs[size].value = kv[i].value;
                // memcpy(kvs[size].key, key[i], (KEY_LEN)*sizeof(char));
                // memcpy(kvs[size].value, value[i], (VAL_LEN)*sizeof(char));
                size++;
                
            }
        }
    }

    bool query_key_in_bucket(T key_, Value_t value_ = NULL){
#ifdef USING_SIMD
        const __m256i item = _mm256_set1_epi32(*(int*)key);
        __m256i *keys_p = (__m256i *)(bucket[d].key);
        int matched = 0;

        __m256i a_comp = _mm256_cmpeq_epi32(item, keys_p[0]);
        matched = _mm256_movemask_ps((__m256)a_comp);
        

        if(matched != 0){
            int matched_lowbit = matched & (-matched);
            int matched_index = __tzcnt_u32((uint32_t)matched_lowbit);
            if(matched_index < bucket[d].used){
                if(value != NULL){
                    memcpy(value, bucket[d].value[matched_index], VAL_LEN*sizeof(char));
                }
                return true;
            }
            return false;
        }
        return false;

#else   
        int mask = 0;
        uint8_t meta_hash = key_ & fp_mask;
        int i;
        #ifdef USING_SSE
            SSE_CMP8(fp, key_ & fp_mask);
            mask &= (1 << BUCKET_CAPACITY) - 1;
        #else
            for(i = 0; i < BUCKET_CAPACITY;  i ++)
            {
                if(fp[i] == meta_hash)
                {
                    mask |= (1 << i);
                }
            }
        #endif
        mask &= bitmap;
        if(mask == 0)
            return false;
        #ifndef LOOP_UNROLLING
        for(i = 0; i < BUCKET_CAPACITY; i += 1){
            // int valid = bitmap & (1 << i);
            if(CHECK_BIT(mask, i) && kv[i].key == key_)
            {
                if(value_ != NULL)
                    value_ = kv[i].value;
                return true;
            }
        }
        #else
        for(i = 0; i < BUCKET_CAPACITY; i += 4){
            if(CHECK_BIT(mask, i) && kv[i].key == key_)
            {
                return true;
            }

            if(CHECK_BIT(mask, i + 1) && kv[i + 1].key == key_)
            {
                return true;
            }

            if(CHECK_BIT(mask, i + 2) && kv[i + 2].key == key_)
            {
                return true;
            }

            if(CHECK_BIT(mask, i + 3) && kv[i + 3].key == key_)
            {
                return true;
            }
        }
        
        if(CHECK_BIT(mask, 8) && kv[8].key == key_)
        {
            return true;
        }

        if(CHECK_BIT(mask, 9) && kv[9].key == key_)
        {
            return true;
        }

        if(CHECK_BIT(mask, 10) && kv[10].key == key_)
        {
            return true;
        }
        #endif
        
        return false;
#endif
    }

    bool remove_record_from_bucket(T key_)
    {
        int mask = 0;
        uint8_t meta_hash = key_ & fp_mask;
        int i;
        #ifdef USING_SSE
            SSE_CMP8(fp, key_ & fp_mask);
            mask &= (1 << BUCKET_CAPACITY) - 1;
        #else
            for(i = 0; i < BUCKET_CAPACITY;  i ++)
            {
                if(fp[i] == meta_hash)
                {
                    mask |= (1 << i);
                }
            }
        #endif
        mask &= bitmap;
        if(mask == 0)
            return false;
        int index = 0;
        for(i = 0; i < BUCKET_CAPACITY; i++)
        {
            index = __builtin_ctz(mask);
            if(CHECK_BIT(mask, index) && key_ == kv[index].key)
            {
                used --;
                bitmap &= ~(1 << index);
                PMAllocator::Persist(&bitmap, sizeof(bitmap) + sizeof(used));
                return true;
            }
            mask &= ~(1 << index);
            if(mask == 0)
                return false;
        }
        return false;
    }
};


struct Cell{
    uint32_t version_lock;
    uint8_t guide_offset[16];
    uint32_t cell_bucket[16];
    #ifndef DIRECTORY
        uint8_t local_depth;
    #endif
    /* lock operations, same as bucekt lock operations*/
    // get the lock: try get the lock util locked
    inline void get_lock()
    {
        uint32_t new_version = 0;
        uint32_t old_version = 0;
        do
        {
            while(true)
            {
                old_version = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
                if(!(old_version & lock_pos))
                {
                    old_version &= lock_mask;
                    break;
                }
            }
            new_version = old_version | lock_pos;
        } while(!CAS(&version_lock, &old_version, new_version));
    }

    // try get the lock: try once, return the result
    inline bool try_get_lock()
    {
        uint32_t version = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
        if(version & lock_pos) // if lock pos is set, return false.
            return false;
        auto old_version = version & lock_mask;
        auto new_version = version | lock_pos;
        return CAS(&version_lock, &old_version, new_version);
    }

    // release the lock: set lock pos and version num ++
    inline void release_lock()
    {
        uint32_t v = version_lock;
        __atomic_store_n(&version_lock, v + 1 - lock_pos, __ATOMIC_RELEASE);
    }

    // test lock set: test if the lock pos is set, return lock pos
    inline bool test_lock_set(uint32_t &version)
    {
        version = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
        return (version & lock_pos) != 0;
    }

};

#ifdef DIRECTORY
struct Directory{
    uint16_t local_depth;
    Cell *cell_ptr;
};
Directory *dir;
#else
Cell *cell;
#endif

template<class T>
struct Layer
{
    typedef Bucket<T> * bucket_p;
    uint32_t global_depth;
    bucket_p bp[0];

    static void New(PMEMoid *layer, uint32_t gd, size_t arg)
    {
        auto init_callback = [](PMEMobjpool *pool, void *layer_p, void *arg)
        {
            auto layer_p_ = reinterpret_cast<Layer<T> *>(layer_p);
            uint8_t gd = *(reinterpret_cast<uint8_t *>(arg));
            uint32_t layer_num_ = pow(2, gd);
            layer_p_->global_depth = gd;
            pmemobj_persist(pool, layer_p, sizeof(Layer<T>) + sizeof(bucket_p) * layer_num_);
            return 0;
        };
        uint32_t layer_num = pow(2, gd);
        PMAllocator::Allocate(layer, CACHELINESIZE, sizeof(Layer<T>) + sizeof(bucket_p) * layer_num,
                            init_callback, reinterpret_cast<void *>(&gd));
    }

    static void New(Layer<T> **layer, uint32_t gd, size_t arg)
    {
        auto init_callback = [](PMEMobjpool *pool, void *layer_p, void *arg)
        {
            auto layer_p_ = reinterpret_cast<Layer<T> *>(layer_p);
            uint8_t gd = *(reinterpret_cast<uint8_t *>(arg));
            uint32_t layer_num_ = pow(2, gd);
            layer_p_->global_depth = gd;
            pmemobj_persist(pool, layer_p, sizeof(Layer<T>) + sizeof(bucket_p) * layer_num_);
            return 0;
        };
        uint32_t layer_num = pow(2, gd);
        PMAllocator::Allocate(layer, CACHELINESIZE, sizeof(Layer<T>) + sizeof(bucket_p) * layer_num,
                            init_callback, reinterpret_cast<void *>(&gd));
    }
};

template <class T>
class EEPH : public Hash<T>{
public:
    uint32_t init_bucket_number;   // initial bucket number
    uint32_t init_cell_number;     // initial cell number
    uint8_t cell_hash;             // number of hash functions in each layer
#ifdef USING_SIMD
    __attribute__((aligned(32))) Bucket<T> *bucket;
#else
    // Bucket<T> *bucket;
#endif
    uint32_t seed_hash_to_cell[MAX_LAYER];
    uint32_t seed_hash_to_bucket[M];
    uint32_t seed_hash_to_guide;
    PMEMoid layers_;
    Layer<T> *layers;
    uint32_t lock; // lock bit + version number
    // interface
    inline int Insert(T key, Value_t value);
    int Insert(T key, Value_t value, bool is_in_epoch);
    inline bool Delete(T);
    bool Delete(T, bool);
    inline Value_t Get(T);
    Value_t Get(T key, bool is_in_epoch);
    void recovery();
    bool find_kv(T key);
    inline int Insert_seq(T key, Value_t value, bool is_in_epoch, uint64_t cur_cnt);
    inline int Insert_rnd(T key, Value_t value, bool is_in_epoch);
    inline Value_t Get_seq(T key, bool is_in_epoch, uint64_t cur_cnt);
    inline Value_t Get_rnd(T key, bool is_in_epoch);
    void Recovery()
    {

    }
    void getNumber()
    {
        int items_num = statistics_bucket_items();
        std::cout << "the size of bucket is " << sizeof(Bucket<T>) << std::endl;
        std::cout << "the number of bucket items is " << items_num << std::endl;
        for(int i  = 0; i < 1 << layers->global_depth; i++)
            std::cout << "level" << i << " items is " << bucket_items_per_level[i] << std::endl;
        std::cout << "the load factor is " << (double)items_num / capacity << std::endl;
        
        std::cout << "===================== Insert Time =======================" << std::endl;
        std::cout << "hash time: " << (double)ihash_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "read cell time: "<< (double)iread_cell_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "cal time: "<< (double)iread_cal_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "calculate bucket time: "<< (double)icalculate_bucket_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "insert time: " << (double)insert_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "threshold time: " << (double)threshold_time / CLOCKS_PER_SEC<< std::endl;
        std::cout << "move time: " << (double)move_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "extend time: " << (double)extend_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "  double bucket time: " << (double)double_bucket_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "  copy bucket time: " << (double)copy_bucket_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "  double cell time: " << (double)double_cell_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "collect time: " << (double)collect_time / CLOCKS_PER_SEC<< std::endl;
        std::cout << "===================== Insert Time =======================" << std::endl;

        std::cout << "===================== Search Time =======================" << std::endl;
        std::cout << "cell time: " << (double)cell_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "init guide time: "<< (double)init_guide_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "bucket time: "<< (double)bucket_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "search time: " << (double)search_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "hash time: " << (double)hash_time / CLOCKS_PER_SEC<< std::endl;
        std::cout << "guide 1 time: " << (double)guide_1_time / CLOCKS_PER_SEC<< std::endl;
        std::cout << "level time: "<< (double) level_time / CLOCKS_PER_SEC<< std::endl;
        std::cout << "level 1 time: " << (double)level_1_time / CLOCKS_PER_SEC<< std::endl;
        std::cout << "read bucket time: " << (double)read_bucket_time / CLOCKS_PER_SEC<< std::endl;
        std::cout << "===================== Search Time =======================" << std::endl;

        std::cout << "===================== Delete Time =======================" << std::endl;
        std::cout << "cell time: " << (double)dcell_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "init guide time: "<< (double)dinit_guide_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "bucket time: "<< (double)dbucket_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "search time: " << (double)dsearch_time / CLOCKS_PER_SEC << std::endl;
        std::cout << "layer time: " << (double)dohter_time / CLOCKS_PER_SEC<< std::endl;
        std::cout << "===================== Delete Time =======================" << std::endl;

        std::cout << "===================== Count =======================" << std::endl;
        std::cout << "cnt: " << cnt << std::endl;
        std::cout << "insert cnt: " << ins_cnt << std::endl;
        std::cout << "inert directly cnt: " << insert_dire_cnt << std::endl;
        std::cout << "threshold cnt: " << threshold_cnt << std::endl;
        std::cout << "move cnt: " << move_cnt << std::endl;
        std::cout << "extend cnt: " << extend_cnt << std::endl;
        std::cout << "rehash cnt: " << rehash_cnt << std::endl;
        std::cout << "rehash muti cnt: " << rehash_muti_cnt << std::endl;
        std::cout << "===================== Count =======================" << std::endl;
    }

    inline bool get_cells_version_lock()
    {
        uint32_t v = __atomic_load_n(&lock, __ATOMIC_ACQUIRE);
        uint32_t old_version = v & lock_mask;
        uint32_t new_version = (old_version + 1) & lock_mask;
        return CAS(&lock, &old_version, new_version);
    }

    inline void release_cells_version_lock()
    {
        SUB(&lock, 1);
    }

    inline void lock_cells()
    {
        uint32_t v = __atomic_load_n(&lock, __ATOMIC_ACQUIRE);
        uint32_t old_version = v & lock_mask;
        uint32_t new_version = old_version | lock_pos;

        // wait for lock == old_lock
        while(!CAS(&lock, &old_version, new_version))
        {
            // do nothing
            old_version = old_version & lock_mask;
            new_version = new_version | lock_pos;
        }

        // wait until release all version lock
        v = __atomic_load_n(&lock, __ATOMIC_ACQUIRE);
        while(v & lock_mask)
            v = __atomic_load_n(&lock, __ATOMIC_ACQUIRE);
    }

    inline void unlock_cells()
    {
        __atomic_store_n(&lock, 0, __ATOMIC_RELEASE);
    }

private:
#ifdef PMEM
    static int ld_callback(PMEMobjpool *pool, void *ptr, void *arg)
    {
       return 0; 
    }
#endif
    void initialize_hash_functions(){
        std::set<int> seeds;
        uint32_t seed = rand()%MAX_PRIME32;
        seed_hash_to_guide = seed;
        seeds.insert(seed);
        for(int i = 0; i < 1; ++i){
            seed = rand()%MAX_PRIME32;
            while(seeds.find(seed) != seeds.end())
                seed = rand()%MAX_PRIME32;
            seed_hash_to_cell[i] = seed;
            seeds.insert(seed);
        }
        for(int i = 0; i < M; ++i){
            seed = rand()%MAX_PRIME32;
            while(seeds.find(seed) != seeds.end())
                seed = rand()%MAX_PRIME32;
            seed_hash_to_bucket[i] = seed;
            seeds.insert(seed);
        }
    }

    // calculate cell index
    uint32_t calculate_cell(const void* key){
        // cost: calculate 1 hash
        #ifdef XXHASH
        return XXH32(key, KEY_LEN, seed_hash_to_cell[SEED]);
        #else
        return MurmurHash3_x86_32(key, KEY_LEN, seed_hash_to_cell[SEED]);
        #endif
    }

    // calculate virtual bucket index
    int calculate_guide(const void* key){
        // cost: calculate 1 hash
        return ((*(uint64_t*)key) & key_mask) % cell_hash;
        int ret;
        if(flag)
        {
            timeval start, stop;
            gettimeofday(&start, NULL);
            #ifdef XXHASH
            tmp = XXH32(key, KEY_LEN, seed_hash_to_guide) % cell_hash;
            #else
            ret = MurmurHash3_x86_32(key, KEY_LEN, seed_hash_to_guide) % cell_hash;
            #endif
            gettimeofday(&stop, NULL);
            guide_1_time += (double)(stop.tv_usec - start.tv_usec) +
                            (double)(stop.tv_sec - start.tv_sec) * 1000000;
        }else
            #ifdef XXHASH
            tmp = XXH32(key, KEY_LEN, seed_hash_to_guide) % cell_hash;
            #else
            ret = MurmurHash3_x86_32(key, KEY_LEN, seed_hash_to_guide) % cell_hash;
            #endif
        return ret;
    }

    int calculate_level(int k, int *level, uint8_t local_depth, uint8_t global_depth)
    {
        int ret;
        int cell_numbers = init_cell_number << global_depth;
        int initial_cell = init_cell_number;
        int delta_depth = global_depth - local_depth;
        if(local_depth == global_depth)
        {
            *level = k / initial_cell;
            return 0;
        }
        else if(local_depth < global_depth)
        {
            if(k < cell_numbers >> delta_depth)
            {
                *level = k / initial_cell;
                return 0;
            }else
            {
                *level = (k % (cell_numbers >> delta_depth)) / initial_cell;
                return 0;
            }
        }
        return -1;
    }

    /*
    discribe: calculate physical bucket index
    param:
        guide      : virtual bucket index
        k          : cell index
        local_depth: local depth
    cost: calculate hash once
    */
    int calculate_bucket(int guide, int k, uint8_t local_depth){
        int level, cell_ID = k % init_cell_number, ret;
        int bkt_num_per_layer;
        int page_size;
        timeval start, stop;

        bkt_num_per_layer = init_bucket_number;
        page_size = bkt_num_per_layer / cell_hash;
        if(flag)
        {
            gettimeofday(&start, NULL);
            #ifdef XXHASH
            ret = XXH32((const void *)&cell_ID, sizeof(int), 
                        seed_hash_to_bucket[guide]) % page_size + guide * page_size;
            #else
            ret = MurmurHash3_x86_32((const void *)&cell_ID, sizeof(int),
                    seed_hash_to_bucket[guide]) % page_size + guide * page_size;
            #endif
            gettimeofday(&stop, NULL);
            hash_time += (double)(stop.tv_usec - start.tv_usec) +
                            (double) (stop.tv_sec - start.tv_sec) * 1000000;
        }else
            #ifdef XXHASH
            ret = XXH32((const void *)&cell_ID, sizeof(int),
                seed_hash_to_bucket[guide]) % page_size + guide * page_size;
            #else
            ret = MurmurHash3_x86_32((const void *)&cell_ID, sizeof(int),
                seed_hash_to_bucket[guide]) % page_size + guide * page_size;
            #endif
        return ret;
    }

    inline bool is_movable(int k, uint32_t d, uint32_t iguide, uint8_t local_depth, 
                        uint8_t global_depth, uint32_t initial_cell, uint32_t *max_pos, uint32_t *max_guide)
    {
        int level, size = 1, max_empty = 0;
        calculate_level(k, &level, local_depth, global_depth);
        Bucket<T> *bucket = layers->bp[level];
        initial_cell = initial_cell << local_depth;
        for(int guide = 0; guide < cell_hash; ++ guide){
            int tmp_d = calculate_bucket(guide, k, local_depth);
            // int tmp_d = cell_bucket_[guide];
            if(tmp_d == d)
            {
                int source_pre_version = __atomic_load_n(&bucket[d].version_lock, __ATOMIC_ACQUIRE);
                for(int i = 0; i < BUCKET_CAPACITY; i++)
                {
                    int valid = bucket[d].bitmap & (1 << i);
                    if(valid == 0)
                        continue;
                    if(bucket[d].k[i] % initial_cell == k % initial_cell
                        && iguide == bucket[d].init_guide[i])
                    {
                        size++;
                    }
                }
            }else
            {
                int rest = BUCKET_CAPACITY - bucket[tmp_d].used;
                if(rest > max_empty)
                {
                    max_empty = rest;
                    *max_pos = tmp_d; 
                    *max_guide = guide;
                }
            }
        }
        if(max_empty >= size)
        {
            // #ifdef DIRECTORY
            //     cell_->guide_offset[init_guide] += (max_guide > guide)? 
            //         (max_guide - guide) : (cell_hash - (guide - max_guide));
            // #else
            //     cell[k].guide_offset[init_guide] += (max_guide > guide)? 
            //         (max_guide - guide) : (cell_hash - (guide - max_guide));
            // #endif
            return true;
        }
        else
            return false;
    }

    bool rehash_a_cell(int k, uint32_t *cell_bucket_, uint8_t *guide_offset_, 
            uint8_t local_depth, uint8_t global_depth, uint32_t cell_number, uint8_t level){   
    REHASH:    
        bool ret = true;
        int guide;
        Bucket<T> *bucket = layers->bp[level];
        int initial_cell = init_cell_number;
        int initial_bucket = init_bucket_number;
        for(guide = 0; guide < cell_hash; ++guide){
            int guide_ = (guide + guide_offset_[guide]) % cell_hash;
            int d = cell_bucket_[guide_];
            bucket[d].get_lock();
            // int d = calculate_bucket(guide, k, local_depth);
            for(int j = 0; j < BUCKET_CAPACITY; ++j){
                int valid = bucket[d].bitmap & (1 << j);
                if(valid == 0)
                    continue;
                if(bucket[d].k[j] % initial_cell == k % initial_cell){
                    int cell_index = calculate_cell(&bucket[d].kv[j].key) % cell_number;
                    // int cell_index = bucket[d].k[j];
                    int cell_number_1 = cell_number >> (global_depth - local_depth);
                    int cell_number_2 = (cell_number * 2) >> (global_depth - local_depth);
                    int pos = cell_index % cell_number_2;
                    if(pos == pos % cell_number_1)
                        continue;
                    ret &= insert(bucket[d].kv[j].key, bucket[d].kv[j].value);
                    if(ret == false)
                    {
                        return false;
                    }
                    bucket[d].bitmap &= ~(1 << j);
                    bucket[d].used --;
                    bucket_items --;
                    // flush metadata
                    PMAllocator::Persist(&bucket[d].bitmap, sizeof(bucket[d].bitmap) + sizeof(bucket[d].used));
                }
            }
            bucket[d].release_lock();
        }
        return true;
    }

    bool fine_grained_circular_move(int cell_k, uint8_t local_depth, int iguide, 
                                uint32_t initial_cell, Bucket<T> *source_bucket, Bucket<T> *target_bucket, T key, Value_t value)
    {
        int ret = 0;
        initial_cell = initial_cell << local_depth;
        for(int i = 0; i < BUCKET_CAPACITY; i++)
        {
            int valid = source_bucket->bitmap & (1 << i);
            if(valid == 0)
                continue;
            if(source_bucket->k[i] % initial_cell == cell_k % initial_cell
                && iguide == source_bucket->init_guide[i])
            {
                // set valid bitmap
                source_bucket->bitmap &= ~(1 << i);
                source_bucket->used --;
                bucket_items --;

                // flush metadata
                PMAllocator::Persist(&source_bucket->bitmap, sizeof(source_bucket->bitmap) + sizeof(source_bucket->used));

                // insert into target bucket
                ret = target_bucket->insert_kv_to_bucket(source_bucket->kv[i].key, source_bucket->kv[i].value, cell_k, iguide);
                if(ret <= 0)
                    return false;
            }
        }
        if(target_bucket->insert_kv_to_bucket(key, value, cell_k, iguide) <= 0)
            return false;
        return true;
    }

    int before_dynamic_extend(Bucket<T> *source_bucket, uint8_t global_depth, 
                                uint32_t cell_number, uint8_t level)
    {
        for(int i  = 0; i < BUCKET_CAPACITY; i ++)
        {
            uint32_t k = source_bucket->k[i];
            uint8_t local_depth = cell[k].local_depth;
            if(local_depth == global_depth)
                continue;
            // int cell_index = source_bucket.k[i];
            int cell_index = calculate_cell(&source_bucket->kv[i].key) % cell_number;
            int cell_number_1 = cell_number >> (global_depth - local_depth);
            int cell_number_2 = (cell_number * 2) >> (global_depth - local_depth);
            int pos = cell_index % cell_number_2;
            if(pos != pos % cell_number_1)
            {
                uint32_t *cell_bucket = cell[cell_index].cell_bucket;
                uint8_t *guide_offset = cell[cell_index].guide_offset;
                
                // rehash a cell
                int diff_depth = global_depth - local_depth;
                int diff_cell = cell_number >> diff_depth;
                int first_offset = k % diff_cell;
                lock_cells();
                uint8_t new_global_depth = layers->global_depth;
                uint8_t new_local_depth = cell[cell_index].local_depth;
                if(new_global_depth != global_depth || new_local_depth != local_depth)
                {
                    unlock_cells();
                    return 0;
                }
                for(int i = 0; i < 1 << diff_depth; i ++)
                    cell[first_offset + i * diff_cell].local_depth ++;
                unlock_cells();
                rehash_a_cell(cell_index, cell_bucket, guide_offset, local_depth, 
                            global_depth, cell_number, level);
                return 1;                
            }
        }
        return -1;
    }

    bool dynamic_extension(uint8_t global_depth)
    {
        timeval start, stop;
        gettimeofday(&start, NULL);
        #ifdef INSERT_CLOCK
            timeval start, stop;
            gettimeofday(&start, NULL);
        #endif
        std::cout << "load factory: " << (double)bucket_items / capacity << std::endl; 
        std::cout << "dynamiccccccccccc" << std::endl;
        bool    ret            = true;
        uint32_t pow_gd_= pow(2, global_depth);
        uint32_t cell_number   = init_cell_number << global_depth;
        uint32_t pow_gd = pow(2, global_depth + 1);
        PMEMoid oid;

        // double cells
         #ifdef DIRECTORY
            Directory *old_dir = dir;

            dir = (Directory *)malloc(sizeof(Directory) * cell_number * 2);
            // PMAllocator::Allocate((void **)&cell, CACHELINESIZE, sizeof(size_t) * cell_number, ld_callback,
            //                 NULL);
            memcpy(dir, old_dir, sizeof(Directory) * cell_number);
            memcpy(dir + cell_number, old_dir, sizeof(Directory) * cell_number);
            free(old_dir);
        #else
            Cell *old_cell = cell;
            cell = (Cell *)malloc(sizeof(Cell) * cell_number * 2);
            // PMAllocator::Allocate((void **)&cell, CACHELINESIZE, sizeof(size_t) * cell_number, ld_callback,
            //                 NULL);
            memcpy(cell, old_cell, sizeof(Cell) * cell_number);
            memcpy(cell + cell_number, old_cell, sizeof(Cell) * cell_number);
            free(old_cell);
        #endif
        #ifdef INSERT_CLOCK
            gettimeofday(&stop, NULL);
            double_cell_time += (double)(stop.tv_usec - start.tv_usec) + 
                                (double)(stop.tv_sec - start.tv_sec) * 1000000;
        #endif

        // double layers and bucket
        Layer<T>::New(&layers_, global_depth + 1, 0);
        Layer<T> *tmp_layers = (Layer<T> *)pmemobj_direct(layers_);
        for(int i = pow_gd_; i < pow_gd; i ++)
        {
            Bucket<T>::New(&oid, init_bucket_number, 0);
            tmp_layers->bp[i] = (Bucket<T> *)pmemobj_direct(oid);
        }
        #ifdef INSERT_CLOCK
            gettimeofday(&stop, NULL);
            double_bucket_time += (double)(stop.tv_usec - start.tv_usec) +
                                (double)(stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        // copy old buckets, low performance, to be optimized
        for(int i = 0; i < pow_gd_; i ++)
        {
            tmp_layers->bp[i] = layers->bp[i];
        }
        PMAllocator::Free_PM(layers_);
        layers = tmp_layers;
        //persist
        PMAllocator::Persist(layers, sizeof(Layer<T>) + sizeof(Bucket<T> *) * pow_gd);
        capacity *= 2;
        #ifdef INSERT_CLOCK
            gettimeofday(&stop, NULL);
            copy_bucket_time += (double)(stop.tv_usec - start.tv_usec) + 
                                (double)(stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        gettimeofday(&stop, NULL);
        std::cout << "extend time: " << (double)(stop.tv_usec - start.tv_usec) +
                                        (stop.tv_sec - start.tv_sec) * 1000000 << std::endl;
        return ret;
    }

    bool rehash_a_cell_1(_Pair<T> *kvs, int size, int k, uint8_t local_depth, uint32_t global_depth, 
                        uint32_t cell_number)
    {
        bool ret = true;
        int i;
        int diff_depth = global_depth - local_depth;
        int diff_cell = cell_number >> diff_depth;
        int first_offset = k % diff_cell;
        #ifdef DIRECTORY
            Cell *cell_ = (Cell *)malloc(sizeof(Cell));
            memcpy(cell_, dir[k].cell_ptr, sizeof(Cell));
            for(i = 0; i < 1 << diff_depth; i ++)
            {
                dir[first_offset + i * diff_cell].local_depth ++;
                if(i % 2 == 1)
                    dir[first_offset + i * diff_cell].cell_ptr = cell_;
        }
        #else
            for(i = 0; i < 1 << diff_depth; i ++)
                cell[first_offset + i * diff_cell].local_depth ++;
        #endif
        for(i = 0; i < size; i++)
            ret &= insert(kvs[i].key, kvs[i].value);
        return ret;
    }

public:
    EEPH(void);
    EEPH(uint32_t bucket_number_, uint32_t cell_number_, uint8_t cell_hash_);

    ~EEPH(void);

    int insert_post(T key, Value_t value, uint32_t k, uint32_t d, uint8_t init_guide, 
            uint8_t guide, uint8_t local_depth)
    {
        int      i;
        /*--------------------------------stage 1-------------------------------*/
        uint8_t  global_depth  = layers->global_depth;
        uint32_t cell_number   = init_cell_number << global_depth;
        uint32_t bucket_number = init_bucket_number << global_depth;
        bool ret = true;

        int level;
        Cell cell_ = cell[k];
        uint32_t *cell_bucket_ = cell_.cell_bucket;
        uint8_t  *guide_offset_ = cell_.guide_offset;
        // int guide = (init_guide + guide_offset_[init_guide]) % cell_hash;
        calculate_level(k, &level, local_depth, global_depth);
        Bucket<T> *bucket = layers->bp[level];
        bucket[d].get_lock();
        if(bucket[d].insert_kv_to_bucket(key, value, k, init_guide) == 1)
        {
            bucket[d].release_lock();
            return true;
        }
        bucket[d].release_lock();
        uint32_t max_pos = 0;
        uint32_t max_guide = 0;

        if(local_depth == global_depth)
        {
            if(is_movable(k, d, init_guide, local_depth, global_depth, 
                    init_cell_number, &max_pos, &max_guide))
            {
                // circular move
                cell[k].guide_offset[init_guide] += (max_guide > guide)? 
                    (max_guide - guide) : (cell_hash - (guide - max_guide));
                bucket[max_pos].get_lock();
                bucket[d].get_lock();
                ret &= fine_grained_circular_move(k, local_depth, init_guide, init_cell_number,
                                                bucket + d, bucket + max_pos, key, value);
                bucket[d].release_lock();
                bucket[max_pos].release_lock();
                return true;
            }else
            {
                // dynamic extend
                if(before_dynamic_extend(bucket + d, global_depth, cell_number, level))
                {
                    ret &= insert(key, value);// goto
                    return ret;
                }
                std::cout << "load factory: " << (double)bucket_items / capacity << std::endl; 
                std::cout << "dynamiccccccccccc" << std::endl;
                lock_cells();
                ret &= dynamic_extension(global_depth);
                unlock_cells();
                global_depth += 1;
                cell_number *= 2;
                std::cout << "extend time: " <<  extend_time/CLOCKS_PER_SEC << std::endl;
                ret &= rehash_a_cell(k, cell_bucket_, guide_offset_, 
                        local_depth, global_depth, cell_number, level);
                ret &= insert(key, value);// goto
                return ret;
            }
        }else if(local_depth < global_depth)
        {
            ret &= rehash_a_cell(k, cell_bucket_, guide_offset_, 
                    local_depth, global_depth, cell_number, level);
            ret &= insert(key, value);// goto
            return ret;
            
        }
        return false;
    }

    bool insert(T key, Value_t value){
         int      i;
        _Pair<T> kvs[M * BUCKET_CAPACITY + 5];
        int      size = 0;
        cnt++;
        // if(key == 646091557254529188)
        //     printf("herrrrrrrr\n");
        /*--------------------------------stage 1-------------------------------*/
        #ifdef INSERT_CLOCK
            timeval start, stop;
            gettimeofday(&start, NULL);
        #endif
        // k: cell index
        int  init_guide = calculate_guide(&key);
        int  k = calculate_cell(&key);
    REDO:
        uint8_t  global_depth  = layers->global_depth;
        uint32_t cell_number   = init_cell_number << global_depth;
        uint32_t bucket_number = init_bucket_number << global_depth;
        bool ret = true;
        int result = 0;
        k %= cell_number;
        
        #ifdef INSERT_CLOCK
            gettimeofday(&stop, NULL);
            ihash_time += (double)(stop.tv_usec - start.tv_usec) +
                                (double)(stop.tv_sec - start.tv_sec) * 1000000;
            
            gettimeofday(&start, NULL);
        #endif
        #ifdef DIRECTORY
            Cell *cell_ = dir[k].cell_ptr; 
            uint8_t local_depth = dir[k].local_depth;
            uint32_t *cell_bucket_ = cell_->cell_bucket;
            uint8_t  *guide_offset_ = cell_->guide_offset;
        #else
    REINSERT:
            Cell cell_ = cell[k];
            uint8_t local_depth = cell_.local_depth;
            uint32_t *cell_bucket_ = cell_.cell_bucket;
            uint8_t  *guide_offset_ = cell_.guide_offset;
        #endif
        #ifdef INSERT_CLOCK
            gettimeofday(&stop, NULL);
            iread_cell_time += (double)(stop.tv_usec - start.tv_usec) +
                                (double)(stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        int guide = (init_guide + guide_offset_[init_guide]) % cell_hash;
        int d = cell_bucket_[guide];
        // int d = calculate_bucket(guide, k, cell_.local_depth);
        #ifdef INSERT_CLOCK
            gettimeofday(&stop, NULL);
            iread_cal_time += (double)(stop.tv_usec - start.tv_usec) +
                                (double)(stop.tv_sec - start.tv_sec) * 1000000;
            
            gettimeofday(&start, NULL);
        #endif
        int level;
        calculate_level(k, &level, local_depth, global_depth);
        #ifdef INSERT_CLOCK
            gettimeofday(&stop, NULL);
            icalculate_bucket_time += (double)(stop.tv_usec - start.tv_usec) +
                                (double)(stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        // uint8_t local_depth = cell_.local_depth;
        Bucket<T> *bucket = layers->bp[level];
        Bucket<T> &bucket_ = bucket[d];
        #ifdef THRESHOLD
            if(local_depth < global_depth)
            {
                float bucket_capacity_ratio = (float) bucket_.used / BUCKET_CAPACITY; // bucket capacity
                float table_capacity_ratio = (float) bucket_items / capacity;
                if(bucket_capacity_ratio > table_capacity_ratio)
                {
                    // rehash items in this cell
                    threshold_cnt ++;
                    ret &= rehash_a_cell(k, cell_bucket_, guide_offset_, 
                                    local_depth, global_depth, cell_number, level);
                    if(ret == false)
                        return false;
                    ret &= insert(key, value);
                    #ifdef INSERT_CLOCK
                        gettimeofday(&stop, NULL);
                        threshold_time += (double)(stop.tv_usec - start.tv_usec) +
                                    (double)(stop.tv_sec - start.tv_sec) * 1000000;
                    #endif
                    return ret;
                }
            }
        #endif
        uint8_t new_local_depth = cell[k].local_depth;
        // if local depth does not consist, reinsert.
        if(new_local_depth != local_depth)
            goto REINSERT;
        #ifdef INSERT_CLOCK
            gettimeofday(&stop, NULL);
            insert_time += (double)(stop.tv_usec - start.tv_usec) +
                                (double)(stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        bucket_.get_lock();
        #ifdef INSERT_CLOCK
            gettimeofday(&stop, NULL);
            collect_time += (double)(stop.tv_usec - start.tv_usec) +
                                (double)(stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        result = bucket_.insert_kv_to_bucket(key, value, k, init_guide);
        if(result == 1)
        {
            // if local depth does not consist, reinsert.
            uint8_t new_local_depth = cell[k].local_depth;
            if(new_local_depth != local_depth)
                goto REINSERT;
            insert_dire_cnt ++;
            #ifdef INSERT_CLOCK
                gettimeofday(&stop, NULL);
                insert_time += (double)(stop.tv_usec - start.tv_usec) +
                                (double)(stop.tv_sec - start.tv_sec) * 1000000;
                gettimeofday(&start, NULL);
            #endif
            bucket_.release_lock();
            #ifdef INSERT_CLOCK
                gettimeofday(&stop, NULL);
                collect_time += (double)(stop.tv_usec - start.tv_usec) +
                                    (double)(stop.tv_sec - start.tv_sec) * 1000000;
                gettimeofday(&start, NULL);
            #endif
            return true;
        }else if(result == 0)
        {
            bucket_.release_lock();
            return false; // redundant insertion
        }
        #ifdef INSERT_CLOCK
            gettimeofday(&start, NULL);
        #endif
        bucket_.release_lock();
        #ifdef INSERT_CLOCK
            gettimeofday(&stop, NULL);
            collect_time += (double)(stop.tv_usec - start.tv_usec) +
                                (double)(stop.tv_sec - start.tv_sec) * 1000000;
        #endif
        /*--------------------------------stage 2-------------------------------*/
        // directly insert failed, collect items and adjust cell state
        uint32_t max_pos = 0;
        uint32_t max_guide = 0;
        #ifdef INSERT_CLOCK
            gettimeofday(&start, NULL);
        #endif
        if(local_depth == global_depth)
        {
            if(is_movable(k, d, init_guide, local_depth, global_depth, 
                    init_cell_number, &max_pos, &max_guide))
            {
                // circular move
                #ifdef DIRECTORY
                    cell_->guide_offset[init_guide] += (max_guide > guide)? 
                        (max_guide - guide) : (cell_hash - (guide - max_guide));
                #else
                    lock_cells();
                    uint8_t new_offset = (cell[k].guide_offset[init_guide] + init_guide) % cell_hash;
                    uint8_t new_global_depth = layers->global_depth;
                    if(new_offset == max_guide || new_global_depth != global_depth)
                    {
                        unlock_cells();
                        goto REDO;
                    }
                    cell[k].guide_offset[init_guide] += (max_guide > guide)? 
                        (max_guide - guide) : (cell_hash - (guide - max_guide));
                    unlock_cells();
                #endif
                bucket_.get_lock();
                bucket[max_pos].get_lock();
                ret &= fine_grained_circular_move(k, local_depth, init_guide, init_cell_number,
                                                bucket + d, bucket + max_pos, key, value);
                bucket_.release_lock();
                bucket[max_pos].release_lock();
                move_cnt ++;
                #ifdef INSERT_CLOCK
                    gettimeofday(&stop, NULL);
                    move_time += (double)(stop.tv_usec - start.tv_usec) +
                                        (double)(stop.tv_sec - start.tv_sec) * 1000000;
                #endif
                return true;
            }else
            {
                // dynamic extension
                int r = before_dynamic_extend(bucket + d, global_depth, cell_number, level);
                if(r == 1)
                {
                    threshold_cnt ++;
                    result = bucket_.insert_kv_to_bucket(key, value, k, init_guide);
                    if(result <=0 )
                        return false;
                    else 
                        return true;
                }else if(r == 0)
                    goto REDO;
                
                lock_cells();
                uint8_t new_global_depth = layers->global_depth;
                if(new_global_depth != global_depth)
                {
                    unlock_cells();
                    goto REDO;
                }
                dynamic_extension(global_depth);
                #ifdef INSERT_CLOCK
                    gettimeofday(&stop, NULL);
                    extend_time += (double)(stop.tv_usec - start.tv_usec) +
                                        (double)(stop.tv_sec - start.tv_sec) * 1000000;
                #endif
                unlock_cells();
                global_depth += 1;
                cell_number *= 2;
                std::cout << "extension time: " <<  extend_time/CLOCKS_PER_SEC << std::endl;
                int diff_depth = global_depth - local_depth;
                int diff_cell = cell_number >> diff_depth;
                int first_offset = k % diff_cell;
                lock_cells();
                new_global_depth  = layers->global_depth;
                new_local_depth = cell[k].local_depth;
                if(new_global_depth != global_depth || new_local_depth != local_depth)
                {
                    unlock_cells();
                    goto REDO;
                }
                for(int i = 0; i < 1 << diff_depth; i ++)
                    cell[first_offset + i * diff_cell].local_depth ++;
                unlock_cells();
                ret &= rehash_a_cell(k, cell_bucket_, guide_offset_, 
                        local_depth, global_depth, cell_number, level);
                ret &= insert(key, value);// goto
                extend_cnt ++;
                return ret;
            }
        }else if(local_depth < global_depth)
        {
            int diff_depth = global_depth - local_depth;
            int diff_cell = cell_number >> diff_depth;
            int first_offset = k % diff_cell;
            lock_cells();
            uint8_t new_global_depth = layers->global_depth;
            new_local_depth = cell[k].local_depth;
            if(new_global_depth != global_depth || new_local_depth != local_depth)
            {
                unlock_cells();
                goto REDO;
            }
            for(int i = 0; i < 1 << diff_depth; i ++)
                cell[first_offset + i * diff_cell].local_depth ++;
            unlock_cells();
            ret &= rehash_a_cell(k, cell_bucket_, guide_offset_, 
                    local_depth, global_depth, cell_number, level);
            ret &= insert(key, value);// goto
            rehash_cnt ++;
            #ifdef INSERT_CLOCK
                gettimeofday(&stop, NULL);
                threshold_time += (double)(stop.tv_usec - start.tv_usec) +
                                    (double)(stop.tv_sec - start.tv_sec) * 1000000;
            #endif
            return ret;
            
        }
        return false;
    }

    bool query(T key, Value_t value = NULL){
        #ifdef SEARCH_CLOCK
            flag = 1;
            query_cnt ++;
        #endif
        int k, init_guide, guide, d;
        #ifdef SEARCH_CLOCK
            timeval start, stop;
            gettimeofday(&start, NULL);
        #endif
        uint8_t global_depth = layers->global_depth;
        uint32_t cell_number = init_cell_number << global_depth;
        // hash twice
        k = calculate_cell(&key) % cell_number;
        init_guide = calculate_guide(&key);
        #ifdef SEARCH_CLOCK
            gettimeofday(&stop, NULL);
            cell_time += (double)(stop.tv_usec - start.tv_usec) +
                        (double) (stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        // read NVM/DRAM once
        #ifdef DIRECTORY
            Cell *cell_ = dir[k].cell_ptr;
            uint8_t local_depth = dir[k].local_depth;
            guide = (init_guide + cell_->guide_offset[init_guide]) % cell_hash;
            d = cell_->cell_bucket[guide];

        #else
            Cell cell_ = cell[k];
            uint8_t local_depth = cell_.local_depth;
            guide = (init_guide + cell_.guide_offset[init_guide]) % cell_hash;
            d = cell_.cell_bucket[guide];
        #endif
        #ifdef SEARCH_CLOCK
            gettimeofday(&stop, NULL);
            init_guide_time += (double)(stop.tv_usec - start.tv_usec) +
                        (double) (stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        // hash once and read NVM/DRAM once
        int level;
        calculate_level(k, &level, local_depth, global_depth);
        #ifdef SEARCH_CLOCK
            gettimeofday(&stop, NULL);
            cell_time += (double)(stop.tv_usec - start.tv_usec) +
                        (double) (stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        auto bp = layers->bp;
        auto cur_bp = bp[level];
        Bucket<T> *bucket = reinterpret_cast<Bucket<T> *>(reinterpret_cast<uint64_t>(cur_bp));

        #ifdef SEARCH_CLOCK
            gettimeofday(&stop, NULL);
            hash_time += (double)(stop.tv_usec - start.tv_usec) + 
                        (double)(stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        // Bucket<T> *bucket = layers->bp[level];
QUERY:
        uint32_t  pre_version, cur_version;
        bool ret;
        Bucket bucket_ = bucket[d];
        // Bucket<T> *bucket_ = bucket + d;
        #ifdef SEARCH_CLOCK
            gettimeofday(&stop, NULL);
            bucket_time += (double)(stop.tv_usec - start.tv_usec) +
                        (double) (stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        pre_version = __atomic_load_n(&bucket_.version_lock, __ATOMIC_ACQUIRE);
        #ifdef SEARCH_CLOCK
            gettimeofday(&stop, NULL);
            level_time += (double)(stop.tv_usec - start.tv_usec) +
                        (double) (stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        // read NVM/DRAM once
        ret = bucket_.query_key_in_bucket(key, value);
        #ifdef SEARCH_CLOCK
            gettimeofday(&stop, NULL);
            search_time += (double)(stop.tv_usec - start.tv_usec) +
                        (double) (stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        cur_version = __atomic_load_n(&bucket_.version_lock, __ATOMIC_ACQUIRE);
        if(pre_version != cur_version)
        {
            printf("reach here.\n");
            goto QUERY;
        }
        #ifdef SEARCH_CLOCK
            gettimeofday(&stop, NULL);
            level_time += (double)(stop.tv_usec - start.tv_usec) +
                        (double) (stop.tv_sec - start.tv_sec) * 1000000;
        #endif
        // if(ret == false)
        //     printf("herrrrr\n");
        return ret;
FINAL:
        return NONE;
    }

    bool remove(T key)
    {
        int k, init_guide, guide, d;
        bool ret = true;
        #ifdef DELETE_CLOCK
            timeval start, stop;
            gettimeofday(&start, NULL);
        #endif
        uint8_t global_depth = layers->global_depth;
        uint32_t cell_number = init_cell_number << global_depth;
        k = calculate_cell(&key) % cell_number;
        init_guide = calculate_guide(&key);
        #ifdef DELETE_CLOCK
            gettimeofday(&stop, NULL);
            dcell_time += (double)(stop.tv_usec - start.tv_usec)
                        + (double)(stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        #ifdef DIRECTORY
            Cell *cell_ = dir[k].cell_ptr;
            uint8_t local_depth = dir[k].local_depth;
        #else
            Cell cell_ = cell[k];
            uint8_t local_depth = cell_.local_depth;
            guide = (init_guide + cell_.guide_offset[init_guide]) % cell_hash;
            d = cell_.cell_bucket[guide];
            // d = calculate_bucket(guide, k, cell_.local_depth);
        #endif
        #ifdef DELETE_CLOCK
            gettimeofday(&stop, NULL);
            dinit_guide_time += (double)(stop.tv_usec - start.tv_usec)
                            +   (double)(stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        int level;
        calculate_level(k, &level, local_depth, global_depth);
        #ifdef DELETE_CLOCK
            gettimeofday(&stop, NULL);
            dcell_time += (double)(stop.tv_usec - start.tv_usec)
                            +   (double)(stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        Bucket<T> *bucket = layers->bp[level];
        #ifdef DELETE_CLOCK
            gettimeofday(&stop, NULL);
            dohter_time += (double)(stop.tv_usec - start.tv_usec)
                            +(double)(stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        // Bucket<T> *bucket_ = bucket + d;
        Bucket<T> bucket_ = bucket[d];
        bucket_.get_lock();
        #ifdef DELETE_CLOCK
            gettimeofday(&stop, NULL);
            dbucket_time += (double)(stop.tv_usec - start.tv_usec)
                            +(double)(stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        ret &= bucket_.remove_record_from_bucket(key);
        #ifdef DELETE_CLOCK
            gettimeofday(&stop, NULL);
            dsearch_time += (double)(stop.tv_usec - start.tv_usec)
                            +(double)(stop.tv_sec - start.tv_sec) * 1000000;
            gettimeofday(&start, NULL);
        #endif
        bucket_.release_lock();
        #ifdef DELETE_CLOCK
            gettimeofday(&stop, NULL);
            dbucket_time += (double)(stop.tv_usec - start.tv_usec)
                            +(double)(stop.tv_sec - start.tv_sec) * 1000000;
        #endif
        return ret;
    }
    
    int insert_seq(T key, Value_t value, uint64_t cur_cnt)
    {
        int d = cur_cnt / BUCKET_CAPACITY;
        Bucket<T> *bucket = layers[0];
        int ret = bucket[d].insert_kv_to_bucket(key, value, 0, 0);
        return ret;
    }

    int insert_rnd(T key, Value_t value)
    {
        int d = rand() % 500000;
        Bucket<T> *bucket = layers[0];
        int ret = bucket[d].insert_kv_to_bucket(key, value, 0, 0);
        return ret;
    }

    Value_t get_seq(T key, uint64_t cur_cnt)
    {
        int d = cur_cnt / BUCKET_CAPACITY;
        Bucket<T> *bucket = layers[0];
        bool ret = bucket[d].query_key_in_bucket(key, NULL);
        if(ret)
            return (Value_t) 1;
        else
            return (Value_t) 0;
    }

    Value_t get_rnd(T key)
    {
        int d = rand() % 500000;
        Bucket<T> *bucket = layers[0];
        if(bucket[d].query_key_in_bucket(key, NULL))
            return (Value_t) 1;
        else
            return (Value_t) 0;
    }
    
// below are statistics function
public:
    int statistics_bucket_items(){
        int i, j;
        uint8_t global_depth = layers->global_depth;
        uint32_t items_num = 0;
        for(i = 0; i < 1 << global_depth; i ++)
        {
            Bucket<T> *bucket = layers->bp[i];
            for(j = 0; j < init_bucket_number; j ++){
                bucket_items_per_level[i] += bucket[j].used;
                items_num += bucket[j].used;
            }
        }
        return items_num;
    }

    double statistics_bit_per_item(uint32_t cell_number){
        int items_num = statistics_bucket_items();
        int cell_bit_sum = 0;
        for(int i = 0; i < 1; ++i)
            cell_bit_sum += cell_number * log(cell_hash);
        return (double)cell_bit_sum / items_num;
    }

    double statistics_load_factor(uint32_t bucket_number){
        int items_num = statistics_bucket_items();
        int bucket_slots = bucket_number * BUCKET_CAPACITY;
        return (double)items_num / bucket_slots;
    }
};

template<class T>
bool EEPH<T>::find_kv(T key)
{
    uint8_t global_depth = layers->global_depth;
    uint32_t bucket_number = init_bucket_number << global_depth;
    uint8_t layer_number = 1 << global_depth;
    for(int i = 0; i < layer_number; i++)
    {
        Bucket<T> *bp = layers->bp[i];
        for(int j = 0; j < init_bucket_number; j++)
        {
            Bucket<T> bucket = bp[j];
            for(int k = 0; k < BUCKET_CAPACITY; k++)
            {
                int valid = bucket.bitmap & (1 << k);
                if(valid && bucket.kv[k].key == key)
                    return true;
            }
        }
    }
    return false;
}

template<class T>
void EEPH<T>::recovery()
{
    #ifdef RECOVERY
    timeval start, stop;
    gettimeofday(&start, NULL);
    #endif
    initialize_hash_functions();
    // if(find_kv((T)15336469859637867841))
    //     printf("herrrrrrr\n");
    uint8_t global_depth = layers->global_depth;
    uint32_t cell_number_ = init_cell_number << global_depth;
    cell = (Cell *)malloc(sizeof(Cell) * cell_number_);
    memset(cell, 0, sizeof(Cell) * cell_number_);
    for(int i = 0; i < cell_number_; i++)
        cell[i].local_depth = 255;
    for(int i = 0; i < cell_number_; i++)
    {
        // if(i == 544)
        //     printf("herrrr\n");
        Cell *cell_ = cell + i;
        uint8_t layer_num = 1 << global_depth;
        uint8_t ld[16] = {0};
        // uint8_t *ld = (uint8_t *)malloc(sizeof(uint8_t) * layer_num);
        // memset(ld, 0, sizeof(uint8_t) * layer_num);
        if(cell_->local_depth != 255)
            continue;
        for(int j = 0; j < cell_hash; j++)
        {
            uint8_t c = 0;
            uint32_t d = calculate_bucket(j, i, 0);
            cell_->cell_bucket[j] = d;
            uint8_t layer = i / init_cell_number;
            Bucket<T> *bucket  = layers->bp[layer];
            Bucket<T>  bucket_ = bucket[d]; 

            for(int k = 0; k < BUCKET_CAPACITY; k ++)
            {
                // if(bucket_.kv[k].key == 15336469859637867841 && i == 1006)
                //     printf("herrrrr\n");
                int valid = bucket_.bitmap & (1 << k);
                if(valid && bucket_.k[k] % init_cell_number == i % init_cell_number)
                {
                    int cell_index = calculate_cell(&bucket_.kv[k].key) % cell_number_;
                    layer = cell_index / init_cell_number;
                    ld[layer] = 1;

                    // recovery offset array
                    uint8_t virtual_bucket_index = bucket_.init_guide[k];
                    if(virtual_bucket_index == j)
                        cell_->guide_offset[j] = 0;
                    else
                        cell_->guide_offset[virtual_bucket_index] = j > virtual_bucket_index? 
                            j - virtual_bucket_index : cell_hash + j - virtual_bucket_index;
                }
            }
            
        }
        // recovery local depth
        uint32_t x = 0;
        for(int j = 0; j < layer_num; j++)
            x += ld[j];
        uint8_t local_depth;
        if(x == 0)
            local_depth = 0;
        else if(x == 1)
            local_depth = global_depth;
        else
            local_depth = global_depth - ceil(log2(x));
        // local_depth = log2(x);
        assert(local_depth <= global_depth);
        uint8_t diff_depth  = global_depth - local_depth;
        uint32_t diff_cell  = cell_number_ >> diff_depth;
        cell_->local_depth = local_depth;
        if(diff_depth > 0 &&  i + 1 * diff_cell >= cell_number_)
            continue;
        for(int j = 1; j < pow(2, diff_depth); j++)
        {
            assert(i + j * diff_cell < cell_number_);
            // memcpy(cell+ i + k * diff_cell, cell_, sizeof(Cell));
            // memcpy(cell[i + k * diff_cell].cell_bucket, cell_->cell_bucket, sizeof(cell_->cell_bucket));
            memcpy(cell[i + j * diff_cell].cell_bucket, cell_->cell_bucket, sizeof(uint32_t) * cell_hash);
            memcpy(cell[i + j * diff_cell].guide_offset, cell_->guide_offset, sizeof(uint8_t) * cell_hash);
            cell[i + j * diff_cell].local_depth = local_depth;
        }
    }

    #ifdef RECOVERY
    gettimeofday(&stop, NULL);
    double recovery_time = (double)(stop.tv_usec - start.tv_usec)
                            +   (double)(stop.tv_sec - start.tv_sec) * 1000000;
    std::cout << "recovery time(s): " << (double) recovery_time / CLOCKS_PER_SEC << "s" << std::endl;
    #endif
}

template<class T>
EEPH<T>::EEPH(void){
    std::cout << "Rebuild Cells in DRAM..."<< std::endl; 
    recovery();
};

template<class T>
EEPH<T>::EEPH(uint32_t bucket_number_, uint32_t cell_number_, uint8_t cell_hash_){
    std::cout << "size of bucket " << sizeof(Bucket<T>) << std::endl;
    init_bucket_number = bucket_number_;
    init_cell_number   = cell_number_;
    cell_hash          = cell_hash_;

    PMEMoid oid;
    Layer<T>::New(&layers_, 0, 0);
    layers = (Layer<T> *)pmemobj_direct(layers_);
    
    Bucket<T>::New(&oid, bucket_number_, 0);
    layers->bp[0] = (Bucket<T> *)pmemobj_direct(oid);
    initialize_hash_functions();
    #ifdef DIRECTORY
        dir = (Directory *)malloc(sizeof(Directory) * cell_number);
        for(int i = 0; i < cell_number; i++)
        {
            dir[i].cell_ptr = (Cell *)malloc(sizeof(Cell));
            dir[i].local_depth = 0;
            Cell *cell = dir[i].cell_ptr;
            memset(cell->guide_offset, 0, sizeof(uint8_t) * 16);
            for(int j = 0; j < cell_hash; j++)
            {
                cell->cell_bucket[j] = calculate_bucket(j, i, 0);
            }
        }
    #else
        cell = (Cell *)malloc(sizeof(Cell) * cell_number_);
        for(int i = 0; i < cell_number_; i++)
        {
            memset(cell[i].guide_offset, 0, sizeof(uint8_t) * 16);
            for(int j = 0; j < cell_hash; j++)
            {
                cell[i].cell_bucket[j] = calculate_bucket(j, i, 0);
            }
            cell[i].local_depth = 0;
        }
    #endif
    // #endif

    capacity = bucket_number_ * BUCKET_CAPACITY;
    bucket_items = 0;
    memset(bucket_items_per_level, 0, sizeof(bucket_items_per_level));
}

template <class T>
EEPH<T>::~EEPH(void){
    PMAllocator::Close_pool();
    free(cell);
}

template <class T>
inline int EEPH<T>::Insert(T key, Value_t value)
{
    // struct KV_entry *kv_entry = generate_kv_pair(key, value);
    int ret = insert(key, value);
    return ret;
}

template <class T>
int EEPH<T>::Insert(T key, Value_t value, bool is_in_epoch)
{
    return Insert(key, value);
}

template <class T>
inline bool EEPH<T>::Delete(T key)
{
    return remove(key);
}

template <class T>
bool EEPH<T>::Delete(T key, bool is_in_epoch)
{
    return remove(key);
}

template <class T>
inline Value_t EEPH<T>::Get(T key)
{
    if(query(key, nullptr))
        return (Value_t) 1;
    else    
        return (Value_t) 0;
}

template <class T>
Value_t EEPH<T>::Get(T key, bool is_in_epoch)
{
    return Get(key);
}

template <class T>
int EEPH<T>::Insert_seq(T key, Value_t value, bool is_in_epoch, uint64_t cur_cnt)
{
    return insert_seq(key, value, cur_cnt);
}
template <class T>
int EEPH<T>::Insert_rnd(T key, Value_t value, bool is_in_epoch)
{
    return insert_rnd(key, value);
}

template <class T>
Value_t EEPH<T>::Get_seq(T key, bool is_in_epoch, uint64_t cur_cnt)
{
    return get_seq(key, cur_cnt);
}
template <class T>
Value_t EEPH<T>::Get_rnd(T key, bool is_in_epoch)
{
    return get_rnd(key);
}


}// namespace eeph

#endif
