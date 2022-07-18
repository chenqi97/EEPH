#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <sys/time.h>
#include <fstream>
#include <cstring>
#include <gflags/gflags.h>

#include "Hash.h"
#include "EEPH.h"
#include "Dash/ex_finger.h"
#include "CCEH/CCEH_baseline.h"
#include "Level/level_baseline.h"
#include "libpmemobj.h"
#include "PMAllocator.h"
typedef struct Record{
    enum Op : uint32_t {INSERT = 0, GET = 1, UPDATA = 2};
    Op op = INSERT;
    uint64_t key;
    uint64_t value;
} Record;

struct range {
  int index;
  uint64_t begin;
  uint64_t end;
  int length;
//   void *workload;
  std::vector<Record> *workload;
  uint64_t random_num;
  struct timeval tv;
};
 
static const char pool_name_eeph[] = "/mnt/pmem0/cq-hash/pmem_eeph.data";
static const char pool_name_dash[] = "/mnt/pmem0/cq-hash/pmem_dash.data";
static const char pool_name_cceh[] = "/mnt/pmem0/cq-hash/pmem_cceh.data";
static const char pool_name_level[] = "/mnt/pmem0/cq-hash/pmem_level.data";
static const char data_base_path[] = "/home/cq/datasets/ycsb/";
// std::filesystem::path dataset_path = "/home/cq/datasets/ycsb/raw_ycsb_wl_1090_uniform.dat";
static const size_t pool_size = 1024ul * 1024ul * 1024ul * 30ul;

DEFINE_string(index, "dash", "");
DEFINE_string(dataset, "ycsb_wl_1090_uniform.dat", "dataset name");
DEFINE_string(fv, "fixed", "fixed kv or variable kv");
DEFINE_uint32(thread, 1, "number of thread");

namespace param
{
    std::string index_name, dataset_name, fv;
    uint32_t thread_num;
}

void set_affinity(uint32_t idx) {
  cpu_set_t my_set;
  CPU_ZERO(&my_set);
  if (idx < 24) {
    CPU_SET(idx, &my_set);
  } else {
    CPU_SET(idx + 24, &my_set);
  }
  sched_setaffinity(0, sizeof(cpu_set_t), &my_set);
}

static bool file_exists(const char *pool_path)
{
    struct stat buffer;
    return (stat(pool_path, &buffer) == 0);
}

void load_dataset(std::string wl_path, std::vector<Record> *data)
{
    std::cout << "Loading Dataset ..." << std::endl;
    std::cout << "path: " << wl_path << std::endl;
    std::ifstream ifs{wl_path, std::ios::binary | std::ios::ate};
    if(!ifs)
        std::cout << "Cannot open file: " << wl_path << std::endl;
    auto end = ifs.tellg();
    ifs.seekg(0, std::ios::beg);
    auto size = std::size_t(end - ifs.tellg());

    if(size == 0)
        std::cout << "empty file" << std::endl;
    const uint64_t number_records = size / sizeof(Record);
    data->resize(number_records);

    if(!ifs.read((char *)data->data(), size))
        std::cout << "Error reading from " << wl_path << std::endl;
    std::cout << "Loading Done. " << data->size() << " kvs are loaded." << std::endl;
    // std::cout << "records: " << cnt << std::endl;
}

template <class T>
void concurr_operation(struct range *range, Hash<T> *index)
{
    set_affinity(range->index);
    uint64_t begin = range->begin;
    uint64_t end = range->end;
    std::vector<Record> dataset = *range->workload;
    int insert_failed = 0, insert_success = 0, not_found = 0;
    int insert_cnt = 0, get_cnt = 0, update_cnt = 0, default_cnt = 0;
    int ret;
    struct timespec start = {0, 0};
    struct timespec stop  = {0, 0};
    for(uint64_t i = begin; i < end; i++)
    {
        const Record& record = dataset[i];
        switch(record.op)
        {
            case Record::INSERT:
                ret = index->Insert(record.key, (Value_t)1);
                if(ret == 0)
                    insert_success ++;
                else
                    insert_failed ++;
                insert_cnt ++;
                break;
            case Record::GET:
                clock_gettime(CLOCK_REALTIME, &start);
                if(index->Get(record.key, true) == NONE)
                    not_found ++;
                clock_gettime(CLOCK_REALTIME, &stop);
                if(i % 150000 == 0)
                {
                    double latency = (double)(stop.tv_nsec - start.tv_nsec) +
                            (double)(stop.tv_sec - start.tv_sec) * 1000000000;
                    // std::cout << "Latency: " << (double)latency << "ns" << std::endl;
                }
                get_cnt ++;
                break;
            case Record::UPDATA:
                update_cnt ++;
                break;
            default:
                default_cnt ++;
        }
    }
    gettimeofday(&range->tv, NULL);
    std::cout << "insert failed: " << insert_failed << " insert success: " << insert_success 
        << std::endl;
    std::cout << "search not found: " << not_found << " search cnt: " << get_cnt << std::endl; 
    std::cout << "update cnt: " << update_cnt << " default cnt: " << default_cnt << std::endl;
}

template <class T>
void run_benchmark(std::vector<Record> *data, Hash<T> *index, uint64_t dataset_size,
        void (*test_func)(struct range *, Hash<T> *))
{
    uint32_t thread_num = param::thread_num;
    uint32_t chunk_size = dataset_size / thread_num;
    std::thread *thread_array[128];
    struct range *rarray = reinterpret_cast<range *>(malloc(thread_num * sizeof(struct range)));
    timeval tv1;
    double duration;
    for(uint32_t i = 0; i < thread_num; i++)
    {
        rarray[i].index = i;
        rarray[i].random_num = rand();
        rarray[i].begin = i * chunk_size;
        rarray[i].end = (i + 1) * chunk_size;
        rarray[i].length = 8;
        rarray[i].workload = data;
    }
    rarray[thread_num - 1].end = dataset_size;
    for(uint32_t i = 0; i < thread_num; i++)
    {
        thread_array[i] = new std::thread(*test_func, rarray + i, index);
    }
    gettimeofday(&tv1, NULL);
    for(uint32_t i = 0; i < thread_num; i ++)
    {
        thread_array[i]->join();
        delete thread_array[i];
    }
    double longest = (double)(rarray[0].tv.tv_usec - tv1.tv_usec) / 1000000 +
                   (double)(rarray[0].tv.tv_sec - tv1.tv_sec);
    double shortest = longest;
    duration = longest;

    for (int i = 1; i < thread_num; ++i) {
        double interval = (double)(rarray[i].tv.tv_usec - tv1.tv_usec) / 1000000 +
                        (double)(rarray[i].tv.tv_sec - tv1.tv_sec);
        duration += interval;
        if (shortest > interval) shortest = interval;
        if (longest < interval) longest = interval;
    }
    duration = duration / thread_num;
    std::cout << thread_num << " threads, Time = " << duration << " s, throughput = " << dataset_size / duration / 1000000
        << " Mops/s, fastest = " << dataset_size / shortest / 1000000 << ", slowest = " << dataset_size / longest / 1000000 
         << " Mops/s" << std::endl;
    std::cout << "==========================================================" << std::endl;
}

template <class T>
void init_index(std::string index_name, std::string dataset_name, std::string fv, uint64_t thread_num){
    std::vector<Record> dataset;
    Hash<uint64_t> *eh;
    bool pool_exist = false;
    // int insert_failed = 0, insert_success = 0, not_found = 0;
    // int insert_cnt = 0, get_cnt = 0, update_cnt = 0, default_cnt = 0;

    // load dataset
    const std::string tmp_path = data_base_path + std::string{dataset_name};
    load_dataset(tmp_path, &dataset);
    int dataset_size = dataset.size();
    std::cout << "============================= "<< index_name << " =============================" << std::endl;
    // load index
    if(index_name == "eeph")
    {
        if(fv == "fixed")
        {
            if(file_exists(pool_name_eeph)) 
                pool_exist = true;
            PMAllocator::Initialize(pool_name_eeph, pool_size);
            eh = reinterpret_cast<Hash<T> *>(PMAllocator::GetRoot(sizeof(eeph::EEPH<T>)));
            if(!pool_exist)
            {
                int bucket_number = 5000000;
                int cell_number = 10000000;
                int cell_hash = 16;
                new (eh) eeph::EEPH<T>(bucket_number, cell_number, cell_hash);  
            }else
            {
                std::cout << "get existed eeph object" << std::endl;
                new (eh) eeph::EEPH<T>();
            }
        }else
        {
            /* variable kv */
        }
        
    }else if(index_name == "dash")
    {
        if(fv == "fixed")
        {
            pool_exist = false;
            if(file_exists(pool_name_dash)) pool_exist = true;
            PMAllocator::Initialize(pool_name_dash, pool_size);
            eh = reinterpret_cast<Hash<T> *>(PMAllocator::GetRoot(sizeof(extendible::Finger_EH<T>)));
            if(!pool_exist)
                new (eh) extendible::Finger_EH<T>(65536, PMAllocator::Get()->pm_pool_);
            else
                new (eh) extendible::Finger_EH<T>();
        }else
        {
            /* variable kv */
        }
    }else if(index_name == "cceh")
    {
        if(fv == "fixed")
        {
            pool_exist = false;
            if (FileExists(pool_name_cceh)) pool_exist = true;
            PMAllocator::Initialize(pool_name_cceh, pool_size);
            eh = reinterpret_cast<Hash<T> *>(PMAllocator::GetRoot(sizeof(cceh::CCEH<T>)));
            if (!pool_exist) {
                new (eh) cceh::CCEH<T>(262144, PMAllocator::Get()->pm_pool_);
            } else {
                new (eh) cceh::CCEH<T>();
            }
        }else
        {
            /* variable kv */
        }
    }else if(index_name == "level")
    {
        if(fv == "fixed")
        {
            pool_exist = false;
            if (FileExists(pool_name_level)) pool_exist = true;
            PMAllocator::Initialize(pool_name_level, pool_size);
            eh = reinterpret_cast<Hash<T> *>(
                PMAllocator::GetRoot(sizeof(level::LevelHashing<T>)));
            if (!pool_exist) {
                new (eh) level::LevelHashing<T>();
                int level_size = 22;
                level::initialize_level(PMAllocator::Get()->pm_pool_,
                                    reinterpret_cast<level::LevelHashing<T> *>(eh),
                                    &level_size);
            } else {
                new (eh) level::LevelHashing<T>();
            }
        }else
        {
            /* variable kv */
        }
    }
    run_benchmark(&dataset, eh, dataset_size, &concurr_operation);
    // insert operations
    // gettimeofday(&start, NULL);
    // for(int op_num = 0; op_num < dataset_size; op_num ++)
    // {
    //     const Record& record = dataset[op_num];
    //     switch(record.op)
    //     {
    //         case Record::INSERT:
    //             ret = eh->Insert(record.key, (Value_t)1);
    //             if(ret == 0)
    //                 insert_success ++;
    //             else
    //                 insert_failed ++;
    //             insert_cnt ++;
    //             break;
    //         case Record::GET:
    //             if(eh->Get(record.key, true) == NONE)
    //                 not_found ++;
    //             get_cnt ++;
    //             break;
    //         case Record::UPDATA:
    //             update_cnt ++;
    //             break;
    //         default:
    //             default_cnt ++;
    //     }
    // }
    // gettimeofday(&stop, NULL);
    // double insert_interval = (double)(stop.tv_usec - start.tv_usec) + 
    //                         (double)(stop.tv_sec - start.tv_sec) * 1000000;
    // std::cout << "throughput(Mops): " << (dataset_size - default_cnt) / insert_interval 
    // << " toal time: " << insert_interval / CLOCKS_PER_SEC << std::endl;
    // std::cout << "insert failed: " << insert_failed << " insert success: " << insert_success 
    //     << " total records: " << dataset.size() << std::endl;
    // std::cout << "search not found: " << not_found << " search cnt: " << get_cnt << std::endl; 
    // std::cout << "update cnt: " << update_cnt << " default cnt: " << default_cnt << std::endl;
    // std::cout << "============================= "<< index_name << " =============================" << std::endl;
}

int main(int argc, char* argv[]){
    google::ParseCommandLineFlags(&argc, &argv, true);
    param::index_name = FLAGS_index;
    param::dataset_name = FLAGS_dataset;
    param::fv = FLAGS_fv;
    param::thread_num = FLAGS_thread;
    init_index<uint64_t>(param::index_name, param::dataset_name, param::fv, param::thread_num);
    return 0;
}