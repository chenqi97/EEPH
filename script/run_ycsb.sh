thread_num=(0 1 4 8 16 24 32)

workload=(0 dns_out ycsb_wl_1090_uniform.dat ycsb_wl_100_uniform.dat ycsb_wl_5050_zipf.dat ycsb_wl_100_zipf.dat)

index_type=(0 eeph dash cceh level)

for i in {1..4}  # index type
do
    for j in 3 # workload
    do
        for k in 4  # thread_num
        do
            echo "TEST STARTING: ${index_type[i]} ${workload[j]} ${thread[k]}"
            rm -f /mnt/pmem0/cq-hash/*
            ../build/ycsb_benchmark \
            -index ${index_type[i]} \
            -dataset ${workload[j]} \
            -thread ${thread_num[k]}
        done
    done
done
