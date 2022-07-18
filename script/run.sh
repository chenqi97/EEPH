#!/bin/bash

# number of threads to run
thread_num=(0 1 4 8 16 24 32)
# benckmark workload, number of opeartions to run
workload=(0 19000000 2000000 40000000 60000000 80000000)
# warm-up workload, number of key-value to insert for warm-up
base=(0 1000000)
# type of keys
key_type=(0 fixed variab1e)
# which index to evaluate
index_type=(0 eeph dash cceh level)
# whether to use to epcoh manager
epoch=(0 1 1 1 1)
# operations
ops=(full insert pos delete neg mixed)

for k in 1  #index
do
	for i in 1  #only fixed kv for eeph
	do
		for j in 1  #thread number
    do
      for p in 1 #operation number
      do
        echo "Begin: ${base[1]} ${workload[p]} ${thread_num[${j}]}"
        rm -f /mnt/pmem0/cq-hash/*
        LD_PRELOAD="../pmdk/src/PMDK/src/nondebug/libpmemobj.so.1 \
        ../pmdk/src/PMDK/src/nondebug/libpmem.so.1" \
        ../build/EEPH \
        -n ${base[1]} \
        -loadType 0 \
        -p ${workload[p]} \
        -t ${thread_num[$j]} \
        -k ${key_type[$i]} \
        -distribution "uniform" \
        -index ${index_type[$k]} \
        -e ${epoch[$k]} \
        -ed 1000 \
        -op "full" \
        -r 0.8 \
        -s 0.2 \
        -ms 100 \
        -ps 60
      done
    done
	done
done

