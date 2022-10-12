# PH_on_PM

### Introduction

Perfect Hashing on Persistent Memory

### Hardware Environment

Before running the project, you should create a filesystem using the /dev/pmem device and mount it by the following four steps. Finally, you can verify the device is mounted successfully.

```
create a region
$ ipmctl create -goal PersistentMemoryType=AppDirect(create region in interleave mode)
```

```
create a namespace
$ ndctl create-namespace --mode=fsdax(devdax) --size=768G --region=region0 --force
```

```
create a filesystem for the device
$ mkfs.ext4 /dev/pmem0
```

```
mount it to a mount point
$ mount -o dax /dev/pmem0 /mnt/pmem0
```

```
verify the result
$ mount -v | grep pmem  
/dev/pmem0 on /mnt/pmem0 type ext4 (rw,relatime,seclabel,dax,data=ordered)
```

### Software Environment

1. GCC
2. pmdk
3. other dependency library

### Install

```
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release -DUSE_PMEM=ON .. or cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_PMEM=ON ..
$ make
```

### Excute

If the previous steps are executed successfully, you can get two executable files named `EEPH` and `ycsb_benchmark`. `EEPH` runs in Pibench, and `ycsb_benchmark` runs in YCSB datasets. Except for the executable files, `../script/run.sh` is a simple script to run in Pibench, and `../script/run_ycsb.sh` is running in YCSB. Note that you can download CAIDA datasets from https://www.caida.org/catalog/datasets/ipv6_dnsnames_dataset/ for further testing.

### Structure

**Header files:** In this project, we have included three comparison hash indexes(Dash, CCEH, and Level) whose header file is in the folder with the same name. For Dash-hybrid and CCEH-hybrid versions, you can set `DASH_HYBRID` (lines 87 in `ex_finger.h`) or `CCEH_HYBRID` (lines 41 in `CCEH.h`) in the source code. `EEPH.h` is the header file of our proposal. You can download the header file in your own codes for testing. Also, you can check the implements of the algorithm mentioned in the paper (Section IV): complement move (implemented in function `is_movable`), local rehashing (implemented in function `rehash_a_cell`), and extension (implemented in function `dynamic_extension`). `xxhash.h` and `murmur3.h` are two hash function headers we use in EEPH (`murmur3` is used default). You can switch to xxhash by turning on the `XXHASH` switch (lines 51 in `EEPH.h`).

**CPP files:** `main.cpp` and  `test_pmem.cpp` are the two main files to run in YCSB and Pibench, respectively.
