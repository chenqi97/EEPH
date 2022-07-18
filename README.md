# PH_on_PM

### Introduction
Perfect Hashing on Persist Memory

### Hardware Environment
1. create region  
    ipmctl create -goal PersistentMemoryType=AppDirect(create region in interleave mode)
2. create namespace  
    ndctl create-namespace --mode=fsdax(devdax) --size=768G --region=region0 --force
3. create file system  
    mkfs.ext4 /dev/pmem0
4. mount  
    mount -o dax /dev/pmem0 /mnt/pmem0
5. verify  
    mount -v | grep pmem  
    /dev/pmem0 on /mnt/pmem0 type ext4 (rw,relatime,seclabel,dax,data=ordered)

### Software Environment
1. c++
2. pmdk
3. other dependency library

### Install

1.  mkdir build
2.  cd build
3.  cmake -DCMAKE_BUILD_TYPE=Release -DUSE_PMEM=ON .. or cmake -DCMAKE_BUILD_TYPE=Debug-DUSE_PMEM=ON ..
4.  make

### Excute

1.  If successfully, you can get two executable files named 'EEPH' and 'ycsb_benchmark'. 'EEPH'  runs in Pibench, 'ycsb_benchmark' runs in YCSB datasets.
2.  Except the executable files, '../build/run.sh' is a simple script for testing.
3.  xxxx
