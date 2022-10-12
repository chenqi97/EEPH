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
3.  cmake -DCMAKE_BUILD_TYPE=Release -DUSE_PMEM=ON .. or cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_PMEM=ON ..
4.  make

### Excute

1.  If successfully, you can get two executable files named 'EEPH' and 'ycsb_benchmark'. 'EEPH'  runs in Pibench, 'ycsb_benchmark' runs in YCSB datasets.
2.  Except the executable files, '../script/run.sh' is a simple script for testing.
3.  xxxx

### Structure
We have included three comparison hash indexes(Dash, CCEH, and Level) in this project, whose header file is in the folder with the same name. For Dash-hybrid and CCEH-hybrid version, you can set "DASH_HYBRID" (lines 87 in ex_finger.h) or "CCEH_HYBRID" (lines 41 in CCEH.h) in the source code. EEPH.h is the header file of our proposal. You can download the header files in your codes for testing. Also, you can check the implements of algorithm mentioned in paper (Section IV): complement move (impemented in function "is_movable"), local rehashing (implemented in function "rehash_a_cell"), and extension (implemented in function "dynamic_extension").
