## Initial development log

1. [Comparison of speeds [2020-06-17]](#log_swe9n_v0002_1)
1. [OpenACC for Time-Loop Attempt 1 [2020-06-17][2020-06-18]](#log_swe9n_v0002_2)
1. [OpenACC for GWCErh2() [2020-06-18]](#log_swe9n_v0002_3)
1. [Tuning num_workers and vector_length [2020-06-17]](#log_swe9n_v0002_4)

### Attempting
- OpenACC directives based GPU parallelisation


### List of Work
- [ ] To be added

-----------------------------------------------

<a name = 'log_swe9n_v0002_4' ></a>

### Tuning num_workers and vector_length [2020-06-17]
- A lot of gains can be made by playing around with gangs, workers and vectors.
- This subroutine gives a good test case.
- The basic thing i tried was which takes 4.56s for TestCase1.

```
Code:
  !$acc parallel loop default(present) gang vector &
  !$acc   private(k, j, neid, etaDx, etaDy, lDep) 
  do k = 1, npt    

    !jxTil(k) = 0d0
    !jyTil(k) = 0d0
    lDep = dep(k)
    etaDx = 0d0    
    etaDy = 0d0    
    do j = 1, pObj(k)%nn
      neid = pObj(k)%neid(j)
      etaDx = etaDx + pObj(k)%phiDx(j)*eta(neid)
      etaDy = etaDy + pObj(k)%phiDy(j)*eta(neid)      
    enddo
    jxTil(k) = jx(k) - grav*lDep*etaDx
    jyTil(k) = jy(k) - grav*lDep*etaDy    
  enddo

PGI_ACC_TIME Output:
    obccalcjxtil  NVIDIA  devicenum=0
    time(us): 0
    801: compute region reached 1800 times
        801: kernel launched 1800 times
            grid: [203]  block: [128]
            elapsed time(us): total=4,560,568 max=3,551 min=2,502 avg=2,533
    801: data region reached 3600 times
```

- I believe the inner loop in sequential in the above form.
- However, on modifying it to use `gang, worker` for outer and `vector` for inner, huge gains could be made depending on the num of workers and vectors. The example below is the best time I could get so far, of 1.02s.
- The number of gangs is auto-calculated.
- I have documented the influence of num_workers and vector_length in the following table.

```
Code:
  !$acc parallel loop default(present) gang worker &
  !$acc   num_workers(4) vector_length(16) &
  !$acc   private(k, j, neid, etaDx, etaDy, lDep) 
  do k = 1, npt    

    !jxTil(k) = 0d0
    !jyTil(k) = 0d0
    lDep = dep(k)
    etaDx = 0d0    
    etaDy = 0d0    
    !$acc loop vector
    do j = 1, pObj(k)%nn
      neid = pObj(k)%neid(j)
      etaDx = etaDx + pObj(k)%phiDx(j)*eta(neid)
      etaDy = etaDy + pObj(k)%phiDy(j)*eta(neid)      
    enddo
    jxTil(k) = jx(k) - grav*lDep*etaDx
    jyTil(k) = jy(k) - grav*lDep*etaDy    
  enddo

PGI_ACC_TIME Output:
    obccalcjxtil  NVIDIA  devicenum=0
    time(us): 0
    801: compute region reached 1800 times
        801: kernel launched 1800 times
            grid: [6481]  block: [16x4]
            elapsed time(us): total=1,025,984 max=626 min=562 avg=569
    801: data region reached 3600 times
```


- Note the `pgaccelinfo` output in following.
- Here 'Maximum Threads per Block' is the max of [x, y, z], which is the block dimension
- Additionally 'Maximum Block Dimensions:' are the limits of [x, y, z]. 
- As per my understanding from this [link](https://forums.developer.nvidia.com/t/maximum-number-of-threads-on-thread-block/46392), both the condition should be satisfied.
    - [1024, 1, 1] Legal
    - [1024, 2, 1] Illegal
    - [1, 1024, 1] Legal
    - [1, 1, 64] Legal
    - [2, 1, 64] Legal
- Note from the table of times above, the block size of 64 for num_wrokers x vector_length gave us the fastest results.
- **Therefore, this tuning maybe specific to this GPU only!**
- Note that outer loop is `do k = 1, npt`, where npt = 25921.
- The inner loop is `do j = 1, pObj(k)%nn`, where nn can be anything maybe 30.

| Block | Grid | num_workers | vector_length | Time (s) |
| ----- | ---- | ----------- | ------------- | -------- |
| [128] | 25921 | NA | 128 | 5.50 |
| [64] | 25921 | 1 | 64 | 3.06 |
| [32x2] | 12961 | 2 | 32 | 1.87 |
| [16x2] | 12961 | 2 | 16 | 1.07 |
| [16x4] | 6481 | 4 | 16 | 1.02 |
| [8x8] | 3241 | 8 | 8 | 1.52 |

```

Device Number:                 0
Device Name:                   Quadro P620
Device Revision Number:        6.1
Global Memory Size:            2089025536
Number of Multiprocessors:     4
Concurrent Copy and Execution: Yes
Total Constant Memory:         65536
Total Shared Memory per Block: 49152
Registers per Block:           65536
Warp Size:                     32
Maximum Threads per Block:     1024
Maximum Block Dimensions:      1024, 1024, 64
Maximum Grid Dimensions:       2147483647 x 65535 x 65535
Maximum Memory Pitch:          2147483647B
Texture Alignment:             512B
Clock Rate:                    1354 MHz
Execution Timeout:             Yes
Integrated Device:             No
Can Map Host Memory:           Yes
Compute Mode:                  default
Concurrent Kernels:            Yes
ECC Enabled:                   No
Memory Clock Rate:             2505 MHz
Memory Bus Width:              128 bits
L2 Cache Size:                 524288 bytes
Max Threads Per SMP:           2048
Async Engines:                 2
Unified Addressing:            Yes
Managed Memory:                Yes
Concurrent Managed Memory:     Yes
Preemption Supported:          Yes
Cooperative Launch:            Yes
  Multi-Device:                Yes
PGI Default Target:            -ta=tesla:cc60
```

#### Lessons
Some of the lessons may be:
- In single loop divide the threads between gang and vector. Dont touch workers
- In nested loop of two levels, divide the outer loop between gand and workers, and inner loop in vectors. Ensure the vector_length is less than the number of iterations in the inner loop to reduce the overheads.

These lessons were applied also to loops inside _swe9n.f90_. Time improvements for TestCase1 are:

- swe9n.f90  799 : Solving Jx and Jy <br> 1.60s to 1.05s
- swe9n.f90  856 : Predictor Eta P Q <br> 2.30s to 1.26s
- swe9n.f90  960 : Corrector Eta P Q <br> 5.78s to 2.15s
- eqnMatrices.f90  555 : Main loop <br> 17.28s to 15.00s

So at this stage the final improvement is from **TestCase1 wallTime = 41.16s** to **TestCase1 wallTime = 29.80s**. Finally it barely crosses OpenMP.

-----------------------------------------------

<a name = 'log_swe9n_v0002_3' ></a>

### OpenACC for GWCErh2() [2020-06-17]

- It appears that update directives on line 787 and 935 (for output from _GWCErh2()_) take quite a long time. The other data updates are quite small.
- So the aim is to make _GWCErh2()_ GPU parallel.
	- This will minimise the run time for this, the slowest function
	- It will also remove the 'update' time for transfer from CPU to GPU being done twice every loop and is the slowest data update.
- Will have to use **atomic** to ensure that the matrices are updated correcttly.
	- Similar to CRITICAL in OpenMP
- **TestCase1 wallTime = 56.5s** improved to **TestCase1 wallTime = 41.16s**.
- The results match with Serial and OpenMP versions.
- The speed issue is not because of atomic. I checked.
- The subroutine _GWCErh2()_ takes 17.8s.
- If I remove the loops updating gD21, gD25, gD35, then this subroutine takes 3.03s only. So why is this simple loop taking 14.7s?!
- Similar odd behaviour is seen in _obcCalcJxTil()_.
	- Here the loop takes 4.56s.
	- But if I change `etaDx = etaDx + pObj(k)%phiDx(j)*eta(neid)` to `etaDx = etaDx + 0.01` then it takes 0.078s!!
-----------------------------------------------

<a name = 'log_swe9n_v0002_2' ></a>

### OpenACC for Time-loop Attempt 1 [2020-06-17]
- So far **TestCase1 wallTime = 56.5s** for serial version compiled on pgfortran.
- Most of the time-loop was made parallel with a lot of 'data update' statements to shift between GPU and CPU
- Shifted to `-ta=tesla:cc60,pinned`
- So far **TestCase1 wallTime = 83.67s**
- Not parallelised _GWCErh2()_. But beause of data transfers slowing down the full time loop its cost changed from 61% to 40% of total time
- Issues in making _inletBC()_ and _openBC2()_ parallel.
- **The positive thing is that results are point to point match**
- It appears that update directives on line 787 and 935 (for output from _GWCErh2()_) take quite a long time. The other data updates are quite small.
- **So the next aim is to make _GWCErh2()_ GPU parallel**.
	- This will minimise the run time for this, the slowest function
	- It will also remove the 'update' time for transfer from CPU to GPU being done twice every loop and is the slowest data update.

#### Update [2020-06-18]
- Significant gains were made just by using `gang vector` clause in `!$acc parallel loop` inside _swe9n.f90_ and one function inside _eqnMatrices.f90_
- Wall run time changed from **TestCase1 wallTime = 83.67s** to **TestCase1 wallTime = 59.67s**, a gain of 24s!
- Now the CPU run _GWCErh2()_ takes up 56% of the execution time, approx 33s.
- Also the data transfer due to CPU functioning of _GWCErh2()_ takes 3s.
- A large portion of this 36s can be cut down by making _GWCErh2()_ run on GPU.

-----------------------------------------------

<a name = 'log_swe9n_v0002_1' ></a>

### Comparison of speeds [2020-06-17]

|  |  |
| :------- | :------------ |
| **Test** | **TestCase1** |
| Mesh | sqB |
| T | 18000s |
| dt | 10s |
| OutputInterval | 1800s |
| Num of Timesteps | 1800 |

All cases run in Dhruv's system.<br>
CUDA is the GPU parallel code implemented by Dhruv.


| Compiler | Options | Parallel | CPU Cores | GPU | Time (s) | Relative |Remark |
| :------- |  :----- | :------- | :-------- | :--- | :------- | :------ | :---- |
| gfortran | -O3 | OpenMP | 4 | No | 30.11 | 0.53x |  |
| pgfortran | -fast |  | 1 | No | 56.65 | 1.00x |  |
| pgfortran |  |  | 1 | No| 335.56 | 5.92x | Prob Dhruv used in thesis as "serial" |
| pgfortran |  | CUDA | 1 | Yes | 338.32 | 5.97x | Dhruv thesis GPU reported | 
| pgfortran | -fast  | CUDA | 1 | Yes| 83.24 | 1.47x |  | 



-----------------------------------------------

## References
