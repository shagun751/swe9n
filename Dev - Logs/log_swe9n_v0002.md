## Initial development log

1. [Comparison of speeds [2020-06-17]](#log_swe9n_v0002_1)
1. [OpenACC for Time-Loop Attempt 1 [2020-06-17][2020-06-18]](#log_swe9n_v0002_2)
1. [OpenACC for GWCErh2() [2020-06-18]](#log_swe9n_v0002_3)

### Attempting
- OpenACC directives based GPU parallelisation


### List of Work
- [ ] To be added

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
