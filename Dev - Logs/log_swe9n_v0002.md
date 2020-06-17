## Initial development log

1. [Comparison of speeds [2020-06-17]](#log_swe9n_v0002_1)

### Attempting
- OpenACC directives based GPU parallelisation


### List of Work
- [ ] To be added

-----------------------------------------------


<a name = 'log_swe9n_v0002_1' ></a>

### Comparison of speeds [2020-06-17]
Test case _sqB_ run for <br> 
T = 18000s <br> 
dt = 10s <br>
OutputInterval = 1800s <br>
Number of Timesteps = 1800

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
