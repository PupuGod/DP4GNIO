# DPGNIO--dynamic programming algorithm for solving generalized nearly isotonic optimization (GNIO) problem

Author: Xuyu Chen and Xudong Li 

Introduction: the DPGNIO software is implemented with C/C++, and it can be used to solve l1-GNIO problem and l2-GNIO problem. The dynamic programming principle is fully exploited in the software and it is very fast and robust. For  detaild introduction of the GNIO model, please see our paper ( ?? https://arxiv.org/pdf/2011.03305.pdf ?? ).

------------------------------------------------------------------------------------------------
The following two functions are designed to solve GNIO problems: 
1. `void l1gnio(double* data, double* w, double* l_read, double* m_read, const int n, double* solution)` in `l1gnio.cpp`
2. `void l2gnio(double* data, double* w, double* l_read, double* m_read, const int n, double* solution )` in `l2gnio.c`


In both of the functions, the inputs are:
1. `data`: an array of real numbers which needs to be optimized, and it is the same as `y` in paper. <br>
2. `w`: an array of postive real numbers, which contains the weight of loss functions. <br>
3. `l_read`: an array of non-negative real numbers with its i-th entry to be $\lambda_i$. <br>
4. `m_read`: an array of non-negative real numbers with its i-th entry to be $\mu_i$. <br>  
5. `n`: problem size, also length of `data` or `w`, `n` should be larger than 1.<br>
6. `solution`: an array to store the optimal solution of the GNIO problem. 

-----------------------------------------------------------------------------------------------

For any question, please contact chenxy18@fudan.edu.cn. 

