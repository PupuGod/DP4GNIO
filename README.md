# DPGNIO--dynamic programming algorithm for solving generalized nearly isotonic optimization (GNIO) problem

Author: Xuyu Chen and Xudong Li 

Introduction: the DPGNIO softwares are C/C++ implementations of the dynamic programming algorithm (https://arxiv.org/pdf/2011.03305.pdf) designed for solvinng l1-GNIO or l2-GNIO problems 

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
To use the softwares, please


1. input data and parameters from desired sources
2. construct an empty double array `solution` 
3. run function `l1gnio` and `l2gnio`


-------------------------------------------------------------------------------------------------
In our source files, we provide a way to input data from txt files in the required format:

Required data format: <br>
Assume the problem size is $n$, then the first n rows should be the data, i.e., $y_1,....,y_n$, the ($n$+1)-th to 2$n$-th rows are the weights, i.e., $w_1,...,w_n$, the (2$n$+1)-th to (3$n$-1)-th should be $\lambda$s, i.e., $\lambda_1,....,\lambda_{n-1}$, and the 3$n$-th to (4$n$-2)-th should be $\mu$s, i.e., $\mu_1,....,\mu_{n-1}$. For example, if $n=2$, the txt file should be:

`data_file.txt`: <br>
$y_1$ <br>
$y_2$ <br>
$w_1$ <br>
$w_2$ <br>
$\lambda_1$<br>
$\mu_1$<br>


One can always choose an alternative way to input the data as he/she wants
------------------------------------------------------------------------------------------------------



For any question, please contact chenxy18@fudan.edu.cn. 

