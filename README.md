# DP4GNIO
A dynamic programming algorithm for solving generalized nearly isotonic optimization (GNIO) problems.
It currently contains  C/C++ implementations of the dynamic programming algorithm for solving 
$\ell_1$−GNIO and $\ell_2$-GNIO problems.

Authors: Xuyu Chen and Xudong Li.





<!--
The DPGNIO softwares are C/C++ implementations of the dynamic programming algorithm (https://arxiv.org/pdf/2011.03305.pdf) designed for solving l1-GNIO or l2-GNIO problems 
-->

------------------------------------------------------------------------------------------------
The following two functions are designed to solve GNIO problems: 
1. `void l1gnio(double* data, double* w, double* l_read, double* m_read, const int n, double* solution)` in `l1gnio.cpp`
2. `void l2gnio(double* data, double* w, double* l_read, double* m_read, const int n, double* solution )` in `l2gnio.c`


In both of the functions, the inputs are:
1. `data`: an array of real input numbers, and it is the same as `y` in the paper. <br>
2. `w`: an array of positive real numbers, which contains the weight of loss functions. <br>
3. `l_read`: an array of non-negative real numbers with its i-th entry to be $\lambda_i$. <br>
4. `m_read`: an array of non-negative real numbers with its i-th entry to be $\mu_i$. <br>  
5. `n`: problem size, also the length of `data` or `w`, `n` should be larger than 1.<br>
6. `solution`: an array to store the optimal solution.

-----------------------------------------------------------------------------------------------
To use the software, please


1. input data and parameters from desired sources
2. construct an empty double array `solution` 
3. run function `l1gnio` and `l2gnio`


-------------------------------------------------------------------------------------------------
In our source files, we provide the following way to input data from the txt file:

Required data format: <br>
Assume the problem size is $n$, then the first $n$ rows should be $y_1,....,y_n$, the ($n$+1)-th to 2$n$-th rows are the weights, i.e., $w_1,...,w_n$, the (2$n$+1)-th to (3$n$-1)-th are $\lambda$s, i.e., $\lambda_1,....,\lambda_{n-1}$, and the 3$n$-th to (4$n$-2)-th are $\mu$s, i.e., $\mu_1,....,\mu_{n-1}$. For example, if $n=2$, the txt file reads:

`data_file.txt`: <br>
$y_1$ <br>
$y_2$ <br>
$w_1$ <br>
$w_2$ <br>
$\lambda_1$<br>
$\mu_1$<br>


Of course, one can always choose an alternative way to input the data as he/she wants.


------------------------------------------------------------------------------------------------------

**Citation Information**:

If you find the software DP4GNIO
useful, please cite it in your publication as follows:
*Zhensheng Yu, Xuyu Chen, and Xudong Li, A dynamic programming approach for generalized nearly isotonic optimization, arXiv:2011.03305, 2020*


For any other questions, please contact chenxy18@fudan.edu.cn. 

