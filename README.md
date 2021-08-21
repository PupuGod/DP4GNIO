# DPGNIO--dynamic programming algorithm for solving generalized nearly isotonic optimization (GNIO) problem

Author: Xuyu Chen, Xudong Li and Zhensheng Yu.

Introduction: the DPGNIO software is implemented with C/C++, and it can be used to solve l1-GNIO problem and l2-GNIO problem. The dynamic programming principle is fully exploited in the software and it is very fast and robust. For  detaild introduction of the GNIO model, please see our paper ( https://arxiv.org/pdf/2011.03305.pdf ).

------------------------------------------------------------------------------------------------

How to use:
1.  build up a project with `l1gnio.cpp` (if you want to use l1-GNIO model) or `l2gnio.c` (if you want to use l2-GNIO model). 
2.  prepare data and parameters including weight, lambda and mu in a txt file with required format.
3. check the `filepath` and problem size `n` in the source files (`l1gnio.cpp` or `l2gnio.c`).
4. start the computation, and the optimal solution is stored in the variable `solution`.


--------------------------------------------------------------------------------------------------
Required data format:
in our software, we shall input path of a txt file and the problem size, the txt file should obey a required data format. Assume the problem size is $n$, then the first n rows should be the data, i.e., $y_1,....,y_n$, the ($n$+1)-th to 2$n$-th rows are the weights, i.e., $w_1,...,w_n$, the (2$n$+1)-th to (3$n$-1)-th should be $\lambda$s, i.e., $\lambda_1,....,\lambda_{n-1}$, and the 3$n$-th to (4$n$-2)-th should be $\mu$s, i.e., $\mu_1,....,\mu_{n-1}$. For example, if $n=2$, the txt file should be:

`data_file.txt`: <br>
$y_1$ <br>
$y_2$ <br>
$w_1$ <br>
$w_2$ <br>
$\lambda_1$<br>
$\mu_1$<br>

-------------------------------------------------------------------------------------------


Notice: for huge data set, the user may encounter `StackOverflow Error`, the reason is that default stack memory is two small, and this error can be avoided by increasing the stack memory. 



-----------------------------------------------------------------------------------------------

For any question, please contact (chenxy18@fudan.edu.cn)
