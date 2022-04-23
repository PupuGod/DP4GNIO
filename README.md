# DP4GNIO
A dynamic programming algorithm for solving generalized nearly isotonic optimization (GNIO) problems.
It currently contains  C/C++ implementations (with a Python Wrapper) of the dynamic programming algorithm for solving $\ell_1$−GNIO and $\ell_2$-GNIO problems.


Authors: Xuyu Chen and Xudong Li.





<!--
The DPGNIO softwares are C/C++ implementations of the dynamic programming algorithm (https://arxiv.org/pdf/2011.03305.pdf) designed for solving l1-GNIO or l2-GNIO problems 
-->

------------------------------------------------------------------------------------------------
The following two functions are designed to solve GNIO problems: 
1. `l1gnio(data, w, lbd, mu, check = True)`
2. `l2gnio(data, w, lbd, mu, check = True)`

In both of the functions, the parameters are:

1. `data`: a ndarray of float;
2.  `w`: a ndarray of positive float numbers;
3. `lbd`: a ndarray of nonnegative float numbers, with  i-th entry to be $\lambda_i$;
4. `mu`: a ndarray of nonnegative float numbers, with  i-th entry to be $\mu_i$;
5. `check`: a boolean value. If it is set to be `true`, the software verifies whether the inputs are legal before the computations. It might be time-consuming for large-scale problems.

The returns are:
`result`: A `2d-tuple` with first element to be the solution (also a ndarray), and 
the second element to be the objective value.

------------------------------------
To use the software, please
1. find the correct version corresponds to you OS (Windows, MacOS, Linux).
2. make sure the python libraries `numpy` and `ctypes` are avaiable.
3. excute `runfile_OS.py` (OS = Win, Mac, Liunx).



------------------------------------------------------------------------------------------------------

**Citation Information**:

If you find the software DP4GNIO
useful, please cite it in your publication as follows:
*Zhensheng Yu, Xuyu Chen, and Xudong Li, A dynamic programming approach for generalized nearly isotonic optimization, arXiv:2011.03305, 2020*


For any other questions, please contact chenxy18@fudan.edu.cn. 

