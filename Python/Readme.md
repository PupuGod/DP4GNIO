# How to use the Python version of DPGNIO


Please collect all the files (l1gnio.py, l2gnio.py, runfile_gnio.py) in the same path and, <br>


(1) open runfile_gnio.py <br>
(2) read in data and set optimization parameters <br>
(3) to use l1-GNIO model, please replace line 52 in runfile_gnio.py with 'f = l1.func()' <br>
(4) to use l2-GNIO model, please replace line 52 in runfile_gnio.py with 'f = l2.func()' <br>
(5) run runfile_gnio.py <br>


Moreover, to change the inputing data, we shall introduce meanings of the variables: <br>


(1) data: a n dimensional list of floating numbers, which contains the series to be optimized <br>
(2) w: a n dimenisonal positive np.array which contains weight in the loss function <br>
(3) lbd: a n-1 dimensional non-negative np.array which controls (increasing) monotonicity <br>
(4) mu: a n-1 dimensional non-negative np.array which controls (decreasing) monotonicity <br>
(5) n : the dimensional (or length) of the data



Directly run the runfile_gnio.py without any modification, one may obtain a example of appling l2-GNIO model on random data.



The users can also read data from excels or other soruces into the memory if necessary, but this is not concluded in our runfile.


Any questions, please contact chenxy18@fudan.edu.cn


