# SK-mROCK
A C++ implementation of the SK-mROCK multirate numerical integrator for stiff stochastic differential equations.


Repository under construction :)

```
./MultirateIntegrators -test 6 -mc iter 1e4 -rk SKmROCK -ofreq 0 -ofile refsol
--------------- SDE Integration ---------------
Solver: SKmROCK.
Step size: 1.000000e-02.
Monte Carlo iterations: 1000.
Exact Brownian motion: yes.
Seed: -1.
Output file name: refsol
Output frequency: -1
Verbose: no.
-----------------------------------------------
-------------- Solver info -------------
Solver:             SKmROCK
Avg dt :            1.0000e-02
Avg f eval :        0.0000e+00
Avg ρ/Nit :         0.0000e+00
Postprocessing:     no
Time to solution:   4.9680e-02 sec
Monte Carlo iter:   1000
Successful iter:    1000
Succeed ratio:      100 %
----------------------------------------
--------- Statistics ---------
|Avg | 9.9444e-01| 8.6592e-01|
------------------------------
|Std | 4.0543e-03| 8.3347e-02|
------------------------------
```

```
./MultirateIntegrators -test 6 -mc iter 1e4 -rk SKmROCK -refsol refsol
--------------- SDE Integration ---------------
Solver: SKmROCK.
Step size: 1.000000e-02.
Monte Carlo iterations: 1000.
Exact Brownian motion: yes.
Seed: -1.
Output file name: sol
Output frequency: -1
Verbose: no.
-----------------------------------------------

----------------- SOLUTION DATA -----------------

-------------- Solver info -------------
Solver:             SKmROCK
Avg dt :            1.0000e-02
Avg f eval :        0.0000e+00
Avg ρ/Nit :         0.0000e+00
Postprocessing:     no
Time to solution:   4.7467e-02 sec
Monte Carlo iter:   1000
Successful iter:    1000
Succeed ratio:      100 %
----------------------------------------
--------- Statistics ---------
|Avg | 9.9456e-01| 8.6842e-01|
------------------------------
|Std | 4.2373e-03| 8.6982e-02|
------------------------------

------------ REFERENCE SOLUTION DATA ------------

-------------- Solver info -------------
Solver:             SKmROCK
Avg dt :            1.0000e-02
Avg f eval :        0.0000e+00
Avg ρ/Nit :         0.0000e+00
Postprocessing:     no
Time to solution:   4.9680e-02 sec
Monte Carlo iter:   1000
Successful iter:    1000
Succeed ratio:      100 %
----------------------------------------
--------- Statistics ---------
|Avg | 9.9444e-01| 8.6592e-01|
------------------------------
|Std | 4.0543e-03| 8.3347e-02|
------------------------------

-------------------- ERRORS --------------------

------- Relative Error -------
|Avg | 1.1631e-04| 2.8974e-03|
------------------------------
|Std | 4.5133e-02| 4.3612e-02|
------------------------------

--- Densisty distance area ---
|nbin|         10|         10|
------------------------------
|SDX | 1.1284e-01| 1.1284e-01|
------------------------------
|SDrX| 1.1284e-01| 1.1284e-01|
------------------------------
|DDA | 0.0000e+00| 1.4000e-02|
------------------------------
```

```
./MultirateIntegrators -test 2 -mciter 1e5 -rk SKmROCK -convtest -minpow 7 -maxpow 10 -safe_add 1
--------------- SDE Integration ---------------
Solver: SKmROCK.
Step size: 1.000000e-02.
Monte Carlo iterations: 100000.
Exact Brownian motion: yes.
Seed: -1.
Convergence test parameters: min_pow = 7, max_pow = 10.
-----------------------------------------------
Elapsed time: 172.7
----------------- Convergence Results ------------------
|S Err | 1.9509e+02| 1.3861e+02| 9.2612e+01| 6.0840e+01|
--------------------------------------------------------
|S Rate| 0.0000e+00| 4.9314e-01| 5.8175e-01| 6.0617e-01|
--------------------------------------------------------
|W Err | 9.8073e-02| 4.5637e-02| 2.0932e-02| 8.8923e-03|
--------------------------------------------------------
|W Rate| 0.0000e+00| 1.1037e+00| 1.1245e+00| 1.2351e+00|
--------------------------------------------------------
```