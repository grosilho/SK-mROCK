# SK-mROCK
This code is a C++ implementation of the multirate explicit stabilized method SK-mROCK for solving stochastic differential equations (SDEs)
$$\mathrm{d}X=f(X)\mathrm{d}t+g(X)\mathrm{d}W_t,\qquad X(0)=X_0,$$
where
$$f(X)=f_F(X)+f_S(X),$$
with $f_F$ a stiff but cheap term and $f_S$ is a mildly stiff but expensive term. 

Such multiscale equations are found in chemistry, physics, and finance, for instance. We display here below an example of a multiscale population-dynamics model and of a reaction-diffusion equation solved by the SK-mROCK method.
<p align="center">
  <img src="./docs/img/population_dynamics.png" height="220"/>
  <img src="./docs/img/reac_diff.png" height="220"/>
</p>

In addition to the multirate method SK-mROCK, we implement the explicit stabilized method SK-ROCK and the standard methods Euler-Maruyama and Platen for non-multirate SDEs.


## Explicit stabilized methods for SDEs
Explicit stabilized methods for SDEs are based on explicit stabilized methods for ODEs.
These methods use an increased number of stages to increase stability, in contrast to standard methods that use more stages to increase accuracy. Due to this different strategy, the stability domain grows quadratically along the negative real axis and the methods have no step size restriction (despite being explicit).
In this code, we implement the explicit stabilized method [SK-ROCK](https://doi.org/10.1137/17M1145859) for non-multirate SDEs.

## Multirate explicit stabilized methods for SDEs
When solving a multirate equation with a non-multirate method, such as SK-ROCK, the two terms $f_F,f_S$ are evaluated together. Therefore, the number of function evaluations might depend on the severe stiffness of very few terms in $f_F$. This destroys efficiency since the expensive term $f_S$ is evaluated as many times as $f_F$. Hence, the evaluation of $f_F$ and $f_S$ must be decoupled.

The SK-mROCK method is a multirate version of SK-ROCK. It is based on and extends the methodology that we derived to design the multirate [mRCK](https://github.com/grosilho/mRKC) scheme for ODEs.

### Modified equation
Instead of solving the original multirate problem
$$\mathrm{d}X=f_F(X)\mathrm{d}t+f_S(X)\mathrm{d}t+g(X)\mathrm{d}W_t,\qquad X(0)=X_0,$$
the idea of the multirate explicit stabilized method SK-mROCK is to solve a _modified equation_
$$\mathrm{d}X_\eta=f_\eta(X_\eta)\mathrm{d}t+g_\eta(X_\eta)\mathrm{d}W_t,\qquad X_\eta(0)=X_0,$$
with the SK-ROCK method. The modified equation is an approximation to the original multirate problem, but with the advantage that it is built in a way that the stiffness of $f_\eta$ depends on the slow terms $f_S$ only, and therefore solving the modified equation is cheaper than the multirate problem. Moreover, evaluating $f_\eta$ and $g_\eta$ has a similar cost as $f_F+f_S$ and $g$, respectively.

### The auxiliary problems
The _averaged force_ $f_\eta(x)$ is defined by 
$$f_\eta(x)=\frac{1}{\eta}(u(\eta)-x),$$
where $u(\eta)$ is computed by solving an _auxiliary problem_
$$u'=f_F(u)+f_S(x) \quad t\in (0,\eta), \quad u(0)=x.$$
Since the expensive term $f_S$ is frozen, solving the auxiliary problem is comparable to evaluating $f_F+f_S$. The value of $\eta>0$ depends on the stiffness of $f_S$ and in general satisfies $\eta\ll\Delta t$, with $\Delta t$ the step size used to solve the modified equation. This deterministic auxiliary problem is solved using an explicit stabilized method for ODEs, in this case, RKC.

Since the drift term $f_F+f_S$ has been replaced by the damped averaged force $f_\eta$, then the diffusion term $g$ has to be damped as well. Indeed, the averaged force isn't able to control the original noise term $g$.

The _damped diffusion_ $g_\eta(x)$ is defined by 
$$
g_\eta(x)=\frac{1}{\eta}(v(\eta)-\bar v(\eta)),
$$
where $v(\eta)$ and $\bar v(\eta)$ are computed by solving the auxiliary problems 
$$
v'=\frac{1}{2}f_F(v)+g(x), \qquad \bar v'=\frac{1}{2}f_F(\bar v),\qquad v(0)=\bar v(0)=x.
$$
Again, computing $g_\eta$ is relatively cheap since only the cheap term $f_F$ is employed. The role of $f_F$ here is to dampen the original diffusion $g$.



## References
For more details on the mathematical background and numerical results, see the following publications:
>Abdulle, A., & Rosilho de Souza, G. (2022). _Explicit stabilized multirate method for stiff stochastic differential equations_. SIAM Journal on Scientific Computing, 44(4), A1859–A1883. [DOI:10.1137/21M1439018](https://doi.org/10.1137/21M1439018)
>
>Abdulle, A., Grote, M. J., & Rosilho de Souza, G. (2022). _Explicit stabilized multirate method for stiff differential equations_. Mathematics of Computation, 91(338), 2681–2714, [DOI:10.1090/mcom/3753](http://dx.doi.org/10.1090/mcom/3753).




## Install 
For compilation, [eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), and [GetPot](https://github.com/loumalouomega/getpot-cpp) are needed. However, they're already added as submodules. Hence, if the code has been downloaded with `git clone ...`, then compatible versions of `eigen` and `GetPot` are cloned automatically by the makefile. Hence, the following command will download the dependencies and compile the code:

> mkdir build && cd build && cmake -S .. -B . && make

Run the code from the ``build`` directory as
> ./MultirateIntegrators ARGS_LIST

For examples, see the next section.

## Run
For the full list of arguments use the `--help` option or check out the [`ARGS_LIST.md`](ARGS_LIST.md) file, here we just provide some examples.

Six problems are already hard-coded into the executable. For each one, three types of simulations can be run:

1. Compute only one sample of the solution,
2. Compute several samples of the solution and do some statistics. If a reference solution is provided then the results are compared.
3. Perform strong and weak convergence experiments.

We provide examples for each case below.

### 1. Compute one sample of the solution

- Suppose we want to solve the Population-Dynamics problem (	problem 6): `-test 6`,
- with the SK-mROCK method and step size $\Delta t=0.01$: `-rk SKmROCK -dt 1e-2`,
- compute only one sample (one Monte Carlo iteration): `-mciter 1`,
- we want to name the solution `sol`: `-ofile sol`,
- write output in binary format every 10 time steps: `-bin true -ofreq 10`.

In the `build` directory execute:

```
./MultirateIntegrators -test 6 -rk SKmROCK -dt 1e-2 -mciter 1 -ofile sol -bin true -ofreq 10
```

The output is stored in the `results` folder. Results are plotted with the `matlab/Plot_PopulationDynamics.m` script. Each one of the implemented problems has its plot script. 

Plotting is not possible when more than one sample is computed, since in that case statistics are computed.

### 2. Compute several samples
To perform many Monte Carlo iterations just change the ``-mciter`` option and run

```./MultirateIntegrators -test 6 -rk SKmROCK -dt 1e-2 -mciter 1e4 -ofile many_mciter```

We removed options `-bin true -ofreq 10` since writing the solution is valid only when computing one sample. The output is:

```
--------------- SDE Integration ---------------
Solver: SKmROCK.
Step size: 1.000000e-02.
Monte Carlo iterations: 10000.
Exact Brownian motion: yes.
Seed: -1.
Output file name: many_mciter
Output frequency: -1
Verbose: no.
-----------------------------------------------
-------------- Solver info -------------
Solver:             SKmROCK
Avg dt :            1.0000e-02
Avg f eval :        0.0000e+00
Avg ρ/Nit :         0.0000e+00
Postprocessing:     no
Time to solution:   4.7145e-01 sec
Monte Carlo iter:   10000
Successful iter:    10000
Succeed ratio:      100 %
----------------------------------------
--------- Statistics ---------
|Avg | 9.9468e-01| 8.7124e-01|
------------------------------
|Std | 4.0386e-03| 8.2949e-02|
------------------------------
```

We could compare the results against a reference solution computed with another method or with a smaller step size. To do so, we first compute a reference solution with for instance the SK-ROCK method. Hence,

- name the solution `refsol`: `-ofile refsol`,
- write all samples at the end of the simulation: `-ofreq 0`.

Execute:

```
./MultirateIntegrators -test 6 -rk SKROCK -dt 1e-3 -mciter 1e4 -ofreq 0 -ofile refsol
```

Then compute the solution with SK-mROCK with a larger step size and compare the solutions by passing the name for the reference solution:

```
./MultirateIntegrators -test 6 -rk SKmROCK -dt 1e-2 -mciter 1e4 -ofile many_mciter -refsol refsol
```

The difference between the solutions is computed and displayed:

```
------- Relative Error -------
|Avg | 2.6327e-05| 2.8743e-04|
------------------------------
|Std | 1.0889e-01| 1.1703e-02|
------------------------------

```

### 3. Perform convergence experiments
Here, we solve problem 2 and perform a convergence experiment for the SK-mROCK method. To do so, pass the `-convtest` option with the parameters for the minimal and maximal step sizes to consider. The code computes reference solutions on the fly. 

The needed options are:

- solve problem 2: `-test 2`,
- use $10^5$ Monte Carlo iterations: `-mciter 1e5`,
- perform a convergence test: `-convtest`,
- the convergence experiment is performed considering different number $N$ of steps, that go from $N=2^{7}$ to $N=2^{10}$: `-minpow 7 -maxpow 10`.

Run:
```
./MultirateIntegrators -test 2 -mciter 1e5 -rk SKmROCK -convtest -minpow 7 -maxpow 10
```

The code displays the results as follows, showing the strong and weak errors and the relative convergence rates. 

```
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

If we are interested in weak errors only, then the Brownian motion can be replaced by a discrete noise by passing the `-contW false` option.