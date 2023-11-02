```
./MultirateIntegrators --help
This code implements the SK-ROCK and SK-mROCK explicit stabilized methods for stochastic differential equations.
Standard methods as Euler-Maruyama and Platen scheme are implemented for comparison.
Different problems are hardcoded in the executable. Look into src/DSdeProblems.cpp.
Run the code from the ./build folder as ./MultirateIntegrators OPTIONS, where OPTIONS is a combination of the below.
Results are stored in the ./results folder.

The following options are available:
--- General options:
    -test       : a number in 1-6 specifying the problem that we want to solve. Default: 1
                  The list of problems is given below.
    -mciter     : Number of Monte Carlo iterations. Default: 1e3
    -ofile      : name of output file. Default: sol
    -rk         : The name of the numerical integrator to use. Default: SKROCK
                  The list of integrators is given below.
    -dt         : Time step size. Default: 1e-2
    -rfreq      : Frequency at which the spectral radii are re-estimated. Default: 5
    -safe_add   : Use some additional stages in explicit stabilized integrators. Default: 0
    -contW      : Generate Brownian motion from a continuous uniform distribution (true)
                  or a discrete distribution (false). Default: true
    -seed       : Fix the seed. Default: -1
                  If -1 then a random seed is chosen.
    -verb       : Enables or disables verbosity. Default: true
--- Options for when doing many Monte Carlo iterations
    -ofreq      : Output frequency. Default: -1
                  - If 0 writes solution only at the end of simulation.
                    In general, used to generate a reference solution.
                  - If -1 never writes the solution.
                    In general, used to compute a solution and just compare it against a reference solution.
    -refsol     : name of the reference solution (if available). Default: ""
                  If available, errors are computed at the end of the simulation.
    -convtest   : If provided, performs a time convergence test, i.e. runs several simulations and checks errors. Default: false
                  If a reference solution is not provided, it is computed on the fly.
    -maxpow     : Minimal step size used for the convergence test is dt=tend/2^maxpow. Default: 6
    -minpow     : Maximal step size used for the convergence test is dt=tend/2^minpow. Default: 3
---- Options for when doing one Monte Carlo iteration (-mciter 1)
     Here are the options to write the solution on file, for later display with matlab scripts.
    -bin        : Writes solution in binary format. Default: false
    -matlab     : Writes solution in matlab format. Default: false
    -ofreq      : Output frequency. Default: -1
                  - If >0 writes solution every ofreq time steps,
--- List of problems:
 This is the list of problems hardcoded in src/DSdeProblems.cpp. You choose them via the -test option. To change a parameter go to src/DSdeProblems.cpp, change it and recompile.
 All problems are solvable with both a standard and a multirate solver.
     1 : The Dahlquist test equation.
     2 : A scalar non stiff non linear problem.
     3 : A problem with many diffusion terms.
     4 : A stochastic Reaction-Diffusion problem.
     5 : A heat equation with non uniform mesh.
     6 : A population dynamics problem.
--- List of numerical integrators:
 Those marked with (m) are multirate integrators.
     EM     : Euler Maruyama.
     Platen : Platen Scheme.
     SKROCK : SK-ROCK method.
     SKmROCK: SK-mROCK multirate method.
   
```