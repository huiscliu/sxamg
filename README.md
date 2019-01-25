# SXAMG

SXAMG is an Algebraic Multigrid Solver Library written by C. The initial code was from FASP
(http://fasp.sourceforge.net/), a linear solver package that implemented Krylov solvers, preconditioners and
AMG solvers. The SXAMG was extracted from FASP to serve as a standalone AMG solver.

# How-to
## Configuration
The simplest way to configure is to run command:
```
./configure
```
This command will detect system and set default options.

## Options
Run command:
```
./configure --help
```

## Compilation
After configuration, a simple **make** command can compile the package:
```
make
```

## Installation
Run command:
```
make install
```
The package will be installed to a directory. The default is /usr/local/sxamg/. A different directory can be set by --prefix=DIR.

# Demo

The following codes show how to use SXAMG as a solver.

```
{
    SX_AMG_PARS pars;
    SX_MAT A;
    SX_VEC b, x;
    int verb = 2;
    SX_RTN rtn;
    
    /* default pars */
    sx_amg_pars_init(&pars);

    /* redefine parameters */
    pars.maxit = 1000;
    pars.verb = 2;
    
    /* set A, b, initial x */
    ....
    
    /* solve the system */
    rtn = sx_solver_amg(&A, &x, &b, &pars);
    
    /* free memory */
    sx_mat_destroy(&A);
    sx_vec_destroy(&b);
    sx_vec_destroy(&x);
}
```
