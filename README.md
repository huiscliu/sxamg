# SXAMG

SXAMG is an Algebraic Multigrid Solver Library. It is written by C and is serial.

# How-to
## Configuration
The simplest way to configure is to run command:
```
./configure
```
This command will try to find optional packages from certain directories. Searching details are included by configure.in and some are explained below.

## Options
Run command:
```
./configure --help
```

Most optional packages are enabled by default. However, a package can be disabled when configuring, such as "--disable-itsol" to disable ITSOL. When a package needs an include path and a library path, they can be set by configure, such as --with-itsol-libdir=DIR and --with-itsol-incdir=DIR, which set library and include paths for ITSOL. If configure cannot find correct paths, users can help configure by using options.


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
The package will be installed to a directory. The default is /usr/local/lssp/. A different directory can be set by --prefix=DIR.

