# daskr

Modernized version of DASKR, a differential-algebraic system solver with rootfinding.

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)

## Description

tbd

## History

The first version of the library, named DASSL, was originally released in 1982 [1]. This was followed by the release of DASPK in 1994 [2]. The last official version, named DASKR, was released in 1998 [3]. 

`daskr` is a modernization of the DASKR code [4], intented to make the library easier to use and maintain. The main changes include:

* [ ] Conversion from fixed-form (`.f`) to free-form (`.f90`).
* [ ] Conversion from upper case to lower case.
* [ ] Modularization.
* [ ] Removal of `DATA` statements, labeled do loops, and (most) `goto`s.
* [ ] Addition of `intent(in/out)` to all procedures.
* [ ] Addition of explicit interfaces to BLAS routines.
* [ ] Implementation of a C API.
* [ ] Automatic code documentation with FORD.

**References**

[1] L. Petzold, "A Description of DASSL: A Differential/Algebraic System Solver, 1982.

[2] Brown, Peter N., Alan C. Hindmarsh, and Linda R. Petzold. "Using Krylov methods in the solution of large-scale differential-algebraic systems." SIAM Journal on Scientific Computing 15.6 (1994): 1467-1488. https://doi.org/10.1137/0915088

[3] Brown, Peter N., Alan C. Hindmarsh, and Linda R. Petzold. "Consistent initial condition calculation for differential-algebraic systems." SIAM Journal on Scientific Computing 19.5 (1998): 1495-1512.
https://doi.org/10.1137/S1064827595289996

[4] Original source code from [Netlib](https://www.netlib.org/odrpack/).


## Build instructions

### With fpm

The easiest way to build/test the code and run the examples is by means of [`fpm`](https://fpm.fortran-lang.org/).

To build the library, do:

```sh
fpm build --profile release
```

To run the tests, do:

```sh
fpm test --profile release
```

To run the provided examples, do:

```sh
fpm run --example "example_name" --profile release
```

### With meson

First, setup the build:

```sh
meson setup builddir -Dbuild_tests=true
```

To build the libraries, do:

```sh
meson compile -C builddir
```

To run the tests, do:

```sh
meson test -C builddir
```

## Licence

* The original `daskr` code is covered by this [license](./original/LICENSE).
* Modifications introduced in this project are covered under the MIT license.