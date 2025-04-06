# daskr

Modernized version of DASKR, a differential-algebraic system solver with rootfinding.

[![Test](https://github.com/HugoMVale/daskr/actions/workflows/test.yml/badge.svg)](https://github.com/HugoMVale/daskr/actions)
[![codecov](https://codecov.io/gh/HugoMVale/daskr/graph/badge.svg?token=AgjzeQ1qFL)](https://codecov.io/gh/HugoMVale/daskr)
[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)

## Status

This is a work in progress. First, I'm modernizing the examples and converting some to unit tests. Then, I'll begin work on the library code itself.

## Description

`daskr` is a library for solving systems of differential-algebraic equations of the form:

$$ G(t, y, y') = 0 $$
$$ y(t_0) = y_0 $$
$$ y'(t_0) = y'_0 $$

where $G$, $y$, and $y'$ are $N$-dimensional vectors. The linear systems which arise at each time step can be solved with dense or banded _direct_ methods (Gaussian elimination with partial pivoting) or with _iterative_ Krylov methods (preconditioned [GMRES]). Additionally, it includes the ability to find the roots of a given set of functions while carrying out the integration.

[GMRES]: https://en.wikipedia.org/wiki/Generalized_minimal_residual_method

## History

The first version of the library, named DASSL [1], solved the linear systems arising from the implicit time integration methods at each time step using direct methods. DASPK [2,3] extended the capabilities of DASSL to include iterative methods, which can be significantly more efficient, especially for large-scale problems. Furthermore, DASPK added the ability to initialize $y'_0$ in case it is not known. Lastly, DASKR [4] included the ability to find the roots of a given set of functions while integrating the DAE system.

| Version | Date written  | Last update | Direct solver | Iterative solver | Root finding |    Standard   |
|:-------:|:-------------:|:-----------:|:-------------:|:----------------:|:------------:|:-------------:|
|  daskr  |      2025     |      --     |       ☑       |         ☑       |       ☑      | Fortran 2018  |
|  DASKR  |      2002     |     2011    |       ☑       |         ☑       |       ☑      |   FORTRAN 77  |
|  DASPK  |      1989     |     2000    |       ☑       |         ☑       |       ☐      |   FORTRAN 77  |
|  DASSL  |      1983     |     2000    |       ☑       |         ☐       |       ☐      |   FORTRAN 77  |

`daskr` is a modernization of the DASKR code [4], intended to make the library easier to use and maintain. The main changes include:

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

[4] Original source code from [Netlib](https://www.netlib.org/ode/).

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