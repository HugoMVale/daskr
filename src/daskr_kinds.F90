module daskr_kinds
!! Real kinds and common numeric constants.
    use, intrinsic :: iso_fortran_env, only: real32, real64
    implicit none
    
    integer, parameter :: sp = real32
    integer, parameter :: dp = real64

#ifdef REAL32
    integer, parameter :: rk = sp
#elif REAL64
    integer, parameter :: rk = dp
#else
    integer, parameter :: rk = dp
#endif

    real(rk), parameter :: zero = 0.0_rk
    real(rk), parameter :: half = 0.5_rk
    real(rk), parameter :: one = 1.0_rk
    real(rk), parameter :: two = 2.0_rk
    real(rk), parameter :: ten = 10.0_rk
    real(rk), parameter :: pi = 4*atan(one)
    
end module daskr_kinds