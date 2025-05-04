!----------------------------------------------------------------------------------------------
! Adapted from original Fortran code in `original/examples/dkrdem.f`
!----------------------------------------------------------------------------------------------

module krdem1_m
!! Procedures for [[test_krdem1]].
   use daskr_kinds, only: rk, zero, one
   implicit none
   private
   
   integer, parameter :: neq = 1, nrt = 2

   public :: neq, nrt
   public :: res, rt   

contains

   pure subroutine res(t, y, yprime, cj, delta, ires, rpar, ipar)
   !! Residuals routine.
      real(rk), intent(in):: t
      real(rk), intent(in):: y(neq)
      real(rk), intent(in):: yprime(neq)
      real(rk), intent(in):: cj
      real(rk), intent(out):: delta(neq)
      integer, intent(out) :: ires
      real(rk), intent(in):: rpar
      integer, intent(in) :: ipar

      ! Check Y to make sure that it is valid input.
      ! If Y is less than or equal to zero, this is invalid input.
      if (y(1) <= zero) then
         ires = -1
         return
      else
         call f(t, y, delta)
         delta = yprime - delta
      end if

   end subroutine res
   
   pure subroutine f(t, y, yprime)
   !! dy1/dt routine.
      real(rk), intent(in) :: t
      real(rk), intent(in) :: y(:)
      real(rk), intent(out) :: yprime(:)

      yprime(1) = ((2*log(y(1)) + 8.0_rk)/t - 5.0_rk)*y(1)

   end subroutine f

   pure subroutine rt(neq, t, y, yprime, nrt, rval, rpar, ipar)
     !! Roots routine.
      integer, intent(in) :: neq
      real(rk), intent(in) :: t
      real(rk), intent(in) :: y(neq)
      real(rk), intent(in) :: yprime(neq)
      integer, intent(in) :: nrt
      real(rk), intent(out) :: rval(nrt)
      real(rk), intent(in) :: rpar
      integer, intent(in) :: ipar

      rval(1) = yprime(1)
      rval(2) = log(y(1)) - 2.2491_rk

   end subroutine rt

end module krdem1_m

program test_krdem1
!! Test program for [[daskr]]: nonstiff problem.
!!
!! The initial value problem is:
!!
!! $$\begin{aligned}
!!  y'(t) &= \left(\frac{2 \ln(y) + 8}{t} - 5\right)y \\
!!  y(1)  &= 1, \quad 1 \le t \le 6
!! \end{aligned}$$
!! 
!! The solution is:
!!
!! $$ y(t) = \exp(-t^2 + 5t - 4), \quad y'(1) = 3 $$
!!
!! The two root functions are:
!!   
!! $$\begin{aligned}
!! r_1(t,y,y') &= y'              \quad &&(\text{with root at } t = 2.5), \\
!! r_2(t,y,y') &= \ln(y) - 2.2491 \quad &&(\text{with roots at } t = 2.47 \text{ and } 2.53)
!! \end{aligned}$$
!!
!! If the errors are too large, or other difficulty occurs, a warning message is printed.
!! To run the demonstration problem with full printing, set `kprint=3`.

   use iso_fortran_env, only: stdout => output_unit
   use daskr_kinds, only: rk, zero, one, two
   use krdem1_m, only: res, rt, neq, nrt
   implicit none

   integer, parameter :: lrwork = 76, liwork = 41 ! @note: to be replaced by formula or alloc
   integer :: idid, iout, ipar, jdum, jtype, kprint, lout, nerr, nre, nrea, nrte, nje, nst
   integer :: info(20), iwork(liwork), jroot(2)
   real(rk) :: er, ero, errt, psdum, rpar, t, tout, yt
   real(rk) :: atol(neq), rtol(neq), rwork(lrwork), y(neq), yprime(neq)

  ! Set report options
   lout = stdout
   kprint = 3
   
   ! Initialize variables and set tolerance parameters.
   nerr = 0
   idid = 0
   info = 0
   rtol = zero
   atol = 1e-6_rk

   ! Set INFO(11) = 1 if DASKR is to compute the initial YPRIME, and generate an initial guess
   ! for YPRIME.  Otherwise, set INFO(11) = 0 and supply the correct initial value for YPRIME.
   info(11) = 0
   y(1) = one
   yprime(1) = 3.0_rk

   ! Note: JTYPE indicates the Jacobian type:
   ! JTYPE = 1 ==> Jacobian is dense and user-supplied
   ! JTYPE = 2 ==> Jacobian is dense and computed internally
   jtype = 2
   info(5) = 2 - jtype
   if (kprint >= 2) then
      write (lout, '(/, a, /)') 'DKRDEM-1: Test Program for DASKR'
      write (lout, '(a)') 'Problem is  dY/dT = ((2*LOG(Y)+8)/T - 5)*Y,  Y(1) = 1'
      write (lout, '(a)') 'Solution is  Y(T) = EXP(-T**2 + 5*T - 4)'
      write (lout, '(a)') 'Root functions are:'
      write (lout, '(a)') 'R1 = dY/dT  (root at T = 2.5)'
      write (lout, '(a)') 'R2 = LOG(Y) - 2.2491  (roots at T = 2.47 and T = 2.53)'
      write (lout, '(a, e10.1, a, e10.1, a, i3, /)') &
         'RTOL =', rtol(1), ' ATOL =', atol(1), ' JTYPE =', jtype
   end if

   ! Call DASKR in loop over TOUT values = 2, 3, 4, 5, 6.
   t = one
   tout = two
   ero = zero
   do iout = 1, 5
      do
         call daskr(res, neq, t, y, yprime, tout, info, rtol, atol, idid, &
                    rwork, lrwork, iwork, liwork, rpar, ipar, jdum, psdum, rt, nrt, jroot)

         ! Print Y and error in Y, and print warning if error too large.
         yt = exp(-t**2 + 5*t - 4.0_rk)
         er = y(1) - yt

         if (kprint > 2) then
            write (lout, '(a, e15.7, 5x, a, e15.7, 5x, a, e12.4)') &
               'At t =', t, ' y =', y(1), ' error =', er
         end if

         if (idid < 0) exit

         er = abs(er)/atol(1)
         ero = max(ero, er)
         if (er >= 1e3_rk) then
            nerr = nerr + 1
            if (kprint >= 2) then
               write (lout, '(/, a, /)') 'WARNING: Error exceeds 1e3*tolerance'
            end if
         end if

         ! If no root found, increment TOUT and loop back.
         ! If a root was found, write results and check root location.
         ! Then return to DASKR to continue the integration.
         if (idid /= 5) then
            tout = tout + one
            exit
         else
            if (kprint > 2) then
               write (lout, '(/, 4x, a, e15.7, 5x, a, 2i5)') &
                  'Root found at t =', t, ' JROOT =', jroot(1:2)
            end if

            if (jroot(1) /= 0) errt = t - 2.5_rk
            if (jroot(2) /= 0 .and. t <= 2.5_rk) errt = t - 2.47_rk
            if (jroot(2) /= 0 .and. t > 2.5_rk) errt = t - 2.53_rk
            if (kprint > 2) then
               write (lout, '(4x, a, e12.4, /)') 'Error in t location of root is', errt
            end if

            if (abs(errt) >= 1e-3_rk) then
               nerr = nerr + 1
               if (kprint >= 2) then
                  write (lout, '(/, a, /)') 'WARNING: Root error exceeds 1e-3'
               end if
            end if

         end if
      
      end do
      
      if (idid < 0) exit
   
   end do

   ! Problem complete. Print final statistics.
   if (idid < 0) nerr = nerr + 1
   nst = iwork(11)
   nre = iwork(12)
   nje = iwork(13)
   nrte = iwork(36)
   nrea = nre
   if (jtype == 2) nre = nre + neq*nje

   if (kprint >= 2) then
      write (lout, '(/, a)') 'Final statistics for this run:'
      write (lout, '(a, i5)') 'number of steps =', nst
      write (lout, '(a, i5)') 'number of Gs    =', nre
      write (lout, '(a, i5)') '(excluding Js)  =', nrea
      write (lout, '(a, i5)') 'number of Js    =', nje
      write (lout, '(a, i5)') 'number of Rs    =', nrte
      write (lout, '(a, e10.2)') 'error overrun   =', ero
   end if

   if (kprint >= 1) then
      write (lout, '(/, a, i3)') 'Number of errors encountered =', nerr
   end if
   
   if (nerr == 0) then
      stop '>>> Test passed. <<<'
   else
      error stop '>>> Test failed. <<<'
   end if

end program test_krdem1
