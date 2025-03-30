module krdem1_module
!! Auxiliary module for `example_krdem1`.
   use daskr_kinds, only: wp, zero, one
   implicit none
   
   ! @note: not so happy with this, but it is required for the time being
   integer, parameter :: neq = 1, nrt = 2

contains

   pure subroutine res(t, y, yprime, cj, delta, ires, rpar, ipar)
   !! Residuals routine.
      real(wp), intent(in):: t
      real(wp), intent(in):: y(neq)
      real(wp), intent(in):: yprime(neq)
      real(wp), intent(in):: cj
      real(wp), intent(out):: delta(neq)
      integer, intent(inout) :: ires
      real(wp), intent(in):: rpar
      integer, intent(in) :: ipar

      ! Check Y to make sure that it is valid input.
      ! If Y is less than or equal to zero, this is invalid input.
      if (y(1) <= zero) then
         ires = -1
         return
      else

         ! Call F to obtain F(T,Y)
         call f(t, y, delta)

         ! Form G = Y' - F(T,Y)
         delta = yprime - delta

      end if

   end subroutine res

   pure subroutine f(t, y, ydot)
   !! dy1/dt routine.
      real(wp), intent(in) :: t
      real(wp), intent(in) :: y(:)
      real(wp), intent(out) :: ydot(:)

      ydot(1) = ((2*log(y(1)) + 8.0_wp)/t - 5.0_wp)*y(1)

   end subroutine f

   pure subroutine rt(neq, t, y, yp, nrt, rval, rpar, ipar)
     !! Roots routine.
      integer, intent(in) :: neq
      real(wp), intent(in) :: t
      real(wp), intent(in) :: y(neq)
      real(wp), intent(in) :: yp(neq)
      integer, intent(in) :: nrt
      real(wp), intent(out) :: rval(nrt)
      real(wp), intent(in) :: rpar
      integer, intent(in) :: ipar

      rval(1) = yp(1)
      rval(2) = log(y(1)) - 2.2491_wp

   end subroutine rt

end module krdem1_module

program test_krdem1
!! Test program for `daskr`: nonstiff problem.
!!
!! The initial value problem is:
!!```
!!    dy/dt = ((2*log(y) + 8)/t - 5)*y,  y(1) = 1,  1 <= t <= 6
!!``` 
!! The solution is:
!!```
!!    y(t) = exp(-t**2 + 5*t - 4), y'(1) = 3
!!```
!! The two root functions are:
!!```
!!    r1(t,y,y') = dy/dt  (with root at t = 2.5),
!!    r2(t,y,y') = log(y) - 2.2491  (with roots at t = 2.47 and 2.53)
!!```
!!
!! If the errors are too large, or other difficulty occurs, a warning message is printed.
!! To run the demonstration problem with full printing, set `kprint=3`.

   use iso_fortran_env, only: stdout => output_unit
   use daskr_kinds, only: wp, zero, one, two
   use krdem1_module, only: res, rt, neq, nrt
   implicit none

   integer, parameter :: lrw = 76, liw = 41
   integer :: idid, iout, ipar, jdum, jtype, kprint, lun, nerr, nre, nrea, nrte, nje, nst
   integer :: info(20), iwork(liw), jroot(2)
   real(wp) :: er, ero, errt, psdum, rpar, t, tout, yt
   real(wp) :: atol(neq), rtol(neq), rwork(lrw), y(neq), yprime(neq)

   ! Set all input parameters and print heading.
   info = 0
   lun = stdout
   kprint = 3
   nerr = 0
   idid = 0

   y(1) = one
   t = one
   tout = two
   rtol = zero
   atol = 1e-6_wp

   ! Set INFO(11) = 1 if DASKR is to compute the initial YPRIME, and generate an initial guess
   ! for YPRIME.  Otherwise, set INFO(11) = 0 and supply the correct initial value for YPRIME.
   info(11) = 0
   yprime(1) = 3.0_wp

   ! Note: JTYPE indicates the Jacobian type:
   ! JTYPE = 1 ==> Jacobian is dense and user-supplied
   ! JTYPE = 2 ==> Jacobian is dense and computed internally
   jtype = 2
   info(5) = 2 - jtype
   if (kprint >= 2) then
      write (lun, '(a, /)') 'DKRDEM-1: Test Program for DASKR'
      write (lun, '(a)') 'Problem is  dY/dT = ((2*LOG(Y)+8)/T - 5)*Y,  Y(1) = 1'
      write (lun, '(a)') 'Solution is  Y(T) = EXP(-T**2 + 5*T - 4)'
      write (lun, '(a)') 'Root functions are:'
      write (lun, '(a)') 'R1 = dY/dT  (root at T = 2.5)'
      write (lun, '(a)') 'R2 = LOG(Y) - 2.2491  (roots at T = 2.47 and T = 2.53)'
      write (lun, '(a, e10.1, a, e10.1, a, i3, /)') &
         'RTOL =', rtol(1), ' ATOL =', atol(1), ' JTYPE =', jtype
   end if

   ! Call DASKR in loop over TOUT values = 2, 3, 4, 5, 6.
   ero = zero
   do iout = 1, 5
      do
         call ddaskr(res, neq, t, y, yprime, tout, info, rtol, atol, idid, &
                     rwork, lrw, iwork, liw, rpar, ipar, jdum, psdum, rt, nrt, jroot)

         ! Print Y and error in Y, and print warning if error too large.
         yt = exp(-t**2 + 5*t - 4.0_wp)
         er = y(1) - yt

         if (kprint > 2) then
            write (lun, '(a, e15.7, 5x, a, e15.7, 5x, a, e12.4)') &
               'At t =', t, ' y =', y(1), ' error =', er
         end if

         if (idid < 0) exit

         er = abs(er)/atol(1)
         ero = max(ero, er)
         if ((er >= 1e3_wp) .and. (kprint >= 2)) then
            nerr = nerr + 1
            write (lun, '(/, a, /)') 'WARNING: Error exceeds 1e3*tolerance'
         end if

         ! If no root found, increment TOUT and loop back.
         ! If a root was found, write results and check root location.
         ! Then return to DASKR to continue the integration.
         if (idid /= 5) then
            tout = tout + one
            exit
         else
            if (kprint > 2) then
               write (lun, '(/, 4x, a, e15.7, 5x, a, 2i5)') &
                  "Root found at t =", t, " JROOT =", jroot(1:2)
            end if

            if (jroot(1) /= 0) errt = t - 2.5_wp
            if (jroot(2) /= 0 .and. t <= 2.5_wp) errt = t - 2.47_wp
            if (jroot(2) /= 0 .and. t > 2.5_wp) errt = t - 2.53_wp
            if (kprint > 2) then
               write (lun, '(4x, a, e12.4, /)') 'Error in t location of root is', errt
            end if

            if ((abs(errt) >= 1e-3_wp) .and. (kprint >= 2)) then
               nerr = nerr + 1
               write (lun, '(/, a, /)') 'WARNING: Root error exceeds 1e-3'
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

   if (kprint > 2) then
      write (lun, '(/, a)') "Final statistics for this run:"
      write (lun, '(a, i5)') "number of steps =", nst
      write (lun, '(a, i5)') "number of Gs    =", nre
      write (lun, '(a, i5)') "(excluding Js)  =", nrea
      write (lun, '(a, i5)') "number of Js    =", nje
      write (lun, '(a, i5)') "number of Rs    =", nrte
      write (lun, '(a, e10.2)') "error overrun   =", ero
   end if

   if (kprint >= 2) then
      write (lun, '(/, 80("-"), /, a, i3)') 'Number of errors encountered =', nerr
   end if
   
   if (nerr == 0) then
      stop "Test passed."
   else
      error stop "Test failed."
   end if

end program test_krdem1
