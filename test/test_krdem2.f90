module krdem2_module
!! Auxiliary module for `test_krdem2`.
   use daskr_kinds, only: wp, zero, one
   implicit none

   integer, parameter :: neq = 2, nrt = 1, nrowpd = 2

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

      ! Call F(T,Y)
      call f(t, y, delta)

      ! Form G = Y' - F(T,Y)
      delta = yprime - delta

   end subroutine res

   pure subroutine f(t, y, ydot)
   !! dy/dt routine.
      real(wp), intent(in) :: t
      real(wp), intent(in) :: y(:)
      real(wp), intent(out) :: ydot(:)

      ydot(1) = y(2)
      ydot(2) = 100*(one - y(1)**2)*y(2) - y(1)

   end subroutine f

   pure subroutine jac(t, y, yprime, pd, cj, rpar, ipar)
   !! Jacobian routine.
      real(wp), intent(in) :: t
      real(wp), intent(in) :: y(neq)
      real(wp), intent(in) :: yprime(neq)
      real(wp), intent(out) :: pd(nrowpd, neq)
      real(wp), intent(in) :: cj
      real(wp), intent(in):: rpar
      integer, intent(in) :: ipar

      ! First define the Jacobian matrix for the right-hand side of the ODE:
      ! Y' = F(T,Y) , i.e. dF/dY.
      pd(1, 1) = zero
      pd(1, 2) = one
      pd(2, 1) = -200*y(1)*y(2) - one
      pd(2, 2) = 100*(one - y(1)**2)

      ! Next update the Jacobian with the right-hand side to form the DAE Jacobian:
      ! CJ*dR/dY' + dR/dY = CJ*I - dF/dY.
      pd(1, 1) = cj - pd(1, 1)
      pd(1, 2) = -pd(1, 2)
      pd(2, 1) = -pd(2, 1)
      pd(2, 2) = cj - pd(2, 2)

   end subroutine jac

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

      rval(1) = y(1)

   end subroutine rt

end module krdem2_module

program test_krdem2
!! Test program for `daskr`: intermittently stiff problem.
!!
!! The initial value problem is:
!!```
!!    dy1/dt = y2
!!    dy2/dt = 100*(1 - y1**2)*y2 - y1
!!
!!    y1(0) = 2,  y2(0) = 0,  0 <= t <= 200
!!    y1'(0) = 0, y2'(0) = -2
!!```
!! The root function is:
!!```
!!    r1(t, y, y') = y1
!!```
!! The analytical solution is not known, but the zeros of `y1` are known to 15 figures.
!!
!! If the errors are too large, or other difficulty occurs, a warning message is printed.
!! To run the demonstration problem with full printing, set `kprint=3`.

   use iso_fortran_env, only: stdout => output_unit
   use daskr_kinds, only: wp, zero, one, two
   use krdem2_module, only: res, rt, jac, neq, nrt
   implicit none

   integer, parameter :: lrw = 100, liw = 100
   integer :: idid, iout, ipar, jtype, kprint, lun, nerr, nre, nrea, nrte, nje, nst, kroot
   integer :: info(20), iwork(liw), jroot(nrt)
   real(wp) :: errt, psdum, rpar, t, tout, tzero
   real(wp) :: atol(neq), rtol(neq), rwork(lrw), y(neq), yprime(neq)

   lun = stdout
   kprint = 3
   nerr = 0
   idid = 0

   ! Set tolerance parameters and print heading.
   ! Note that INFO(2) is set to 1, indicating that RTOL and ATOL are arrays. Each entry of
   ! RTOL and ATOL must then be defined.
   info = 0
   info(2) = 1
   rtol(:) = 1e-6_wp
   atol(1) = 1e-6_wp
   atol(2) = 1e-4_wp

   if (kprint >= 2) then
      write (lun, '(/, a, /)') 'DKRDEM-2: Test Program for DASKR'
      write (lun, '(a)') 'Van Der Pol oscillator'
      write (lun, '(a)') 'Problem is dY1/dT = Y2,  dY2/dT = 100*(1-Y1**2)*Y2 - Y1'
      write (lun, '(a)') '            Y1(0) = 2,    Y2(0) = 0'
      write (lun, '(a)') 'Root function is  R(T,Y,YP) = Y1'
      write (lun, '(a, e10.1, a, 2(e10.1))') 'RTOL =', rtol(1), ' ATOL =', atol(1:2)
   end if

   ! Note: JTYPE indicates the Jacobian type:
   ! JTYPE = 1 ==> Jacobian is dense and user-supplied
   ! JTYPE = 2 ==> Jacobian is dense and computed internally

   ! Loop over JTYPE = 1, 2.  Set remaining parameters and print JTYPE.
   do jtype = 1, 2

      ! Set INFO(1) = 0 to indicate start of a new problem
      ! Set INFO(5) = 2-JTYPE to tell DDASKR the Jacobian type.
      info(1) = 0
      info(1) = 0
      info(5) = 2 - jtype
      t = zero
      y(1) = two
      y(2) = zero
      yprime(1) = zero
      yprime(2) = -two
      tout = 20.0_wp

      if (kprint > 2) then
         write (lun, '(/, 70("."), /, a, i2, /)') "Solution with JTYPE =", jtype
      end if

      ! Call DDASKR in loop over TOUT values = 20, 40, ..., 200.
      do iout = 1, 10

         do
            call ddaskr(res, neq, t, y, yprime, tout, info, rtol, atol, idid, &
                        rwork, lrw, iwork, liw, rpar, ipar, jac, psdum, rt, nrt, jroot)

            ! Print Y1 and Y2.
            if (kprint > 2) then
               write (lun, '(a, e15.7, 5x, a, e15.7, 5x, a, e15.7)') &
                  "At t =", t, "y1 =", y(1), "y2 =", y(2)
            end if

            if (idid < 0) exit

            ! If no root found, increment TOUT and loop back.
            ! If a root was found, write results and check root location.
            ! Then return to DDASKR to continue the integration.
            if (idid /= 5) then
               tout = tout + 20.0_wp
               exit
            else

               if (kprint > 2) then
                  write (lun, '(/, a, e15.7, 2x, a, i3)') "Root found at t =", t, "JROOT =", jroot(1)
               end if

               kroot = int(t/81.2_wp + 0.5_wp)
               tzero = 81.17237787055_wp + float(kroot - 1)*81.41853556212_wp
               errt = t - tzero
               if (kprint > 2) then
                  write (lun, '(a, e12.4, /)') "Error in t location of root is", errt
               end if

               if (errt >= one) then
                  nerr = nerr + 1
                  if (kprint >= 2) then
                     write (lun, '(/, a, /)') 'WARNING: Root error exceeds 1.0'
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
         write (lun, '(/, a)') 'Final statistics for this run:'
         write (lun, '(a, i5)') 'number of steps =', nst
         write (lun, '(a, i5)') 'number of Gs    =', nre
         write (lun, '(a, i5)') '(excluding Js)  =', nrea
         write (lun, '(a, i5)') 'number of Js    =', nje
         write (lun, '(a, i5)') 'number of Rs    =', nrte
      end if

   end do

   if (kprint >= 1) then
      write (lun, '(/, a, i3)') 'Number of errors encountered =', nerr
   end if

   if (nerr == 0) then
      stop ">>> Test passed. <<<"
   else
      error stop ">>> Test failed. <<<"
   end if

end program test_krdem2
