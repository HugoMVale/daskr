!----------------------------------------------------------------------------------------------
! Adapted from original Fortran code in `original/examples/dkrdem.f`
!----------------------------------------------------------------------------------------------

module krdem2_m
!! Procedures for [[test_krdem2]].

   use daskr_kinds, only: rk, zero, one
   implicit none
   private

   integer, parameter :: neq = 2, nrt = 1, nrowpd = 2

   public :: neq, nrt, nrowpd
   public :: res, jac, rt

contains

   pure subroutine res(t, y, ydot, cj, delta, ires, rpar, ipar)
   !! Residuals routine.

      real(rk), intent(in):: t
      real(rk), intent(in):: y(neq)
      real(rk), intent(in):: ydot(neq)
      real(rk), intent(in):: cj
      real(rk), intent(out):: delta(neq)
      integer, intent(out) :: ires
      real(rk), intent(inout):: rpar
      integer, intent(inout) :: ipar

      call f(t, y, delta)
      delta = ydot - delta

   end subroutine res
   
   pure subroutine f(t, y, ydot)
   !! dy/dt routine.

      real(rk), intent(in) :: t
      real(rk), intent(in) :: y(:)
      real(rk), intent(out) :: ydot(:)

      ydot(1) = y(2)
      ydot(2) = 100*(one - y(1)**2)*y(2) - y(1)

   end subroutine f

   pure subroutine jac(t, y, ydot, pd, cj, rpar, ipar)
   !! Jacobian routine.

      real(rk), intent(in) :: t
      real(rk), intent(in) :: y(neq)
      real(rk), intent(in) :: ydot(neq)
      real(rk), intent(out) :: pd(nrowpd, neq)
      real(rk), intent(in) :: cj
      real(rk), intent(in):: rpar
      integer, intent(in) :: ipar

      ! First define the Jacobian matrix for the right-hand side of the ODE:
      ! Y' = F(T,Y), i.e. dF/dY.
      pd(1, 1) = zero
      pd(1, 2) = one
      pd(2, 1) = -200*y(1)*y(2) - one
      pd(2, 2) = 100*(one - y(1)**2)

      ! Next update the Jacobian with the right-hand side to form the DAE Jacobian:
      ! cj*dG/dY' + dG/dY = cj*I - dF/dY.
      pd(1, 1) = cj - pd(1, 1)
      pd(1, 2) = -pd(1, 2)
      pd(2, 1) = -pd(2, 1)
      pd(2, 2) = cj - pd(2, 2)

   end subroutine jac

   pure subroutine rt(neq, t, y, ydot, nrt, rval, rpar, ipar)
     !! Roots routine.
   
      integer, intent(in) :: neq
      real(rk), intent(in) :: t
      real(rk), intent(in) :: y(neq)
      real(rk), intent(in) :: ydot(neq)
      integer, intent(in) :: nrt
      real(rk), intent(out) :: rval(nrt)
      real(rk), intent(in) :: rpar
      integer, intent(in) :: ipar

      rval(1) = y(1)

   end subroutine rt

end module krdem2_m

program test_krdem2
!! Test program for [[daskr]]: intermittently stiff problem.
!!
!! The initial value problem is:
!!
!! $$\begin{aligned}
!! \dot{y}_1(t) &= y_2 \\
!! \dot{y}_2(t) &= 100(1 - y_1^2)y_2 - y_1
!! \end{aligned}$$
!! 
!! with initial conditions:
!!      
!! $$\begin{aligned}
!! y_1(0)  &= 2, \quad y_2(0)  = 0, \quad 0 \le t \le 200 \\
!! \dot{y}_1(0) &= 0, \quad \dot{y}_2'(0) = -2
!! \end{aligned}$$
!!
!! The root function is:
!!
!! $$ r_1(t, y, \dot{y}) = y_1 $$
!!
!! The analytical solution is not known, but the zeros of \(y_1\) are known to 15 figures.
!!
!! If the errors are too large, or other difficulty occurs, a warning message is printed.
!! To run the demonstration problem with full printing, set `kprint = 3`.

   use iso_fortran_env, only: stdout => output_unit
   use daskr_kinds, only: rk, zero, one, two
   use krdem2_m, only: res, rt, jac, neq, nrt
   implicit none

   integer, parameter :: lrwork = 100, liwork = 100 ! @note: to be replaced by formula or alloc
   integer :: idid, iout, ipar, jtype, kprint, lout, nerr, nre, nrea, nrte, nje, nst, kroot
   integer :: info(20), iwork(liwork), jroot(nrt)
   real(rk) :: errt, psdum, rpar, t, tout, tzero
   real(rk) :: atol(neq), rtol(neq), rwork(lrwork), y(neq), ydot(neq)

   ! Set report options
   lout = stdout
   kprint = 3

   ! Initialize variables and set tolerance parameters.
   ! Note that INFO(2) is set to 1, indicating that RTOL and ATOL are arrays. Each entry of
   ! RTOL and ATOL must then be defined.
   nerr = 0
   idid = 0
   info = 0
   info(2) = 1
   rtol(:) = 1e-6_rk
   atol(1) = 1e-6_rk
   atol(2) = 1e-4_rk

   if (kprint >= 2) then
      write (lout, '(/, a, /)') "DKRDEM-2: Test Program for DASKR"
      write (lout, '(a)') "Van Der Pol oscillator"
      write (lout, '(a)') "Problem is dy1/dt = y2,  dy2/dt = 100*(1-y1**2)*y2 - y1"
      write (lout, '(a)') "            y1(0) = 2,    y2(0) = 0"
      write (lout, '(a)') "Root function is  r(t,y,y') = y1"
      write (lout, '(a, e10.1, a, 2(e10.1))') 'RTOL =', rtol(1), ' ATOL =', atol(1:2)
   end if

   ! Note: JTYPE indicates the Jacobian type:
   ! JTYPE = 1 ==> Jacobian is dense and user-supplied
   ! JTYPE = 2 ==> Jacobian is dense and computed internally

   ! Loop over JTYPE = 1, 2.  Set remaining parameters and print JTYPE.
   do jtype = 1, 2

      ! Set INFO(1) = 0 to indicate start of a new problem
      ! Set INFO(5) = 2 - JTYPE to tell DASKR the Jacobian type.
      info(1) = 0
      info(1) = 0
      info(5) = 2 - jtype
      t = zero
      y(1) = two
      y(2) = zero
      ydot(1) = zero
      ydot(2) = -two
      tout = 20.0_rk

      if (kprint > 2) then
         write (lout, '(/, 70("."), /, a, i2, /)') 'Solution with JTYPE =', jtype
      end if

      ! Call DASKR in loop over TOUT values = 20, 40, ..., 200.
      do iout = 1, 10

         do
            call daskr(res, neq, t, y, ydot, tout, info, rtol, atol, idid, &
                       rwork, lrwork, iwork, liwork, rpar, ipar, jac, psdum, rt, nrt, jroot)

            ! Print Y1 and Y2.
            if (kprint > 2) then
               write (lout, '(a, e15.7, 5x, a, e15.7, 5x, a, e15.7)') &
                  'At t =', t, 'y1 =', y(1), 'y2 =', y(2)
            end if

            if (idid < 0) exit

            ! If no root found, increment TOUT and loop back.
            ! If a root was found, write results and check root location.
            ! Then return to DASKR to continue the integration.
            if (idid /= 5) then
               tout = tout + 20.0_rk
               exit
            else

               if (kprint > 2) then
                  write (lout, '(/, a, e15.7, 2x, a, i3)') &
                     'Root found at t =', t, 'JROOT =', jroot(1)
               end if

               kroot = int(t/81.2_rk + 0.5_rk)
               tzero = 81.17237787055_rk + (kroot - 1)*81.41853556212_rk
               errt = t - tzero
               if (kprint > 2) then
                  write (lout, '(a, e12.4, /)') 'Error in t location of root is', errt
               end if

               if (errt >= one) then
                  nerr = nerr + 1
                  if (kprint >= 2) then
                     write (lout, '(/, a, /)') 'WARNING: Root error exceeds 1.0'
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
      end if

   end do

   if (kprint >= 1) then
      write (lout, '(/, a, i3)') 'Number of errors encountered =', nerr
   end if

   if (nerr == 0) then
      stop '>>> Test passed. <<<'
   else
      error stop '>>> Test failed. <<<'
   end if

end program test_krdem2
