!----------------------------------------------------------------------------------------------
! Adapted from original Fortran code in `original/examples/dheat.f`
!----------------------------------------------------------------------------------------------

module heat_m
!! Auxiliary module for [[example_heat]].
   use daskr_kinds, only: rk, zero, one
   implicit none

contains

   pure subroutine uinit(u, uprime, rpar, ipar)
   !! This routine computes and loads the vector of initial values.
   !! The initial `u` values are given by the polynomial `u = 16x(1-x)y(1-y)`.
   !! The initial `uprime` values are set to zero. ([[daskr]] corrects these during the first 
   !! time step.)
      real(rk), intent(out) :: u(:)
      real(rk), intent(out) :: uprime(:)
      real(rk), intent(in) :: rpar(2)
      integer, intent(in) :: ipar(4)

      integer :: i, ioff, j, k, m
      real(rk) :: dx, xj, yk

      m = ipar(4)
      dx = rpar(1)

      do k = 0, m + 1
         yk = k*dx
         ioff = (m + 2)*k
         do j = 0, m + 1
            xj = j*dx
            i = ioff + j + 1
            u(i) = 16*xj*(one - xj)*yk*(one - yk)
         end do
      end do

      uprime = zero

   end subroutine uinit

   pure subroutine res(t, u, uprime, cj, delta, ires, rpar, ipar)
   !! User-supplied residuals subroutine.
   !! It computes the residuals for the 2D discretized heat equation, with zero boundary values.
      real(rk), intent(in) :: t
      real(rk), intent(in) :: u(*)
      real(rk), intent(in) :: uprime(*)
      real(rk), intent(in) :: cj
      real(rk), intent(out) :: delta(*)
      integer, intent(inout) :: ires
      real(rk), intent(in) :: rpar(2)
      integer, intent(in) :: ipar(4)

      integer :: i, ioff, j, k, m, m2, neq
      real(rk) :: coeff, temx, temy

      ! Set problem constants using IPAR and RPAR.
      neq = ipar(3)
      m = ipar(4)
      coeff = rpar(2)
      m2 = m + 2

      ! Load U into DELTA, in order to set boundary values.
      delta(1:neq) = u(1:neq)

      ! Loop over interior points, and load residual values.
      do k = 1, m
         ioff = m2*k
         do j = 1, m
            i = ioff + j + 1
            temx = u(i - 1) + u(i + 1)
            temy = u(i - m2) + u(i + m2)
            delta(i) = uprime(i) - (temx + temy - 4*u(i))*coeff
         end do
      end do

   end subroutine res

   pure subroutine rt(neq, t, u, uprime, nrt, rval, rpar, ipar)
   !! Roots routine.
      integer, intent(in) :: neq
      real(rk), intent(in) :: t
      real(rk), intent(in) :: u(neq)
      real(rk), intent(in) :: uprime(neq)
      integer, intent(in) :: nrt
      real(rk), intent(out) :: rval(nrt)
      real(rk), intent(in) :: rpar(2)
      integer, intent(in) :: ipar(4)

      real(rk) :: umax

      umax = maxval(u)
      rval(1) = umax - 0.1_rk
      rval(2) = umax - 0.01_rk

   end subroutine rt

end module heat_m

program example_heat
!! Example program for [[daskr]]:
!! DAE system derived from the discretized heat equation on a square; solution  using the 
!! Krylov option with banded preconditioner.
!!
!! This program solves a DAE system that arises from the heat equation,
!! 
!! $$ \frac{\partial u}{\partial t} = 
!! \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} $$
!!       
!! posed on the 2D unit square \(0 \le x,y \le y\) with zero Dirichlet boundary conditions. 
!! An \((M+2)^2\) mesh is set on the square, with uniform spacing \(\Delta x = \Delta y = 1/(M+1)\).
!! The spatial deriviatives are represented by standard central finite difference approximations.
!! At each interior point of the mesh, the discretized PDE becomes an ODE for the discrete value
!! of \(u\). At each point on the boundary, we pose the equation \(u = 0\). The discrete values 
!! of \(u\) form a vector \(U\), ordered first by \(x\), then by \(y\). The result is a DAE 
!! system \(G(t,U,U') = 0\) of size \(\mathrm{NEQ} = (M+2)^2\).
!!   
!! The initial conditions are posed as:
!!
!! $$ u(x,y,0) = 16x(1-x)y(1-y) $$
!!   
!! The problem is solved by [[daskr]] on the time interval \(0 \le t \le 10.24\).
!!
!! The root functions are:
!!   
!! $$\begin{aligned}
!! r_1(u) &= \max(u) - 0.1 \\
!! r_2(u) &= \max(u) - 0.01
!! \end{aligned}$$
!!
!! The Krylov linear system solution method, with preconditioning, is selected. The 
!! preconditioner is a band matrix with half-bandwidths equal to 1, i.e. a tridiagonal matrix.
!! (The true half-bandwidths are equal to \(M+2\)). This corresponds to ignoring the y-direction
!! coupling in the ODEs, for purposes of preconditioning. The extra iterations resulting from
!! this approximation are offset by the lower storage and linear system solution costs for a
!! tridiagonal matrix.
!!
!! The routines [[DBANJA]] and [[DBANPS]] that generate and solve the banded preconditioner are
!! provided in a separate file for general use.
!!
!! The output times are \(t = 0.01 \times 2^n, (n = 0,..., 10)\). The maximum of \(|u|\) over
!! the mesh and various performance statistics are printed.
!!
!! For details and test results on this problem, see the reference:
!!
!! * Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold, "Using Krylov Methods in the
!!   Solution of Large-Scale Differential-Algebraic Systems", SIAM J. Sci. Comput., 15 (1994),
!!   pp. 1467-1488.
!! 
!! @note
!! [[example_heatilu]] solves the same problem, but using sparse incomplete LU factorization.

   use iso_fortran_env, only: stdout => output_unit
   use daskr_kinds, only: rk, one, zero
   use heat_m
   implicit none

   integer, parameter :: maxm = 10, maxm2 = maxm + 2, mxneq = maxm2**2, &
                         lenrw = 107 + 18*mxneq, leniw = 40, leniwp = mxneq, &
                         lenwp = 4*mxneq + 2*((mxneq/3) + 1)
   integer :: idid, iout, lenpd, liw, liwp, lrw, lwp, m, mband, ml, msave, mu, &
              ncfl, ncfn, neq, nli, nni, nout, npe, nps, nqu, nre, nrt, nrte, nst
   integer :: iwork(leniw + leniwp), info(20), ipar(4), jroot(2)

   real(rk) :: atol, avdim, coeff, dx, hu, rtol, t, tout, umax
   real(rk) :: u(mxneq), uprime(mxneq), rwork(lenrw + lenwp), rpar(2)

   external :: dbanja, dbanps

   ! Here set parameters for the problem being solved. Use RPAR and IPAR to communicate these
   ! to the other routines.
   m = maxm
   dx = one/(m + one)
   neq = (m + 2)*(m + 2)
   coeff = one/dx**2

   ipar(3) = neq
   ipar(4) = m
   rpar(1) = dx
   rpar(2) = coeff

   ! Here set NRT = number of root functions
   nrt = 2

   ! Here set the half-bandwidths and load them into IPAR for use by the preconditioner
   ! routines.
   ml = 1
   mu = 1
   ipar(1) = ml
   ipar(2) = mu

   ! Here set the lengths of the preconditioner work arrays WP and IWP, load them into IWORK,
   ! and set the total lengths of WORK and IWORK.
   lenpd = (2*ml + mu + 1)*neq
   mband = ml + mu + 1
   msave = (neq/mband) + 1
   lwp = lenpd + 2*msave
   liwp = neq
   iwork(27) = lwp
   iwork(28) = liwp
   lrw = lenrw + lwp
   liw = leniw + liwp

   ! Call subroutine UINIT to initialize U and UPRIME.
   call uinit(u, uprime, rpar, ipar)

   ! Here we set up the INFO array, which describes the various options in the way we want
   ! DASKR to solve the problem. In this case, we select the iterative preconditioned Krylov 
   ! method, and we supply the band preconditioner routines DBANJA/DBANPS.
   !
   ! We first initialize the entire INFO array to zero, then set select entries to nonzero
   ! values for desired solution options.
   !
   ! To select the Krylov iterative method for the linear systems, we set INFO(12) = 1.
   !
   ! Since we are using a preconditioner that involves approximate Jacobian elements requiring
   ! preprocessing, we have a JAC routine, namely subroutine DBANJA, and we must set INFO(15) = 1
   ! to indicate this to DASKR.
   !
   ! No other entries of INFO need to be changed for this example.
   info = 0
   info(12) = 1
   info(15) = 1

   ! Here we set tolerances for DASKR to indicate how much accuracy we want in the solution,
   ! in the sense of local error control.
   ! For this example, we ask for pure absolute error control with a tolerance of 1e-5.
   rtol = zero
   atol = 1.0e-5_rk

   ! Here we generate a heading with important parameter values.
   ! stdout is the unit number of the output device.
   write (stdout, '(a, /)') &
      'DHEAT: Heat Equation Example Program for DASKR'
   write (stdout, '(a, i3, a, i4)') &
      'M+2 by M+2 mesh, M =', m, ',  System size NEQ =', neq
   write (stdout, '(a)') &
      'Root functions are: R1 = max(u) - 0.1 and R2 = max(u) - 0.01'
   write (stdout, '(a, i3, a)') &
      'Linear solver method flag INFO(12) =', info(12), '    (0 = direct, 1 = Krylov)'
   write (stdout, '(a, i3, a, i3)') &
      'Preconditioner is a banded approximation with ML =', ml, '  MU =', mu
   write (stdout, '(a, e10.1, a, e10.1, /)') &
      'Tolerances are RTOL =', rtol, '   ATOL =', atol
   write (stdout, '("t", 12x, "UMAX", 9x, "NQ", 6x, "H", 8x, "STEPS", 5x, "NNI", 5x, "NLI")')

   !-------------------------------------------------------------------------------------------
   ! Now we solve the problem.
   !
   ! DASKR will be called to compute 11 intermediate solutions from tout=0.01 to tout=10.24
   ! by powers of 2.
   !
   ! We pass to DASKR the names DBANJA and DBANPS for the JAC and PSOL routines to do the
   ! preconditioning.
   !
   ! At each output time, we compute and print the max-norm of the solution (which should decay
   ! exponentially in t).  We also print some relevant statistics -- the current method order
   ! and step size, the number of time steps so far, and the numbers of nonlinear and
   ! linear iterations so far.
   !
   ! If a root was found, we flag this, and return to the DASKR call.
   !
   ! If DASKR failed in any way (IDID < 0) we print a message and stop the integration.
   !-------------------------------------------------------------------------------------------
   
   t = zero
   tout = 0.01_rk
   nout = 11
   do iout = 1, nout

      do
         call daskr(res, neq, t, u, uprime, tout, info, rtol, atol, idid, &
                    rwork, lrw, iwork, liw, rpar, ipar, dbanja, dbanps, rt, nrt, jroot)

         umax = maxval(abs(u))
         
         hu = rwork(7)
         nqu = iwork(8)
         nst = iwork(11)
         nni = iwork(19)
         nli = iwork(20)
         write (stdout, '(e11.5, e12.4, i5, e14.3, i7, i8, i8)') t, umax, nqu, hu, nst, nni, nli

         if (idid == 5) then
            write (stdout, '(15x, a, 2i3)') '*****   Root found, JROOT =', jroot(1:2)
         else
            exit
         end if
      
      end do

      if (idid < 0) then
         write (stdout, '(/, a, e12.4,/)') ' Final time reached =', t
         exit
      end if

      tout = 2*tout

   end do

   ! Here we display some final statistics for the problem.
   ! The ratio NLI/NNI is the average dimension of the Krylov subspace involved in the Krylov 
   ! linear iterative method.
   nst = iwork(11)
   npe = iwork(13)
   nre = iwork(12) + npe*mband
   liw = iwork(17)
   lrw = iwork(18)
   nni = iwork(19)
   nli = iwork(20)
   nps = iwork(21)
   if (nni /= 0) avdim = real(nli)/real(nni)
   ncfn = iwork(15)
   ncfl = iwork(16)
   nrte = iwork(36)

   write (stdout, '(/, a)') &
      'Final statistics for this run:'
   write (stdout, '(a, i5, a, i4)') &
      'RWORK size =', lrw, '   IWORK size =', liw
   write (stdout, '(a, i5)') &
      'Number of time steps ................ =', nst
   write (stdout, '(a, i5)') &
      'Number of residual evaluations ...... =', nre
   write (stdout, '(a, i5)') &
      'Number of root function evaluations . =', nrte
   write (stdout, '(a, i5)') &
      'Number of preconditioner evaluations  =', npe
   write (stdout, '(a, i5)') &
      'Number of preconditioner solves ..... =', nps
   write (stdout, '(a, i5)') &
      'Number of nonlinear iterations ...... =', nni
   write (stdout, '(a, i5)') &
      'Number of linear iterations ......... =', nli
   write (stdout, '(a, f8.4)') &
      'Average Krylov subspace dimension ... =', avdim
   write (stdout, '(x, i5, x, a, i5, x, a, /)') &
      ncfn, 'nonlinear conv. failures,', ncfl, 'linear conv. failures'

end program example_heat
