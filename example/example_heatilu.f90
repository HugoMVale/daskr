module heatilu_module
        !! Auxiliary module for example_heatilu
   use daskr_kinds, only: wp, zero, one
   implicit none

contains

   pure subroutine uinit(u, uprime, rpar, ipar)
   !! This routine computes and loads the vector of initial values.
   !! The initial U values are given by the polynomial u = 16x(1-x)y(1-y).
   !! The initial UPRIME values are set to zero. (DASKR corrects these during the first time step.)
      real(wp), intent(out) :: u(:)
      real(wp), intent(out) :: uprime(:)
      real(wp), intent(in) :: rpar(4)
      integer, intent(in) :: ipar(34)

      integer :: i, ioff, j, k, m
      real(wp) :: dx, xj, yk

      m = ipar(34)
      dx = rpar(3)

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

   pure subroutine resh(t, u, uprime, cj, delta, ires, rpar, ipar)
   !! This is the user-supplied RES subroutine for this example.
   !! It computes the residuals for the 2-D discretized heat equation, with zero boundary values.
      real(wp), intent(in) :: t
      real(wp), intent(in) :: u(*)
      real(wp), intent(in) :: uprime(*)
      real(wp), intent(in) :: cj
      real(wp), intent(out) :: delta(*)
      integer, intent(in) :: ires
      real(wp), intent(in) :: rpar(4)
      integer, intent(in) :: ipar(34)

      integer :: i, ioff, j, k, m, m2, neq
      real(wp) :: coeff, temx, temy

      ! Set problem constants using IPAR and RPAR.
      neq = ipar(33)
      m = ipar(34)
      coeff = rpar(4)
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

   end subroutine resh

   pure subroutine rtheat(neq, t, u, up, nrt, rval, rpar, ipar)
   !! This routine finds the max of U, and sets RVAL(1) = max(u) - 0.1, RVAL(2) = max(u) - 0.01.
      integer, intent(in) :: neq
      real(wp), intent(in) :: t
      real(wp), intent(in) :: u(neq)
      real(wp), intent(in) :: up(neq)
      integer, intent(in) :: nrt
      real(wp), intent(out) :: rval(nrt)
      real(wp), intent(in) :: rpar(4)
      integer, intent(in) :: ipar(34)

      integer :: i
      real(wp) :: umax

      umax = zero
      do i = 1, neq
         umax = max(umax, u(i))
      end do
      rval(1) = umax - 0.1_wp
      rval(2) = umax - 0.01_wp

   end subroutine rtheat

end module heatilu_module

program example_heatilu
!! Example program for `daskr`:
!! DAE system derived from the discretized heat equation on a square.
!!
!! This program solves a DAE system that arises from the heat equation,
!! $$
!! \frac{\partial u}{\partial t} =
!! \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}
!! $$
!! posed on the 2-D unit square with zero Dirichlet boundary conditions.
!! An M+2 by M+2 mesh is set on the square, with uniform spacing 1/(M+1).
!! The spatial deriviatives are represented by standard central finite
!! difference approximations.  At each interior point of the mesh,
!! the discretized PDE becomes an ODE for the discrete value of u.
!! At each point on the boundary, we pose the equation u = 0.  The
!! discrete values of u form a vector U, ordered first by x, then by y.
!! The result is a DAE system G(t,U,U') = 0 of size NEQ = (M+2)*(M+2).
!!
!! Initial conditions are posed as u = 16x(1-x)y(1-y) at t = 0.
!! The problem is solved by DDASKR on the time interval 0 <= t <= 10.24.
!!
!! The root functions are R1(U) = max(u) - 0.1, R2(U) = max(u) - 0.01.
!!
!! The Krylov linear system solution method, with preconditioning, is
!! selected.  The preconditioner is a sparse matrix with half-bandwidths
!! equal to 1, i.e. a tridiagonal matrix.  (The true half-bandwidths
!! are equal to M+2.)  This corresponds to ignoring the y-direction
!! coupling in the ODEs, for purposes of preconditioning.  The extra
!! iterations resulting from this approximation are offset by the lower
!! storage and linear system solution costs for a tridiagonal matrix.
!!
!! The routines DJACILU and DPSOLILU that generate and solve the sparse
!! preconditioner are provided in a separate file for general use.
!!
!! The output times are t = .01 * 2**n (n = 0,...,10).  The maximum of
!! abs(u) over the mesh, and various performance statistics, are printed.
!!
!! For details and test results on this problem, see the reference:
!!
!!   Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
!!   Using Krylov Methods in the Solution of Large-Scale Differential-
!!   Algebraic Systems, SIAM J. Sci. Comput., 15 (1994), pp. 1467-1488.

! Here are necessary declarations.  The dimension statements use a
! maximum value for the mesh parameter M.
   use iso_fortran_env, only: stdout => output_unit
   use daskr_kinds, only: wp, one, zero
   use heatilu_module
   implicit none

   integer, parameter :: lenpfac = 5, lenplufac = 5, ipremeth = 1, lfililut = 5, &
                         ireorder = 1, isrnorm = 1, normtype = 2, jacout = 0, &
                         jscalcol = 1, maxm = 10, maxm2 = maxm + 2, mxneq = maxm2**2, &
                         lenwp = 2*lenpfac*mxneq + lenplufac*mxneq + isrnorm*mxneq &
                         + 2*(mxneq + 1), &
                         leniwp = 4*(mxneq + 1) + 3*lenpfac*mxneq + 2*lenplufac*mxneq &
                         + ireorder*2*mxneq + (ipremeth - 1)*2*mxneq, &
                         lenrw = 107 + 18*mxneq, leniw = 40

   real(wp), parameter :: permtol = 0.01_wp, tolilut = 0.001_wp

   integer :: i, idid, ierr, iout, liw, liwpmin, lrw, lwpmin, m, mband, ml, mu, ncfl, ncfn, &
              neq, nli, nni, nout, npe, nps, nqu, nre, nrt, nrte, nst
   integer :: info(20), iwork(leniw + leniwp), ipar(34), jroot(2)

   real(wp) :: atol, avdim, coeff, dx, hu, rtol, t, tout, umax
   real(wp) :: rpar(4), rwork(lenrw + lenwp), u(mxneq), uprime(mxneq)

   external :: djacilu, dpsolilu

   ! Open matrix output file if JACOUT==1
   if (jacout == 1) then
      ipar(29) = 1
      open (unit=1, file='Heat_Test_Matrix.dat', status='unknown')
   end if

   ! Here set parameters for the problem being solved. Use RPAR and IPAR to communicate these
   ! to the other routines.
   m = maxm
   dx = one/(m + one)
   neq = (m + 2)*(m + 2)
   coeff = one/dx**2
   ipar(33) = neq
   ipar(34) = m
   rpar(3) = dx
   rpar(4) = coeff

   ! Set NRT = number of root functions.
   nrt = 2

   ! Here set the lengths of the preconditioner work arrays WP and IWP, load them into IWORK,
   !and set the total lengths of WORK and IWORK.
   iwork(27) = lenwp
   iwork(28) = leniwp
   lrw = lenrw + lenwp
   liw = leniw + leniwp

   ! Load values into IPAR and RPAR for sparse preconditioner.
   ml = 1
   mu = 1
   ipar(1) = ml
   ipar(2) = mu
   ipar(3) = lenpfac
   ipar(4) = lenplufac
   ipar(5) = ipremeth
   ipar(6) = lfililut
   ipar(7) = ireorder
   ipar(8) = isrnorm
   ipar(9) = normtype
   ipar(10) = jacout
   ipar(11) = jscalcol
   ipar(30) = 0
   rpar(1) = tolilut
   rpar(2) = permtol

   ! Check IPAR, RPAR, LENWP and LENIWP for illegal entries and long enough work array lengths.
   call dspsetup(neq, lenwp, leniwp, rpar, ipar, ierr, lwpmin, liwpmin)
   if (ierr /= 0) then
      write (stdout, '(a, i5)') ' Error return from DSPSETUP: IERR = ', ierr
      if (lwpmin > lenwp) then
         write (stdout, *) ' More WP work array length needed'
      end if
      if (liwpmin > leniwp) then
         write (stdout, *) ' More IWP work array length needed'
      end if
      stop
   end if

   ! Call subroutine UINIT to initialize U and UPRIME.
   call uinit(u, uprime, rpar, ipar)

   !-----------------------------------------------------------------------
   ! Here we set up the INFO array, which describes the various options
   ! in the way we want DDASKR to solve the problem.
   ! In this case, we select the iterative preconditioned Krylov method,
   ! and we supply the sparse preconditioner routines DJACILU/DPSOLILU.
   !
   ! We first initialize the entire INFO array to zero, then set select
   ! entries to nonzero values for desired solution options.
   !
   ! To select the Krylov iterative method for the linear systems,
   ! we set INFO(12) = 1.
   !
   ! Since we are using a preconditioner that involves approximate
   ! Jacobian elements requiring preprocessing, we have a JAC routine,
   ! namely subroutine DJACILU, and we must set INFO(15) = 1 to indicate
   ! this to DDASKR.
   !
   ! No other entries of INFO need to be changed for this example.
   !-----------------------------------------------------------------------
   info = 0
   info(12) = 1
   info(15) = 1

   ! Here we set tolerances for DDASKR to indicate how much accuracy we want in the solution, 
   ! in the sense of local error control.
   ! For this example, we ask for pure absolute error control with a tolerance of 1e-5.
   rtol = zero
   atol = 1.0e-5_wp

   ! Here we generate a heading with important parameter values.
   write (stdout, '(5x, a, //)') 'HEATILU: Heat Equation Example Program for DDASKR'
   write (stdout, '(5x, a, i3, a, i4)') 'M+2 by M+2 mesh, M =', m, ', System size NEQ =', neq
   write (stdout, '(5x, a)') 'Root functions are: R1 = max(u) - 0.1 and R2 = max(u) - 0.01'
   write (stdout, '(5x, a, i3, a)') 'Linear solver method flag INFO(12) =', info(12), ' (0 = direct, 1 = Krylov)'
   write (stdout, '(5x, a, i3, a, i3)') 'Preconditioner is a sparse approximation with ML =', ml, ' MU =', mu
   write (stdout, '(5x, a, i2, a)') 'Incomplete factorization option =', ipremeth, ' (1 = ILUT, 2 = ILUTP)'
   write (stdout, '(5x, a, e10.1, a, e10.1, //)') 'Tolerances are RTOL =', rtol, ' ATOL =', atol
   write (stdout, "(5x, 't', 12x, 'UMAX', 8x, 'NQ', 8x, 'H', 8x, 'STEPS', 5x, 'NNI', 5x, 'NLI')")

   !-------------------------------------------------------------------------------------------
   ! Now we solve the problem.
   !
   ! DASKR will be called to compute 11 intermediate solutions from tout=0.01 to tout=10.24
   ! by powers of 2.
   !
   ! We pass to DASKR the names DJACILU and DPSOLILU for the JAC and PSOL routines to do the
   ! preconditioning.
   !
   ! At each output time, we compute and print the max-norm of the solution (which should
   ! decay exponentially in t). We also print some relevant statistics -- the current method
   ! order and step size, the number of time steps so far, and the numbers of nonlinear and
   ! linear iterations so far.
   !
   ! If a root was found, we flag this, and return to the DASKR call.
   !
   ! If DASKR failed in any way (IDID < 0) we print a message and stop the integration.
   !-------------------------------------------------------------------------------------------

   t = zero
   tout = 0.01_wp
   nout = 11
   do iout = 1, nout
      do
         call ddaskr(resh, neq, t, u, uprime, tout, info, rtol, atol, &
                     idid, rwork, lrw, iwork, liw, rpar, ipar, djacilu, dpsolilu, &
                     rtheat, nrt, jroot)

         umax = zero
         do i = 1, neq
            umax = max(umax, abs(u(i)))
         end do

         hu = rwork(7)
         nqu = iwork(8)
         nst = iwork(11)
         nni = iwork(19)
         nli = iwork(20)
         write (stdout, '(E15.5, E12.4, I5, E14.3, I7, I9, I8)') t, umax, nqu, hu, nst, nni, nli

         if (idid == 5) then
            write (stdout, '(20X, A, 2I3)') '*****   Root found, JROOT =', jroot(1), jroot(2)
         else
            exit
         end if
      end do

      if (idid < 0) then
         write (stdout, '(/, A, E12.4, /)') ' Final time reached =', t
         exit
      end if

      tout = 2*tout

   end do

   ! Here we display some final statistics for the problem.
   ! The ratio of NLI to NNI is the average dimension of the Krylov subspace involved in the 
   ! Krylov linear iterative method.
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

   write (stdout, '(//, 5x, a)') 'Final statistics for this run:'
   write (stdout, '(5x, a, i5, a, i4)') 'RWORK size =', lrw, ' IWORK size =', liw
   write (stdout, '(5x, a, i5)') 'Number of time steps ................ =', nst
   write (stdout, '(5x, a, i5)') 'Number of residual evaluations ...... =', nre
   write (stdout, '(5x, a, i5)') 'Number of res. evals. for precond.... =', ipar(30)
   write (stdout, '(5x, a, i5)') 'Number of root function evaluations . =', nrte
   write (stdout, '(5x, a, i5)') 'Number of preconditioner evaluations  =', npe
   write (stdout, '(5x, a, i5)') 'Number of preconditioner solves ..... =', nps
   write (stdout, '(5x, a, i5)') 'Number of nonlinear iterations ...... =', nni
   write (stdout, '(5x, a, i5)') 'Number of linear iterations ......... =', nli
   write (stdout, '(5x, a, f8.4)') 'Average Krylov subspace dimension =', avdim
   write (stdout, '(x, i5, a, i5, a)') ncfn, ' nonlinear conv. failures,', ncfl, ' linear conv. failures'
   write (stdout, '(5x, a, i7, 1x, i7)') 'Minimum lengths for work arrays WP and IWP: ', lwpmin, liwpmin

end program example_heatilu
