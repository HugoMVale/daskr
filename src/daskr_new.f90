!----------------------------------------------------------------------------------------------
! Adapted from original Fortran code in `original/solver/ddaskr.f`
!----------------------------------------------------------------------------------------------

subroutine dslvk( &
   neq, y, t, ydot, savr, x, ewt, rwm, iwm, res, ires, psol, &
   iersl, cj, epslin, sqrtn, rsqrtn, rhok, rpar, ipar)
!! This routine uses a restart algorithm and interfaces to [[dspigm]] for
!! the solution of the linear system arising from a Newton iteration.

   use daskr_kinds, only: rk, zero
   use blas_interfaces, only: dscal, dcopy
   implicit none

   integer, intent(in) :: neq
      !! Problem size.
   real(rk), intent(in) :: y(neq)
      !! Current dependent variables.
   real(rk), intent(in) :: t
      !! Current time.
   real(rk), intent(in) :: ydot(neq)
      !! Current derivatives of dependent variables.
   real(rk), intent(in) :: savr(neq)
      !! Current residual evaluated at `(t, y, ydot)`.
   real(rk), intent(inout) :: x(neq)
      !! On entry, the right-hand side vector of the linear system to be solved,
      !! and on exit, the solution vector.
   real(rk), intent(inout) :: ewt(neq) ! @todo: harmonize names
      !! Nonzero elements of the diagonal scaling matrix.
   real(rk), intent(inout) :: rwm(*)
      !! Real work space containing data for the algorithm (Krylov basis vectors,
      !! Hessenberg matrix, etc.).
   integer, intent(inout) :: iwm(*)
      !! Integer work space containing data for the algorithm.
   external :: res
      !! Residuals routine.
   integer, intent(out) :: ires
      !! Error flag from `res`.
   external :: psol
      !! Preconditioner routine.
   integer, intent(out) :: iersl
      !! Error flag.
      !! `iersl = 0` means no trouble occurred (or user `res` routine returned `ires < 0`).
      !! `iersl = 1` means the iterative method failed to converge ([[dspigm]] returned `iflag > 0`).
      !! `iersl = -1` means there was a nonrecoverable error in the iterative solver, and an error exit will occur.
   real(rk), intent(in) :: cj
      !! Scalar used in forming the system Jacobian.
   real(rk), intent(in) :: epslin
      !! Tolerance for linear system solver.
   real(rk), intent(in) :: sqrtn
      !! Square root of `neq`.
   real(rk), intent(in) :: rsqrtn
      !! Reciprocal of square root of `neq`.
   real(rk), intent(out) :: rhok
      !! Weighted norm of the final preconditioned residual.
   real(rk), intent(inout) :: rpar(*)
      !! User real workspace.
   integer, intent(inout) :: ipar(*)
      !! User integer workspace.

   integer, parameter :: lnre = 12, lncfl = 16, lnli = 20, lnps = 21, &
                         llocwp = 29, llciwp = 30, &
                         lmiter = 23, lmaxl = 24, lkmp = 25, lnrmax = 26

   integer, save :: irst = 1
   integer :: i, iflag, kmp, ldl, lhes, lgmr, liwp, lq, lr, lv, lwk, lwp, lz, &
              maxl, maxlp1, miter, ncfl, nli, nps, npsol, nre, nres, nrmax, nrsts

   liwp = iwm(llciwp)
   nli = iwm(lnli)
   nps = iwm(lnps)
   ncfl = iwm(lncfl)
   nre = iwm(lnre)
   lwp = iwm(llocwp)
   maxl = iwm(lmaxl)
   kmp = iwm(lkmp)
   nrmax = iwm(lnrmax)
   miter = iwm(lmiter)
   iersl = 0
   ires = 0

   ! Use a restarting strategy to solve the linear system P*X = -F. Parse the work vector,
   ! and perform initializations. Note that zero is the initial guess for X.
   maxlp1 = maxl + 1
   lv = 1
   lr = lv + neq*maxl
   lhes = lr + neq + 1
   lq = lhes + maxl*maxlp1
   lwk = lq + 2*maxl
   ldl = lwk + min(1, maxl - kmp)*neq
   lz = ldl + neq
   call dscal(neq, rsqrtn, ewt, 1)
   call dcopy(neq, x, 1, rwm(lr), 1)
   x = zero

   ! Top of loop for the restart algorithm. Initial pass approximates X and sets up a
   ! transformed system to perform subsequent restarts to update X. NRSTS is initialized
   ! to -1, because restarting does not occur until after the first pass.
   ! Update NRSTS; conditionally copy DL to R; call the DSPIGM algorithm to solve A*Z = R;
   ! updated counters; update X with the residual solution.
   ! Note: if convergence is not achieved after NRMAX restarts, then the linear solver is
   ! considered to have failed.
   iflag = 1
   nrsts = -1
   do while ((iflag .eq. 1) .and. (nrsts .lt. nrmax) .and. (ires .eq. 0))
      nrsts = nrsts + 1
      if (nrsts .gt. 0) call dcopy(neq, rwm(ldl), 1, rwm(lr), 1)
      call dspigm(neq, t, y, ydot, savr, rwm(lr), ewt, maxl, &
                  kmp, epslin, cj, res, ires, nres, psol, npsol, rwm(lz), rwm(lv), &
                  rwm(lhes), rwm(lq), lgmr, rwm(lwp), iwm(liwp), rwm(lwk), &
                  rwm(ldl), rhok, iflag, irst, nrsts, rpar, ipar)
      nli = nli + lgmr
      nps = nps + npsol
      nre = nre + nres
      do i = 1, neq
         x(i) = x(i) + rwm(lz + i - 1)
      end do
   end do

   ! The restart scheme is finished. Test IRES and IFLAG to see if convergence was not
   ! achieved, and set flags accordingly.
   if (ires .lt. 0) then
      ncfl = ncfl + 1
   elseif (iflag .ne. 0) then
      ncfl = ncfl + 1
      if (iflag .gt. 0) iersl = 1
      if (iflag .lt. 0) iersl = -1
   end if

   ! Update IWM with counters, rescale EWT, and return.
   iwm(lnli) = nli
   iwm(lnps) = nps
   iwm(lncfl) = ncfl
   iwm(lnre) = nre
   call dscal(neq, sqrtn, ewt, 1)

end subroutine dslvk

subroutine dspigm( &
   neq, t, y, ydot, savr, r, wght, maxl, &
   kmp, epslin, cj, res, ires, nres, psol, npsol, z, v, &
   hes, q, lgmr, rwp, iwp, wk, dl, rhok, iflag, irst, nrsts, &
   rpar, ipar)
!! This routine solves the linear system \(A z = r\) using a scaled preconditioned version
!! of the generalized minimum residual method. An initial guess of \(z=0\) is assumed.

   use daskr_kinds, only: rk, zero, one
   use blas_interfaces, only: dnrm2, dcopy, dscal, daxpy
   implicit none

   integer, intent(in) :: neq
      !! Problem size.
   real(rk), intent(in) :: t
      !! Current time.
   real(rk), intent(in) :: y(neq)
      !! Current dependent variables.
   real(rk), intent(in) :: ydot(neq)
      !! Current derivatives of dependent variables.
   real(rk), intent(in) :: savr(neq)
      !! Current residual evaluated at `(t, y, ydot)`.
   real(rk), intent(inout) :: r(neq)
      !! On entry, the right hand side vector. Also used as work space when computing
      !! the final approximation and will therefore be destroyed. `r` is the same as
      !! `v(:,maxl+1)` in the call to this routine.
   real(rk), intent(in) :: wght(neq)
      !! Nonzero elements of the diagonal scaling matrix.
   integer, intent(in) :: maxl
      !! Maximum allowable order of the matrix `hes`.
   integer, intent(in) :: kmp
      !! Number of previous vectors the new vector `vnew` must be made orthogonal to
      !! (`kmp <= maxl`).
   real(rk), intent(in) :: epslin
      !! Tolerance on residuals \(Az - r\) in weighted rms norm.
   real(rk), intent(in) :: cj
      !! Scalar used in forming the system Jacobian.
   external :: res
      !! Residuals routine.
   integer, intent(out) :: ires
      !! Error flag from `res`.
   integer, intent(inout) :: nres
      !! Number of calls to `res`.
   external :: psol
      !! Preconditioner routine.
   integer, intent(inout) :: npsol
      !! Number of calls to `psol`.
   real(rk), intent(out) :: z(neq)
      !! Final computed approximation to the solution of the system \(Az = r\).
   real(rk), intent(out) :: v(neq, *)
      !! Matrix of shape `(neq,lgmr+1)` containing the `lgmr` orthogonal vectors `v(:,1)`
      !! to `v(:,lgmr)`.
   real(rk), intent(out) :: hes(maxl+1, maxl)
      !! Upper triangular factor of the QR decomposition of the upper Hessenberg matrix
      !! whose entries are the scaled inner-products of `A*v(:,i)` and `v(:,k)`.
   real(rk), intent(out) :: q(2*maxl)
      !! Components of the Givens rotations used in the QR decomposition of `hes`. It is
      !! loaded in [[dheqr]] and used in [[dhels]].
   integer, intent(out) :: lgmr
      !! The number of iterations performed and the current order of the upper
      !! Hessenberg matrix `hes`.
   real(rk), intent(inout) :: rwp(*)
      !! Real work array used by preconditioner `psol`.
   integer, intent(inout) :: iwp(*)
      !! Integer work array used by preconditioner `psol`.
   real(rk), intent(inout) :: wk(*)
      !! Real work array used by [[datv]] and `psol`.
   real(rk), intent(inout) :: dl(*)
      !! Real work array used for calculation of the residual norm `rhok` when the
      !! method is incomplete (`kmp < maxl`) and/or when using restarting.
   real(rk), intent(out) :: rhok
      !! Weighted norm of the final preconditioned residual.
   integer, intent(out) :: iflag
      !! Error flag.
      !! 0 means convergence in `lgmr` iterations, `lgmr <= maxl`.
      !! 1 means the convergence test did not pass in `maxl` iterations, but the new
      !! residual norm (`rho`) is less than the old residual norm (`rnrm`), and so `z`
      !! is computed.
      !! 2 means the convergence test did not pass in `maxl` iterations, new residual
      !! norm (`rho`) is greater than or equal to old residual norm (`rnrm`), and the
      !! initial guess, `z = 0`, is returned.
      !! 3 means there was a recoverable error in `psol` caused by the preconditioner
      !! being out of date.
      !! -1 means there was an unrecoverable error in `psol`.
   integer, intent(in) :: irst
      !! Restarting flag. If `irst > 0`, then restarting is being performed.
   integer, intent(in) :: nrsts
      !! Number of restarts on the current call to [[dspigm]]. If `nrsts > 0`, then the
      !! residual `r` is already scaled, and so scaling of `r` is not necessary.
   real(rk), intent(inout) :: rpar(*)
      !! User real workspace.
   integer, intent(inout) :: ipar(*)
      !! User integer workspace.

   integer :: i, ierr, info, ip1, i2, k, ll, llp1, maxlm1, maxlp1
   real(rk) :: c, dlnrm, prod, rho, rnrm, s, snormw, tem

   ierr = 0
   iflag = 0
   lgmr = 0
   npsol = 0
   nres = 0

   ! The initial guess for Z is 0. The initial residual is therefore
   ! the vector R. Initialize Z to 0.
   z = zero

   ! Apply inverse of left preconditioner to vector R if NRSTS .EQ. 0.
   ! Form V(*,1), the scaled preconditioned right hand side.
   if (nrsts .eq. 0) then
      call psol(neq, t, y, ydot, savr, wk, cj, wght, rwp, iwp, r, epslin, ierr, rpar, ipar)
      npsol = 1
      if (ierr .ne. 0) goto 300
      v(:, 1) = r*wght
   else
      v(:, 1) = r
   end if

   ! Calculate norm of scaled vector V(*,1) and normalize it
   ! If, however, the norm of V(*,1) (i.e. the norm of the preconditioned
   ! residual) is .le. EPLIN, then return with Z=0.
   rnrm = dnrm2(neq, v, 1)
   if (rnrm .le. epslin) then
      rhok = rnrm
      return
   end if
   tem = one/rnrm
   call dscal(neq, tem, v(1, 1), 1)

   ! Zero out the HES array.
   hes = zero

   ! Main loop to compute the vectors V(*,2) to V(*,MAXL).
   ! The running product PROD is needed for the convergence test.
   prod = one
   maxlp1 = maxl + 1
   do ll = 1, maxl
      lgmr = ll

      ! Call routine DATV to compute VNEW = ABAR*V(LL), where ABAR is
      ! the matrix A with scaling and inverse preconditioner factors applied.
      call datv(neq, y, t, ydot, savr, v(1, ll), wght, z, &
                res, ires, psol, v(1, ll + 1), wk, rwp, iwp, cj, epslin, &
                ierr, nres, npsol, rpar, ipar)
      if (ires .lt. 0) return
      if (ierr .ne. 0) goto 300

      ! Call routine DORTH to orthogonalize the new vector VNEW = V(*,LL+1).
      call dorth(v(1, ll + 1), v, hes, neq, ll, maxlp1, kmp, snormw)
      hes(ll + 1, ll) = snormw

      ! Call routine DHEQR to update the factors of HES.
      call dheqr(hes, maxlp1, ll, q, info, ll)
      if (info .eq. ll) goto 120

      ! Update RHO, the estimate of the norm of the residual R - A*ZL.
      ! If KMP .LT. MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
      ! necessarily orthogonal for LL .GT. KMP.  The vector DL must then
      ! be computed, and its norm used in the calculation of RHO.
      prod = prod*q(2*ll)
      rho = abs(prod*rnrm)
      if ((ll .gt. kmp) .and. (kmp .lt. maxl)) then
         if (ll .eq. kmp + 1) then
            call dcopy(neq, v(1, 1), 1, dl, 1)
            do i = 1, kmp
               ip1 = i + 1
               i2 = i*2
               s = q(i2)
               c = q(i2 - 1)
               do k = 1, neq
                  dl(k) = s*dl(k) + c*v(k, ip1)
               end do
            end do
         end if
         s = q(2*ll)
         c = q(2*ll - 1)/snormw
         llp1 = ll + 1
         do k = 1, neq
            dl(k) = s*dl(k) + c*v(k, llp1)
         end do
         dlnrm = dnrm2(neq, dl, 1)
         rho = rho*dlnrm
      end if

      ! Test for convergence. If passed, compute approximation ZL.
      ! If failed and LL .LT. MAXL, then continue iterating.
      if (rho .le. epslin) goto 200
      if (ll .eq. maxl) goto 100

      ! Rescale so that the norm of V(1,LL+1) is one.
      tem = one/snormw
      call dscal(neq, tem, v(1, ll + 1), 1)

   end do

100 continue
   if (rho .lt. rnrm) goto 150

120 continue
   iflag = 2
   z = zero
   return

150 continue
   ! The tolerance was not met, but the residual norm was reduced.
   ! If performing restarting (IRST .gt. 0) calculate the residual vector
   ! RL and store it in the DL array.  If the incomplete version is
   ! being used (KMP .lt. MAXL) then DL has already been calculated.
   iflag = 1
   if (irst .gt. 0) then

      if (kmp .eq. maxl) then
         !  Calculate DL from the V(I)'s.
         call dcopy(neq, v(1, 1), 1, dl, 1)
         maxlm1 = maxl - 1
         do i = 1, maxlm1
            ip1 = i + 1
            i2 = i*2
            s = q(i2)
            c = q(i2 - 1)
            do k = 1, neq
               dl(k) = s*dl(k) + c*v(k, ip1)
            end do
         end do
         s = q(2*maxl)
         c = q(2*maxl - 1)/snormw
         do k = 1, neq
            dl(k) = s*dl(k) + c*v(k, maxlp1)
         end do
      end if

      ! Scale DL by RNRM*PROD to obtain the residual RL.
      tem = rnrm*prod
      call dscal(neq, tem, dl, 1)
   end if

   ! Compute the approximation ZL to the solution.
   ! Since the vector Z was used as work space, and the initial guess
   ! of the Newton correction is zero, Z must be reset to zero.
200 continue
   ll = lgmr
   llp1 = ll + 1
   r(1:llp1) = zero
   r(1) = rnrm
   call dhels(hes, maxlp1, ll, q, r)
   z = zero
   do i = 1, ll
      call daxpy(neq, r(i), v(1, i), 1, z, 1)
   end do
   z = z/wght
   ! Load RHO into RHOK.
   rhok = rho
   return

   ! This block handles error returns forced by routine PSOL.
300 continue
   if (ierr .lt. 0) iflag = -1
   if (ierr .gt. 0) iflag = 3

end subroutine dspigm

subroutine datv( &
   neq, y, t, ydot, savr, v, wght, ydottemp, res, &
   ires, psol, z, vtemp, rwp, iwp, cj, epslin, ierr, nres, npsol, rpar, ipar)
!! This routine computes the matrix-vector product:
!!
!! $$ z = D^{-1} P^{-1} (\partial F / \partial y) D v $$
!!
!! where \(F = G(t, y, c_J(y - a))\), \(c_J\) is a scalar proportional to \(1/h\), and \(a\)
!! involves the past history of \(y\). The quantity \(c_J(y - a)\) is an approximation to 
!! the first derivative of \(y\) and is stored in \(\dot{y}\). 
!!
!! \(D\) is a diagonal scaling matrix, and \(P\) is the left preconditioning matrix. \(v\)
!! is assumed to have L-2 norm equal to 1. The product is stored in \(z\) and is computed by
!! means of a difference quotient, a call to `res`, and one call to `psol`.

   use daskr_kinds, only: rk

   integer, intent(in) :: neq
      !! Problem size.
   real(rk), intent(in) :: y(neq)
      !! Current dependent variables.
   real(rk), intent(in) :: t
      !! Current time.
   real(rk), intent(in) :: ydot(neq)
      !! Current derivatives of dependent variables.
   real(rk), intent(in) :: savr(neq)
      !! Current residual evaluated at `(t, y, ydot)`.
   real(rk), intent(in) :: v(neq)
      !! Orthonormal vector (can be the same as `z`).  
   real(rk), intent(in) :: wght(neq)
      !! Scaling factors. `1/wght(i)` are the diagonal elements of the matrix `D`.
   real(rk), intent(out) :: ydottemp(neq)
      !! Work array used to store the incremented value of `ydot`.
   external :: res
      !! Residuals routine. 
   integer, intent(out) :: ires
      !! Error flag from `res`.
   external :: psol
      !! Preconditioner routine.
   real(rk), intent(out) :: z(neq)
      !! Desired scaled matrix-vector product.
   real(rk), intent(out) :: vtemp(neq)
      !! Work array used to store the unscaled version of `v`.
   real(rk), intent(inout) :: rwp(*)
      !! Real work array used by preconditioner `psol`.
   integer, intent(inout) :: iwp(*)
      !! Integer work array used by preconditioner `psol`.
   real(rk), intent(in) :: cj
      !! Scalar used in forming the system Jacobian.
   real(rk), intent(in) :: epslin
      !! Tolerance for linear system.
   integer, intent(out) :: ierr
      !! Error flag from `psol`.
   integer, intent(inout) :: nres
      !! Number of calls to `res`.
   integer, intent(inout) :: npsol
      !! Number of calls to `psol`.
   real(rk), intent(inout) :: rpar(*)
      !! User real workspace.
   integer, intent(inout) :: ipar(*)
      !! User integer workspace.

   ires = 0
   ierr = 0

   ! Set VTEM = D * V.
   vtemp = v/wght

   ! Store Y in Z and increment Z by VTEMP.
   ! Store YDOT in YDOTTEMP and increment YDOTTEMP by VTEM*CJ.
   ydottemp = ydot + vtemp*cj
   z = y + vtemp

   ! Call RES with incremented Y, YDOT arguments
   ! stored in Z, YDOTTEMP. VTEMP is overwritten with new residual.
   call res(t, z, ydottemp, cj, vtemp, ires, rpar, ipar)
   nres = nres + 1
   if (ires .lt. 0) return

   ! Set Z = (dF/dY) * VBAR using difference quotient.
   ! (VBAR is old value of VTEMP before calling RES)
   z = vtemp - savr

   ! Apply inverse of left preconditioner to Z.
   call psol(neq, t, y, ydot, savr, ydottemp, cj, wght, rwp, iwp, z, epslin, ierr, rpar, ipar)
   npsol = npsol + 1
   if (ierr .ne. 0) return

   ! Apply D-inverse to Z
   z = z * wght

end subroutine datv

pure subroutine dorth(vnew, v, hes, n, ll, ldhes, kmp, snormw)
!! This routine orthogonalizes the vector `vnew` against the previous `kmp` vectors in the
!! `v` matrix. It uses a modified Gram-Schmidt orthogonalization procedure with conditional
!! reorthogonalization.

   use daskr_kinds, only: rk, zero
   use blas_interfaces, only: daxpy, ddot, dnrm2
   implicit none

   real(rk), intent(inout) :: vnew(n)
      !! On entry, vector containing a scaled product of the Jacobian and the vector `v(*,ll)`.
      !! On return, the new vector orthogonal to `v(*,i0)`, where `i0 = max(1, ll - kmp + 1)`.
   real(rk), intent(in) :: v(n, ll)
      !! Matrix containing the previous `ll` orthogonal vectors `v(*,1)` to `v(*,ll)`.
   real(rk), intent(inout) :: hes(ldhes, ll)
      !! On entry, an upper Hessenberg matrix of shape `(ll, ll)` containing in `hes(i,k)`,
      !! for `k < ll`, the scaled inner products of `a*v(*,k)` and `v(*,i)`.
      !! On return, an upper Hessenberg matrix with column `ll` filled in with the scaled
      !! inner products of `a*v(*,ll)` and `v(*,i)`.
   integer, intent(in) :: n
      !! Order of the matrix `a`.
   integer, intent(in) :: ll
      !! Current order of the matrix `hes`.
   integer, intent(in) :: ldhes
      !! Leading dimension of `hes`.
   integer, intent(in) :: kmp
      !! Number of previous vectors the new vector `vnew` must be made orthogonal to
      !! (`kmp <= maxl`).
   real(rk), intent(out) :: snormw
      !! L-2 norm of `vnew`.

   integer :: i, i0
   real(rk) :: arg, sumdsq, tem, vnrm

   ! Get norm of unaltered VNEW for later use.
   vnrm = dnrm2(n, vnew, 1)

   ! Do Modified Gram-Schmidt on VNEW = A*V(LL).
   ! Scaled inner products give new column of HES.
   ! Projections of earlier vectors are subtracted from VNEW.
   i0 = max(1, ll - kmp + 1)
   do i = i0, ll
      hes(i, ll) = ddot(n, v(1, i), 1, vnew, 1)
      tem = -hes(i, ll)
      call daxpy(n, tem, v(1, i), 1, vnew, 1)
   end do

   ! Compute SNORMW = norm of VNEW.
   ! If VNEW is small compared to its input value (in norm), then
   ! Reorthogonalize VNEW to V(*,1) through V(*,LL).
   ! Correct if relative correction exceeds 1000*(unit roundoff).
   ! Finally, correct SNORMW using the dot products involved.
   snormw = dnrm2(n, vnew, 1)
   if (vnrm + snormw/1000 .ne. vnrm) return ! @todo: fix this comparison

   sumdsq = zero
   do i = i0, ll
      tem = -ddot(n, v(1, i), 1, vnew, 1)
      if (hes(i, ll) + tem/1000 .eq. hes(i, ll)) cycle ! @todo: fix this comparison
      hes(i, ll) = hes(i, ll) - tem
      call daxpy(n, tem, v(1, i), 1, vnew, 1)
      sumdsq = sumdsq + tem**2
   end do
   if (sumdsq .eq. zero) return

   arg = max(zero, snormw**2 - sumdsq)
   snormw = sqrt(arg)

end subroutine dorth

pure subroutine dheqr(a, lda, n, q, info, ijob)
!! This routine performs a QR decomposition of an upper Hessenberg matrix `a` using Givens
!! rotations. There are two options available:
!!
!! 1. Performing a fresh decomposition;
!! 2. Updating the QR factors by adding a row and a column to the matrix `a`.

   use daskr_kinds, only: rk, zero, one
   implicit none

   real(rk), intent(inout) :: a(lda, n)
      !! On entry, the Hessenberg matrix to be decomposed. On return, the upper triangular
      !! matrix R from the QR decomposition of `a`.
   integer, intent(in) :: lda
      !! Leading dimension of `a`.
   integer, intent(in) :: n
      !! Number of columns in `a` (originally a matrix with shape `(n+1, n)`).
   real(rk), intent(out) :: q(2*n)
      !! Coefficients of the Givens rotations used in decomposing `a`.
   integer, intent(out) :: info
      !! Info flag.
      !! `info = 0`, if successful.
      !! `info = k`, if `a(k,k) = 0`; this is not an error condition for this subroutine, but
      !! it does indicate that [[dhels]] will divide by zero if called.
   integer, intent(in) :: ijob
      !! Job flag.
      !! `ijob = 1`, means that a fresh decomposition of the matrix `a` is desired.
      !! `ijob >= 2`, means that the current decomposition of `a` will be updated by the
      !! addition of a row and a column.

   integer :: i, iq, j, k, km1, kp1, nm1
   real(rk) :: c, s, t, t1, t2

   ! A new factorization is desired.
   if (ijob == 1) then

      ! QR decomposition without pivoting.
      info = 0
      do k = 1, n
         km1 = k - 1
         kp1 = k + 1

         ! Compute Kth column of R.
         ! First, multiply the Kth column of A by the previous K-1 Givens rotations.
         if (km1 .lt. 1) goto 20
         do j = 1, km1
            i = 2*(j - 1) + 1
            t1 = a(j, k)
            t2 = a(j + 1, k)
            c = q(i)
            s = q(i + 1)
            a(j, k) = c*t1 - s*t2
            a(j + 1, k) = s*t1 + c*t2
         end do

         ! Compute Givens components C and S.
20       continue
         iq = 2*km1 + 1
         t1 = a(k, k)
         t2 = a(kp1, k)
         if (t2 .ne. zero) goto 30
         c = one
         s = zero
         goto 50

30       continue
         if (abs(t2) .lt. abs(t1)) goto 40
         t = t1/t2
         s = -one/sqrt(one + t*t)
         c = -s*t
         goto 50

40       continue
         t = t2/t1
         c = one/sqrt(one + t*t)
         s = -c*t

50       continue
         q(iq) = c
         q(iq + 1) = s
         a(k, k) = c*t1 - s*t2
         if (a(k, k) .eq. zero) info = k
      end do

      ! The old factorization of A will be updated.  A row and a column
      ! has been added to the matrix A.
      ! N by N-1 is now the old size of the matrix.
   else if (ijob >= 2) then

      nm1 = n - 1

      ! Multiply the new column by the N previous Givens rotations.
      do k = 1, nm1
         i = 2*(k - 1) + 1
         t1 = a(k, n)
         t2 = a(k + 1, n)
         c = q(i)
         s = q(i + 1)
         a(k, n) = c*t1 - s*t2
         a(k + 1, n) = s*t1 + c*t2
      end do

      ! Complete update of decomposition by forming last Givens rotation,
      ! and multiplying it times the column vector (A(N,N),A(NP1,N)).
      info = 0
      t1 = a(n, n)
      t2 = a(n + 1, n)
      if (t2 .ne. zero) goto 110
      c = one
      s = zero
      goto 130

110   continue
      if (abs(t2) .lt. abs(t1)) goto 120
      t = t1/t2
      s = -one/sqrt(one + t*t)
      c = -s*t
      goto 130

120   continue
      t = t2/t1
      c = one/sqrt(one + t*t)
      s = -c*t

130   continue
      iq = 2*n - 1
      q(iq) = c
      q(iq + 1) = s
      a(n, n) = c*t1 - s*t2
      if (a(n, n) .eq. zero) info = n

   end if

end subroutine dheqr

pure subroutine dhels(a, lda, n, q, b)
!! This routine solves the least squares problem:
!! 
!!  $$ \min{\| b - A x \|^2} $$
!! 
!! using the factors computed by [[dheqr]]. This is similar to the LINPACK routine [[dgesl]]
!! except that \(A\) is an upper Hessenberg matrix.

   use daskr_kinds, only: rk
   use blas_interfaces, only: daxpy
   implicit none

   real(rk), intent(in) :: a(lda, n)
      !! Output from [[dheqr]] which contains the upper triangular factor R in the QR
      !! decomposition of `a`.
   integer, intent(in) :: lda
      !! Leading dimension of `a`.
   integer, intent(in) :: n
      !! Number of columns in `a` (originally a matrix with shape `(n+1, n)`).
   real(rk), intent(in) :: q(2*n)
      !! Coefficients of the Givens rotations used in decomposing `a`.
   real(rk), intent(inout) :: b(n+1)
      !! On entry, the right hand side vector. On return, the solution vector \(x\).

   integer :: iq, k, kb, kp1
   real(rk) :: c, s, t, t1, t2

   ! Minimize (B-A*X,B-A*X).
   ! First form Q*B.
   do k = 1, n
      kp1 = k + 1
      iq = 2*(k - 1) + 1
      c = q(iq)
      s = q(iq + 1)
      t1 = b(k)
      t2 = b(kp1)
      b(k) = c*t1 - s*t2
      b(kp1) = s*t1 + c*t2
   end do

   ! Now solve R*X = Q*B.
   do kb = 1, n
      k = n + 1 - kb
      b(k) = b(k)/a(k, k)
      t = -b(k)
      call daxpy(k - 1, t, a(1, k), 1, b(1), 1)
   end do

end subroutine dhels
