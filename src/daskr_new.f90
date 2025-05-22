!----------------------------------------------------------------------------------------------
! Adapted from original Fortran code in `original/solver/ddaskr.f`
!----------------------------------------------------------------------------------------------

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
!! is assumed to have L-2 norm equal to 1. The product is stored in `z` and is computed by
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
!! 1. performing a fresh decomposition;
!! 2. updating the QR factors by adding a row and a column to the matrix `a`.

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
!! This routine solves the least squares problem
!!
!!       $$ min (b - a x, b - a x) $$
!!
!! using the factors computed by [[dheqr]]. This is similar to the LINPACK routine [[DGESL]]
!! except that `a` is an upper Hessenberg matrix.

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
      !! On entry, the right hand side vector. On return, the solution vector x.

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
