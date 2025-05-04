!----------------------------------------------------------------------------------------------
! Adapted from original Fortran code in `original/preconds/dbanpre.f`
!----------------------------------------------------------------------------------------------

module daskr_banpre
!! Preconditioner Routines for Banded Problems
!!
!! The following pair of subroutines – [[banja]] and [[banps]] – provides a general-purpose 
!! banded preconditioner matrix for use with the [[daskr]] solver, with the Krylov linear system
!! method. When using [[daskr]] to solve a problem \(G(t,y,y') = 0\), whose iteration matrix 
!! (Jacobian) \( J = dG/dy + c_j  dG/dy' \), where \(c\) is a scalar, is either banded or 
!! approximately equal to a banded matrix, these routines can be used to generate a banded 
!! approximation to \(J\) as the preconditioner and to solve the resulting banded linear system,
!! in conjunction with the Krylov method option (`info(12) = 1`) in [[daskr]].
!!
!! Other than the user-supplied residual routine `res` defining \(G(t,y,y')\), the only other 
!! inputs required by these routines are the half-bandwidth parameters \(\mathrm{ml}\) and 
!! \(\mathrm{mu}\) of the approximate banded Jacobian. If the system size is \(\mathrm{neq}\),
!! the half-bandwidths are defined as integers between 0 and \(\mathrm{neq} -1\) such that only
!! elements with indices \((i,j)\) satisfying \( -\mathrm{ml} \le j, i \le \mathrm{mu} \) are
!! to be retained in the preconditioner. For example, if \(\mathrm{ml} = \mathrm{mu} = 0\), a
!! diagonal matrix will be generated as the preconditioner.  The banded preconditioner is 
!! obtained by difference quotient approximations. If the true problem Jacobian is not banded 
!! but is approximately equal to a matrix that is banded, the procedure used here will have the
!! effect of lumping the elements outside of the band onto the elements within the band.
!!
!! To use these routines in conjunction with [[daskr]], the user's calling program should include
!! the following, in addition to setting the other [[daskr]] input parameters:
!!
!! * Dimension the array `ipar` to have length at least 2, and load the half-bandwidths into 
!!   `ipar` as `ipar(1) = ml` and `ipar(2) = mu`. Array `ipar` is used to communicate these 
!!   parameters to [[banja]] and [[banps]]. If the user program also uses `ipar` for communication 
!!   with `res`, that data should be located beyond the first 2 positions.
!!
!! * Import this module. Set `info(15) = 1` to indicate that a `jac` routine exists. Then in the
!!   call to [[daskr]], pass the procedure names `banja` and `banps` as the arguments `jac` and
!!   `psol`, respectively.
!!
!! * The [[daskr]] work arrays `rwork` and `iwork` must include segments `wp` and `iwp` for use
!!   by  [[banja]] and [[banps]]. The lengths of these arrays depend on the problem size and 
!!   half-bandwidths, as follows:
!!```  
!!   lwp  = length of rwork segment wp  = (2*ml + mu + 1)*neq + 2*((neq/(ml + mu + 1)) + 1)
!!   liwp = length of iwork segment iwp = neq
!!```
!!   Note the integer divide in `lwp`. Load these lengths in `iwork` as `iwork(27) = lwp` and
!!   `iwork(28) = liwp`, and include these values in the declared size of `rwork` and `iwork`.
!!
!! The [[banja]] and [[banps]] routines generate and solve the banded preconditioner matrix 
!! \(P\) within the preconditioned Krylov algorithm used by [[daskr]] when `info(12) = 1`. \(P\)
!! is generated and LU-factored periodically during the integration, and the factors are used to
!! solve systems \(Px = b\) as needed.

   use daskr_kinds, only: rk, one
   implicit none
   private

   public :: banja, banps

contains

   subroutine banja(res, ires, neq, t, y, yprime, rewt, savres, wk, h, cj, wp, iwp, ier, rpar, ipar)
   !! This subroutine generates a banded preconditioner matrix \(P\) that approximates the
   !! iteration matrix \(J = dG/dy + c_j dG/dy'\), where the DAE system is \(G(t,y,y') = 0\). 
   !! The band matrix \(P\) has half-bandwidths \(\mathrm{ml}\) and \(\mathrm{mu}\). It is computed
   !! by making \(\mathrm{ml} + \mathrm{mu} + 1\) calls to the user's `res` routine and forming
   !! difference quotients, exactly as in the banded direct method option of [[daskr]]. 
   !! [[banja]] calls the LINPACK routine [[DGBFA]] to do an LU factorization of this matrix.
      external :: res
      integer, intent(out) :: ires
        !! Output flag set by `res`. See `res` description in [[daskr]].
      integer, intent(in) :: neq
        !! Problem size.
      real(rk), intent(in) :: t
        !! Independent variable.
      real(rk), intent(inout) :: y(*)
        !! Current dependent variables.
      real(rk), intent(inout) :: yprime(*)
        !! Current derivatives of dependent variables.
      real(rk), intent(in) :: rewt(*)
        !! Vector of reciprocal error weights, used here for computing increments.
      real(rk), intent(in) :: savres(*)
        !! Current residual evaluated at `(t, y, yprime)`.
      real(rk), intent(in) :: wk(*)
        !! Real work space of length `neq`.
      real(rk), intent(in) :: h
        !! Current step size.
      real(rk), intent(in) :: cj
        !! Scalar proportional to `1/h`.
      real(rk), intent(inout) :: wp(*)
        !! Real work array for P, etc. On output, it contains the LU decomposition of the banded
        !! approximation P.
      integer, intent(inout) :: iwp(*)
        !! Integer work space for matrix pivot information.
      integer, intent(out) :: ier
        !! Output flag: `ier > 0` if P is singular, and `ier = 0` otherwise.
      real(rk), intent(inout) :: rpar(*)
        !! Real array used for communication between the calling program and external user
        !! routines.
      integer, intent(inout) :: ipar(*)
        !! Integer array used for communication between the calling program and external user
        !! routines. `ipar(1)` and `ipar(2)` must contain `ml` and `mu`, respectively.

      external :: dgbfa

      real(rk) :: del, delinv, squround
      integer :: i, i1, i2, ii, ipsave, isave, j, k, lenp, mba, mband, meb1, meband, ml, &
                 msave, mu, n

      ! Set band parameters.
      ml = ipar(1)
      mu = ipar(2)
      mband = ml + mu + 1
      mba = min(mband, neq)
      meband = mband + ml
      meb1 = meband - 1

      ! Set the machine unit roundoff UROUND and SQRT(UROUND), used to set increments in the
      ! difference quotient procedure.
      squround = sqrt(epsilon(one))

      ! Set pointers into WP. LENP is the length of the segment for P.
      ! Following that are two segments of size (NEQ/MBAND), with offsets
      ! ISAVE and IPSAVE, for temporary storage of Y and YPRIME elements.
      lenp = (2*ml + mu + 1)*neq
      msave = (neq/mband) + 1
      isave = lenp
      ipsave = isave + msave

      ! Initialize error flags.
      ier = 0
      ires = 0

      ! Generate the banded approximate iteration matrix P using difference quotients on the
      ! results of calls to RES.
      do j = 1, mba

         do n = j, neq, mband
            k = (n - j)/mband + 1
            wp(isave + k) = y(n)
            wp(ipsave + k) = yprime(n)
            del = squround*max(abs(y(n)), abs(h*yprime(n)), abs(one/rewt(n)))
            del = sign(del, h*yprime(n))
            del = (y(n) + del) - y(n)
            y(n) = y(n) + del
            yprime(n) = yprime(n) + cj*del
         end do
         
         call res(t, y, yprime, cj, wk, ires, rpar, ipar)
         
         if (ires < 0) return
         
         do n = j, neq, mband
            k = (n - j)/mband + 1
            y(n) = wp(isave + k)
            yprime(n) = wp(ipsave + k)
            del = squround*max(abs(y(n)), abs(h*yprime(n)), abs(one/rewt(n)))
            del = sign(del, h*yprime(n))
            del = (y(n) + del) - y(n)
            delinv = one/del
            i1 = max(1, n - mu)
            i2 = min(neq, n + ml)
            ii = n*meb1 - ml
            do i = i1, i2
               wp(ii + i) = (wk(i) - savres(i))*delinv
            end do
         end do
      
      end do

      ! Do LU decomposition of the band matrix P.
      call dgbfa(wp, meband, neq, ml, mu, iwp, ier)

   end subroutine banja

   subroutine banps(neq, t, y, yprime, savres, wk, cj, wght, wp, iwp, b, eplin, ier, rpar, ipar)
   !! This subroutine uses the factors produced by [[banja]] to solve linear systems \(P x = b\)
   !! for the banded preconditioner \(P\), given a vector \(b\). It calls the LINPACK routine
   !! [[DGBSL]] for this.
      integer, intent(in) :: neq
        !! Problem size.
      real(rk), intent(in) :: t
        !! Independent variable (not used).
      real(rk), intent(in) :: y(*)
        !! Current dependent variables (not used).
      real(rk), intent(in) :: yprime(*)
        !! Current derivatives of dependent variables (not used).
      real(rk), intent(in) :: savres(*)
        !! Current residual evaluated at `(t, y, yprime)` (not used).
      real(rk), intent(in) :: wk(*)
        !! Real work space of length `neq` (not used).
      real(rk), intent(in) :: cj
        !! Scalar proportional to `1/h` (not used).
      real(rk), intent(in) :: wght(*)
        !! Error weights for computing norms (not used).
      real(rk), intent(inout) :: wp(*)
        !! Real work array containing the LU decomposition of P.
      integer, intent(inout) :: iwp(*)
        !! Integer array containing matrix pivot information.
      real(rk), intent(inout) :: b(*)
        !! Right-hand side vector on input; solution on output.
      real(rk), intent(in) :: eplin
        !! Tolerance for linear system (not used).
      integer, intent(out) :: ier
        !! Output error flag (not used; assumed 0 on input).
      real(rk), intent(inout) :: rpar(*)
        !! Real array used for communication between the calling program and external user
        !! routines (not used).
      integer, intent(inout) :: ipar(*)
        !! Integer array used for communication between the calling program and external user
        !! routines (not used). `ipar(1)` and `ipar(2)` must contain `ml` and `mu`, respectively.

      external :: dgbsl
      
      integer :: meband, ml, mu

      ml = ipar(1)
      mu = ipar(2)
      meband = 2*ml + mu + 1
      call dgbsl(wp, meband, neq, ml, mu, iwp, b, 0)

   end subroutine banps

end module daskr_banpre
