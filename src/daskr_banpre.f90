!----------------------------------------------------------------------------------------------
! Adapted from original Fortran code in `original/preconds/dbanpre.f`
!----------------------------------------------------------------------------------------------

module daskr_banpre
!! # Preconditioner Routines for Banded Problems
!!
!! This module provides a general-purpose banded preconditioner for use with the [[daskr]] solver,
!! with the Krylov linear system method.
!! 
!! When using [[daskr]] to solve a problem \(G(t,y,\dot{y}) = 0\), whose Jacobian
!! \( J = \partial G/ \partial y + c_J \partial G/ \partial \dot{y} \), where \(c_J\) is a
!! scalar, is either banded or approximately equal to a banded matrix, the routines [[jac_banpre]]
!! and [[psol_banpre]] can be used to generate a banded approximation to \(J\) as the preconditioner
!! and to solve the resulting banded linear system, in conjunction with the Krylov method option.
!!  
!! Other than the user-supplied residual routine `res` defining \(G(t,y,\dot{y})\), the only
!! other inputs required by these routines are the half-bandwidth parameters \(\mathrm{ml}\) and 
!! \(\mathrm{mu}\) of the approximate banded Jacobian. If the system size is \(\mathrm{neq}\),
!! the half-bandwidths are defined as integers between 0 and \(\mathrm{neq} -1\) such that only
!! elements with indices \((i,j)\) satisfying \( -\mathrm{ml} \le j, i \le \mathrm{mu} \) are
!! to be retained in the preconditioner. For example, if \(\mathrm{ml} = \mathrm{mu} = 0\), a
!! diagonal matrix will be generated as the preconditioner. The banded preconditioner is 
!! obtained by difference quotient approximations. If the true problem Jacobian is not banded 
!! but is approximately equal to a matrix that is banded, the procedure used here will have the
!! effect of lumping the elements outside of the band onto the elements within the band.
!!
!! ## Usage
!!  
!! To use these routines in conjunction with [[daskr]], the user's calling program should include
!! the following, in addition to setting the other [[daskr]] input parameters:
!!
!! * Dimension the array `ipar` to have length at least 2, and load the half-bandwidths into 
!!   `ipar` as `ipar(1) = ml` and `ipar(2) = mu`. Array `ipar` is used to communicate these 
!!   parameters to [[jac_banpre]] and [[psol_banpre]]. If the user program also uses `ipar` for 
!!   communication with `res`, that data should be located beyond the first 2 positions.
!!
!! * Set `info(12) = 1` to select the Krylov iterative method and `info(15) = 1` to indicate
!!   that a `jac` routine exists. Then in the call to [[daskr]], pass the procedure names
!!   `jac_banpre` and `psol_banpre` as the arguments `jac` and `psol`, respectively.
!!
!! * The [[daskr]] work arrays `rwork` and `iwork` must include segments `rwp` and `iwp` for
!!   use by  [[jac_banpre]] and [[psol_banpre]]. The lengths of these arrays depend on the 
!!   problem size and half-bandwidths, as as shown in the table below. Note the integer divide
!!   in `lrwp`. Load these lengths in `iwork` as `iwork(27) = lrwp` and `iwork(28) = liwp`, and
!!   include these values in the declared size of `rwork` and `iwork`, respectively.
!!
!!    | Variable | Length                                              |
!!    |----------|-----------------------------------------------------|     
!!    | `lrwp`   | `(2*ml + mu + 1)*neq + 2*((neq/(ml + mu + 1)) + 1)` |
!!    | `liwp`   | `neq`                                               |
!!
!! ## Example
!!
!! The program [[example_heat]] demonstrates the use of this preconditioner.  

   use daskr_kinds, only: rk, one
   use daskr, only: res_t
   implicit none
   private

   public :: jac_banpre, psol_banpre

contains

   subroutine jac_banpre( &
      res, ires, neq, t, y, ydot, rewt, savr, wk, h, cj, rwp, iwp, ierr, rpar, ipar)
   !! This routine generates a banded preconditioner matrix \(P\) that approximates the Jacobian
   !! matrix \(J\). The banded matrix \(P\) has half-bandwidths \(\mathrm{ml}\) and \(\mathrm{mu}\),
   !! and is computed by making \(\mathrm{ml} + \mathrm{mu} + 1\) calls to the user's `res`
   !! routine and forming difference quotients, exactly as in the banded direct method option
   !! of [[daskr]]. Afterwards, this matrix is LU factorized by [[dgbfa]] from LINPACK and the
   !! factors are stored in the work arrays `rwp` and `iwp`.

      use dlinpack, only: dgbfa

      procedure(res_t) :: res
        !! User-defined residuals routine.
      integer, intent(out) :: ires
        !! Error flag set by `res`.
      integer, intent(in) :: neq
        !! Problem size.
      real(rk), intent(in) :: t
        !! Current independent variable.
      real(rk), intent(inout) :: y(*)
        !! Current dependent variables.
      real(rk), intent(inout) :: ydot(*)
        !! Current derivatives of dependent variables.
      real(rk), intent(in) :: rewt(*)
        !! Reciprocal error weights for scaling `y` and `ydot`.
      real(rk), intent(inout) :: savr(*)
        !! Current residual evaluated at `(t, y, ydot)`.
      real(rk), intent(inout) :: wk(*)
        !! Real work space available to this subroutine.
      real(rk), intent(in) :: h
        !! Current step size.
      real(rk), intent(in) :: cj
        !! Scalar used in forming the system Jacobian.
      real(rk), intent(inout) :: rwp(*)
        !! Real work array for \(P\), etc. On output, it contains the LU decomposition of the 
        !! banded approximation \(P\).
      integer, intent(inout) :: iwp(*)
        !! Integer work space for matrix pivot information.
      integer, intent(inout) :: ierr
        !! Error flag: `ierr > 0` if \(P\) is singular, and `ierr = 0` otherwise.
      real(rk), intent(inout) :: rpar(*)
        !! Real array used for communication between the calling program and external user
        !! routines.
      integer, intent(inout) :: ipar(*)
        !! Integer array used for communication between the calling program and external user
        !! routines. `ipar(1)` and `ipar(2)` must contain `ml` and `mu`, respectively.

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
      ! ISAVE and IPSAVE, for temporary storage of Y and YDOT elements.
      lenp = (2*ml + mu + 1)*neq
      msave = (neq/mband) + 1
      isave = lenp
      ipsave = isave + msave

      ! Initialize error flags.
      ierr = 0
      ires = 0

      ! Generate the banded approximate iteration matrix P using difference quotients on the
      ! results of calls to RES.
      do j = 1, mba

         do n = j, neq, mband
            k = (n - j)/mband + 1
            rwp(isave + k) = y(n)
            rwp(ipsave + k) = ydot(n)
            del = squround*max(abs(y(n)), abs(h*ydot(n)), abs(one/rewt(n)))
            del = sign(del, h*ydot(n))
            del = (y(n) + del) - y(n)
            y(n) = y(n) + del
            ydot(n) = ydot(n) + cj*del
         end do
         
         call res(t, y, ydot, cj, wk, ires, rpar, ipar)
         
         if (ires < 0) return
         
         do n = j, neq, mband
            k = (n - j)/mband + 1
            y(n) = rwp(isave + k)
            ydot(n) = rwp(ipsave + k)
            del = squround*max(abs(y(n)), abs(h*ydot(n)), abs(one/rewt(n)))
            del = sign(del, h*ydot(n))
            del = (y(n) + del) - y(n)
            delinv = one/del
            i1 = max(1, n - mu)
            i2 = min(neq, n + ml)
            ii = n*meb1 - ml
            do i = i1, i2
               rwp(ii + i) = (wk(i) - savr(i))*delinv
            end do
         end do
      
      end do

      ! Do LU decomposition of the band matrix P.
      call dgbfa(rwp, meband, neq, ml, mu, iwp, ierr)

   end subroutine jac_banpre

   subroutine psol_banpre( &
      neq, t, y, ydot, savr, wk, cj, wght, rwp, iwp, b, epslin, ierr, rpar, ipar)
   !! This routine solves the linear system \(P x = b\) for the banded preconditioner \(P\),
   !! given a vector \(b\), using the LU decomposition produced by [[jac_banpre]]. The solution
   !! is carried out by [[dgbsl]] from LINPACK.

      use dlinpack, only: dgbsl

      integer, intent(in) :: neq
        !! Problem size.
      real(rk), intent(in) :: t
        !! Current independent variable (not used).
      real(rk), intent(in) :: y(*)
        !! Current dependent variables (not used).
      real(rk), intent(in) :: ydot(*)
        !! Current derivatives of dependent variables (not used).
      real(rk), intent(in) :: savr(*)
        !! Current residual evaluated at `(t, y, ydot)` (not used).
      real(rk), intent(in) :: wk(*)
        !! Real work space available to this subroutine (not used).
      real(rk), intent(in) :: cj
        !! Scalar used in forming the system Jacobian (not used).
      real(rk), intent(in) :: wght(*)
        !! Error weights for computing norms (not used).
      real(rk), intent(inout) :: rwp(*)
        !! Real work array containing the LU decomposition of \(P\).
      integer, intent(inout) :: iwp(*)
        !! Integer array containing matrix pivot information.
      real(rk), intent(inout) :: b(*)
        !! Right-hand side vector on input; solution on output.
      real(rk), intent(in) :: epslin
        !! Tolerance for linear system (not used).
      integer, intent(inout) :: ierr
        !! Error flag.
      real(rk), intent(inout) :: rpar(*)
        !! Real array used for communication between the calling program and user routines
        !! (not used).
      integer, intent(inout) :: ipar(*)
        !! Integer array used for communication between the calling program and user routines 
        !! (not used).
      
      integer :: meband, ml, mu

      ml = ipar(1)
      mu = ipar(2)
      meband = 2*ml + mu + 1
      
      ierr = 0
      call dgbsl(rwp, meband, neq, ml, mu, iwp, b, 0)

   end subroutine psol_banpre

end module daskr_banpre
