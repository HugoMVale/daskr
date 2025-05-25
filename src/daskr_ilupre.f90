!----------------------------------------------------------------------------------------------
! Adapted from original Fortran code in `original/preconds/dilupre.f`
!----------------------------------------------------------------------------------------------

module daskr_ilupre
!! # Preconditioner Routines for Sparse Problems
!!
!! This module provides a general-purpose sparse incomplete LU (ILU) preconditioner for use with 
!! the [[daskr]] solver, with the Krylov linear system method. 
!!
!! When using [[daskr]] to solve a problem \(G(t,y,\dot{y}) = 0\), whose Jacobian
!! \( J = \partial G/ \partial y + c_J \partial G/ \partial \dot{y} \), where \(c_J\) is a
!! scalar, is a general sparse matrix, the routines [[jac_ilupre]] and [[psol_ilupre]] can be
!! used to generate an approximation to \(J\) as the preconditioner and to solve the resulting
!! sparse linear system, in conjunction with the Krylov method option.
!!
!! Internally, the incomplete LU factorization is achieved via one of two routines — [[ilut]]
!! or [[ilutp]] — from the SPARSKIT library. The routine [[ilut]] performs an ILU factorization
!! of a sparse matrix using a dual thresholding technique based on a drop tolerance (`tolilut`)
!! and a level of fill-in parameter (`lfililut`). The parameter `lfililut` controls the amount
!! of fill-in allowed in the factorization (limited to a maximum of `2*lfililut*neq`, but normally
!! much less). Increasing `lfililut` will generally make the ILU factorization more accurate.
!! The parameter `tolilut` also controls the accuracy of the ILU factorization via a drop tolerance
!! based on element size. Decreasing `tolilut` will increase the amount of fill-in and make for
!! a more accurate factorization. The routine [[ilutp]] is a variant of [[ilut]] that in addition
!! performs pivoting based on a tolerance ratio `permtol`.
!!
!! An important aspect of using incomplete factorization techniques is that of reordering the
!! rows and columns in the Jacobian matrix \(J\) before performing the ILU. In this package,
!! this is accomplished via the parameter `ireorder`, which when equal to 1 performs a reverse
!! Cuthill-McKee (RCM) reordering before performing the ILU factorization. Based on the limited
!! amount of testing done so far, RCM seems the best overall choice. It is possible to include
!! a different reordering technique if desired.
!!
!! ## Usage
!!  
!! To use these routines in conjunction with [[daskr]], the user's calling program should include
!! the following, in addition to setting the other [[daskr]] input parameters:
!
!! * Dimension the array `ipar` to have length at least 30, and load the following parameters
!!   into `ipar` as:
!!   
!!    | Index  | Name         | Description                                                      |
!!    |--------|--------------|------------------------------------------------------------------|
!!    | 1      | `ml`         | The lower bandwidth used in calculating \(J\).                   |
!!    | 2      | `mu`         | The upper bandwidth used in calculating \(J\).                   |
!!    | 3      | `lenpfac`    | The average number of nonzeros in a row of \(J\). The maximum of nonzeros allowed in \(J\) is `nnzmx = lenpfac*neq`. `lenpfac >= 2`. |
!!    | 4      | `lenplufac`  | The average amount of fill-in per row in the factored \(J\). The maximum number of nonzeros allowed in the factored \(J\) is `lenplumx = nnzmx + lenplufac*neq`. `lenplufac >=2`. |
!!    | 5      | `ipremeth`   | Preconditioner type flag. `=1` means [[ilut]] and `=2` means [[ilutp]]. |
!!    | 6      | `lfililut`   | Fill-in parameter for [[ilut]] and [[ilutp]]. The largest `lfililut` elements per row of the L and U factors are kept. Each row of L and U will have a maximum of `lfililut` elements in addition to their original number of nonzero elements. |
!!    | 7      | `ireorder`   | Reordering flag. `=0` means no reordering of \(J\) rows/columns before incomplete factorization. `=1` means reverse Cuthill-McKee (RCM) reordering is performed. |
!!    | 8      | `isrnorm`    | Row norm flag. `=1` means compute and use row norms as scalings in the preconditioner system \(P x = b\). `=0` means no row norm scaling. |
!!    | 9      | `normtype`   | Type of row norm scaling for `isrnorm`. `=0` means max-norm. `=1` means 1-norm. `=2` means 2-norm. |
!!    | 10     | `jacout`     | Output Jacobian flag. `=1` means write \(J\) and initial residual \(G\) to a file (logical unit in `ipar(29)`). Integration halts with `ires = -2`. `=0` means no output. Storage format is Boeing-Harwell. |
!!    | 11     | `jscalcol`   | Flag for scaling \(J\) columns by the inverses of elements in the `ewt` array. `=0` means no scaling. `=1` means perform scaling. |
!!    | 21:28  | —            | Used to hold pointer information.                                |
!!    | 29     | `jacout_unit`| Logical unit number for matrix output file. Used only if `jacout = 1`. |
!!    | 30     | `rescalls`   | On return from [[daskr]], holds the number of calls to the `res` routine used in preconditioner evaluations. | 
!!
!! * Dimension the array `rpar` to have length at least 2, and load the parameters shown in the
!!   table below into `rpar`.
!!  
!!    | Index | Name      | Description                                                          |
!!    |-------|-----------|----------------------------------------------------------------------|
!!    | 1     | `tolilut` | Drop tolerance for use by [[ilut]] and [[ilutp]]. `tolilut >= 0`. Larger values cause less fill-in. Good values range from 0.001 to 0.01.                                   |
!!    | 2     | `permtol` | Tolerance ratio used in determining column pivoting by ILUTP. `permtol >= 0`. Good values are from 0.1 to 0.01. Two columns are permuted only if `a(i,j)*permtol > a(i,i)`. |
!!  
!! * The two parameters `tolilut` and `lfililut` give the user a great deal of flexibility. One 
!!   can use `tolilut = 0` to get a strategy based on keeping the largest elements in each row 
!!   of \(L\) and \(U\). Taking `tolilut /= 0` but `lfililut = neq` will give the usual threshold
!!   strategy (however, fill-in is then unpredictable).
!!
!! * Set `info(12) = 1` to select the Krylov iterative method and `info(15) = 1` to indicate
!!   that a `jac` routine exists. Then in the call to [[daskr]], pass the procedure names
!!   `jac_ilupre` and `psol_ilupre` as the arguments `jac` and `psol`, respectively.
!!
!! * The [[daskr]] work arrays `rwork` and `iwork` must include segments `rwp` and `iwp` for use
!!   by [[jac_ilupre]] and [[psol_ilupre]]. The lengths of these depend on the problem size, 
!!   half-bandwidths, and other parameters as shown in the table below. Load these lengths in
!!   `iwork` as `iwork(27) = lrwp` and `iwork(28) = liwp` and include these values in the declared
!!   size of `rwork` and `iwork`, respectively.
!!
!!    | Variable | Length                                                                              |
!!    |----------|-------------------------------------------------------------------------------------|
!!    | `lrwp`   | `2*lenpfac*neq + lenplufac*neq + isrnorm*neq + neq`                                 |
!!    | `liwp`   | `3*neq + 1 + 3*lenpfac*neq + 2*lenplufac*neq + 2*ireorder*neq + 2*(ipremeth-1)*neq` |
!!      
!! ## Example
!!
!! The program [[example_heatilu]] demonstrates the use of this preconditioner.

   use daskr_kinds, only: rk, zero, one
   use daskr, only: res_t
   implicit none
   private

   public :: setup_ilupre, jac_ilupre, psol_ilupre

contains

   pure subroutine setup_ilupre(neq, lrwp, liwp, rpar, ipar, ierr, lrwp_min, liwp_min)
   !! Setup routine for the incomplete LU preconditioner. This routine checks the user input
   !! and calculates the minimum length needed for the preconditioner workspace arrays.

      integer, intent(in) :: neq
         !! Problem size.
      integer, intent(in) :: lrwp
         !! Current length of `rwp`.
      integer, intent(in) :: liwp
         !! Current length of `iwp`.
      real(rk), intent(in) :: rpar(*)
         !! User real workspace.
      integer, intent(in) :: ipar(*)
         !! User integer workspace.
      integer, intent(out) :: ierr
         !! Error flag (0 means success, else failure):
         !! `1 <= ierr <= 11` means there's an illegal value for `ipar(ierr)`;
         !! `ierr = 12` means `ipar(29)` is illegal;
         !! `ierr = 21` means `rpar(1)` is illegal;
         !! `ierr = 22` means `rpar(2)` is illegal;
         !! `ierr = 30` means more `wp` length is needed;
         !! `ierr = 31` means more `iwp` length is needed.
      integer, intent(out) :: lrwp_min
         !! Minimum `rwp` length needed.
      integer, intent(out) :: liwp_min
         !! Minimum `iwp` length needed.

      integer :: lbw, ubw, lenplumx, ljac, ljaci, ljacj, lrownrms, lrwk1, liwk1, lenpfac, & 
                 lenplufac, lfililut, ipremeth, neqp1, nnzmx, lplu, lju, ljlu, lperm, lqperm, &
                 llevels, lmask
      integer :: isrnorm  ! =1 causes row normalization of JAC.
      integer :: normtype ! =0,1,2 for max-norm, 1-norm, or 2-norm row scaling
      integer :: ireorder ! =1 causes row and column reordering of JAC.
      integer :: jacout   ! =1 causes the Jacobian matrix and SAVR to be written to a file and
                          ! then exit with ierr = 1 to signal a stop to daskr.
      integer :: jscalcol ! =1 causes the columns of the Jacobian matrix to be scaled by EWT-inverse
      real(rk)  :: tolilut, permtol

      ! Load values from IPAR and RPAR. Check for illegal values.
      lbw = ipar(1)        ! LBW must be > 0
      if (lbw <= 0) then
         ierr = 1
         return
      end if

      ubw = ipar(2)        ! UBW must be > 0
      if (ubw <= 0) then
         ierr = 2
         return
      end if

      lenpfac = ipar(3)    ! LENPFAC must be >= 2
      if (lenpfac <= 1) then
         ierr = 3
         return
      end if

      lenplufac = ipar(4)  ! LENPLUFAC must be >= 2
      if (lenplufac <= 1) then
         ierr = 4
         return
      end if

      ipremeth = ipar(5)   ! IPREMETH must be == 1 or 2 currently
      if (ipremeth /= 1 .and. ipremeth /= 2) then
         ierr = 5
         return
      end if

      lfililut = ipar(6)   ! LFILILUT must be >= 0
      if (lfililut < 0) then
         ierr = 6
         return
      end if

      ireorder = ipar(7)   ! IREORDER must be 0 or 1
      if ((ireorder < 0) .or. (ireorder > 1)) then
         ierr = 7
         return
      end if

      isrnorm = ipar(8)    ! ISRNORM must be 0 or 1
      if ((isrnorm < 0) .or. (isrnorm > 1)) then
         ierr = 8
         return
      end if

      normtype = ipar(9)   ! NORMTYPE must be 0, 1, or 2
      if ((normtype < 0) .or. (normtype > 2)) then
         ierr = 9
         return
      end if

      jacout = ipar(10)    ! JACOUT must be 0 or 1
      if ((jacout < 0) .or. (jacout > 1)) then
         ierr = 10
         return
      end if

      jscalcol = ipar(11)    ! JSCALCOL must be 0 or 1
      if ((jscalcol < 0) .or. (jscalcol > 1)) then
         ierr = 11
         return
      end if

      if (jacout == 1) then  ! IPAR(29) must be > 0
         if (ipar(29) <= 0) then
            ierr = 12
            return
         end if
      end if

      tolilut = rpar(1)    ! TOLILUT must be >= 0.0
      if (tolilut < zero) then
         ierr = 21
         return
      end if

      if (ipremeth == 2) then
         permtol = rpar(2)        ! PERMTOL must be >= 0.0
         if (permtol < zero) then
            ierr = 22
            return
         end if
      end if

      ! Calculate minimum work lengths for WP and IWP arrays.
      neqp1 = neq + 1
      nnzmx = lenpfac*neq
      lenplumx = nnzmx + lenplufac*neq

      ! Set up pointers into WP.
      ljac = 1
      lrownrms = nnzmx + ljac
      if (isrnorm == 1) then
         lplu = lrownrms + neq
      else
         lplu = lrownrms
      end if

      lrwk1 = lplu + lenplumx
      lrwp_min = lrwk1 + neq - 1
      if (lrwp < lrwp_min) then
         ierr = 30  ! more WP length needed.
         return
      end if

      ! Set up pointers into IWP
      ljaci = 1
      ljacj = ljaci + neqp1
      lju = ljacj + nnzmx
      ljlu = lju + max(lenplumx, neqp1)
      if (ireorder /= 0) then
         lperm = ljlu + lenplumx
         lqperm = lperm + neq
         liwk1 = lqperm + neq
         llevels = ljlu + nnzmx          ! assumes that LENPLUFAC >= 2.
         lmask = llevels + neq
      else
         lperm = 0
         lqperm = 0
         llevels = 0
         lmask = 0
         liwk1 = ljlu + lenplumx
      end if

      liwp_min = liwk1 + 2*neq - 1
      if (ipremeth == 2) liwp_min = liwp_min + 2*neq
      if (liwp < liwp_min) then
         ierr = 31 ! more IWP length needed.
         return
      end if

      ierr = 0

   end subroutine setup_ilupre

   subroutine jac_ilupre( &
      res, ires, neq, t, y, ydot, rewt, savr, wk, h, cj, rwp, iwp, ierr, rpar, ipar)
   !! This subroutine uses finite-differences to calculate the Jacobian matrix in sparse format,
   !! and then performs an incomplete LU decomposition using either [[ilut]] or [[ilutp]] from
   !! SPARSKIT.

      use dsparskit, only: amudia, dvperm, prtmt, roscal

      procedure(res_t) :: res
         !! User-defined residuals routine.
      integer, intent(out) :: ires
         !! Error flag set by `res`.
      integer, intent(in) :: neq
         !! Problem size.
      real(rk), intent(in) :: t
         !! Current independent variable.
      real(rk), intent(inout) :: y(neq)
         !! Current dependent variables.
      real(rk), intent(inout) :: ydot(neq)
         !! Current derivatives of dependent variables.
      real(rk), intent(in) :: rewt(neq)
         !! Reciprocal error weights for scaling `y` and `ydot`.
      real(rk), intent(inout) :: savr(neq)
         !! Current residual evaluated at `(t, y, ydot)`.
      real(rk), intent(inout) :: wk(neq)
         !! Real work space available to this subroutine.
      real(rk), intent(in) :: h
         !! Current step size.
      real(rk), intent(in) :: cj
         !! Scalar used in forming the system Jacobian.
      real(rk), intent(inout) :: rwp(*)
         !! Matrix elements of ILU.
      integer, intent(inout) :: iwp(*)
         !! Array indices for elements of ILU.
      integer, intent(out) :: ierr
         !! Error flag (0 means success, else failure).
      real(rk), intent(inout) :: rpar(*)
         !! Real array used for communication between the calling program and user routines.
      integer, intent(inout) :: ipar(*)
         !! Integer array used for communication between the calling program and user routines.

      external ::  xerrwd !@todo: replace by module

      character(len=8), parameter :: PMETH(4) = [character(len=8) :: 'ILUT', 'ILUTP', 'ILU0', 'MILU0']
      real(rk) :: tolilut, permtol, sqrtn
      integer :: i, lbw, ubw, lenplumx, ljac, ljaci, ljacj, liperm, lrownrms, lrwk1, liwk1, &
                 ifmt, lenpfac, lenplufac, lfililut, ipremeth, neqp1, nnzmx, lplu, lju, ljlu, &
                 lperm, lqperm, llevels, lmask
      integer :: isrnorm  ! =1 causes row normalization of JAC.
      integer :: normtype ! =0,1,2 for max-norm, 1-norm, or 2-norm row scaling
      integer :: ireorder ! =1 causes row and column reordering of JAC.
      integer :: jacout   ! =1 causes the Jacobian matrix and SAVR to be written to a file and
                          ! then exit with IRES = -2 to signal a stop to DDASPK.
      integer :: iunit    ! logical unit number to use when JACOUT == 1
      integer :: jscalcol ! =1 causes the columns of the Jacobian matrix to be scaled by REWT-inverse
      integer :: nre      ! number of RES calls needed to evaluate Jacobian. NRE is returned in IPAR(30).
      character(len=8) :: premeth
      character(len=72) :: title
      character(len=80) :: msg

      ! Zero out NRE counter (number of RES calls needed to evaluate Jacobian)
      nre = 0

      ! Load values from IPAR and RPAR.
      lbw = ipar(1)
      ubw = ipar(2)
      lenpfac = ipar(3)
      lenplufac = ipar(4)
      ipremeth = ipar(5)
      lfililut = ipar(6)
      ireorder = ipar(7)
      isrnorm = ipar(8)
      normtype = ipar(9)
      jacout = ipar(10)
      jscalcol = ipar(11)
      tolilut = rpar(1)
      permtol = rpar(2)
      premeth = pmeth(ipremeth)

      ! Set pointers into the WP and IWP arrays.
      neqp1 = neq + 1
      nnzmx = lenpfac*neq
      lenplumx = nnzmx + lenplufac*neq

      ! Set up pointers into WP
      ljac = 1
      lrownrms = nnzmx + ljac
      if (isrnorm == 1) then
         lplu = lrownrms + neq
      else
         lplu = lrownrms
      end if
      lrwk1 = lplu + lenplumx

      ! Set up pointers into IWP
      ljaci = 1
      ljacj = ljaci + neqp1
      lju = ljacj + nnzmx
      ljlu = lju + lenplumx

      ! Calculate Jacobian matrix.
      ierr = 0
      call jcalc(neq, t, y, ydot, savr, lbw, ubw, wk, rewt, res, h, cj, nnzmx, rwp(ljac), &
                 iwp(ljacj), iwp(ljaci), rwp(lplu), iwp(ljlu), iwp(lju), ipar, rpar, ires, &
                 nre, ierr)
      if (ires < 0) return
      if (ierr /= 0) return

      ! Save NRE value for user output.
      ipar(30) = ipar(30) + nre

      ! Modify pointers into IWP
      ljlu = lju + neqp1
      if (ireorder /= 0) then
         lperm = ljlu + lenplumx
         lqperm = lperm + neq
         liwk1 = lqperm + neq
         llevels = ljlu + nnzmx          ! assumes that LENPLUFAC >= 2.
         lmask = llevels + neq
      else
         lperm = 0
         lqperm = 0
         llevels = 0
         lmask = 0
         liwk1 = ljlu + lenplumx
      end if

      if (premeth == 'ILUTP') then
         liperm = liwk1 + 2*neq
      else
         liperm = liwk1
      end if

      ! Multiply Jacobian columns by inverse of scaling vector REWT.
      ! In PSOLILU, the WGHT array equals REWT/SQRT(NEQ), so we must be consistent here.
      if (jscalcol == 1) then
         sqrtn = sqrt(real(neq))
         do i = 1, neq
            wk(i) = sqrtn/rewt(i)
         end do
         call amudia(neq, 0, rwp(ljac), iwp(ljacj), iwp(ljaci), wk, rwp(ljac), iwp(ljacj), iwp(ljaci))
      end if

      ! Normalize Jacobian rows, if desired.
      if (isrnorm == 1) then
         call roscal(neq, 0, normtype, rwp(ljac), iwp(ljacj), iwp(ljaci), &
                     rwp(lrownrms), rwp(ljac), iwp(ljacj), iwp(ljaci), ierr)
         if (ierr /= 0) return
      end if

      ! Reorder Jacobian rows and columns, if desired.
      if (ireorder == 1) then
         call jreord(neq, nnzmx, &
                     rwp(ljac), iwp(ljacj), iwp(ljaci), &
                     rwp(lplu), iwp(ljlu), iwp(lju), &
                     iwp(lperm), iwp(lqperm), iwp(llevels), iwp(lmask))
      end if

      ! Write matrix JAC and scaled RES to file if JACOUT == 1.
      if (jacout == 1) then
         iunit = ipar(29)
         if (isrnorm == 1) then
            do i = 1, neq
               savr(i) = savr(i)*rwp(lrownrms + i - 1)
            end do
         end if
         if (ireorder /= 0) call dvperm(neq, savr, iwp(lperm))
         title = 'DASPK Test Matrix '
         ifmt = 15
         call prtmt(neq, neq, rwp(ljac), iwp(ljacj), iwp(ljaci), savr, &
                    'NN', title, 'SPARSKIT', 'RUA', ifmt, 3, iunit)
         msg = 'DJACILU -- Jacobian Matrix written to file.'
         call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
         ierr = 1
         ires = -2
         return
      end if

      ! Compute ILU decomposition.
      call jilu(neq, nnzmx, rwp(ljac), iwp(ljacj), iwp(ljaci), iwp(lju), rwp(lplu), iwp(ljlu), &
                rwp(lrwk1), iwp(liwk1), lenplumx, tolilut, &
                lfililut, permtol, premeth, iwp(liperm), ierr)
      if ((ierr == -2) .or. (ierr == -3)) then
         ires = -2 ! Stop since more storage needed.
      end if

      ! Save pointers for use in DPSOLILU into IPAR array.
      ipar(21) = lplu
      ipar(22) = lju
      ipar(23) = ljlu
      ipar(24) = lrownrms
      ipar(25) = lperm
      ipar(26) = lqperm

   end subroutine jac_ilupre

   subroutine psol_ilupre( &
      neq, t, y, ydot, r0, wk, cj, wght, rwp, iwp, b, epslin, ierr, rpar, ipar)
   !! This subroutine solves the linear system \(P x = b\) for the banded preconditioner \(P\),
   !! given a vector \(b\), using the LU decomposition produced by [[jac_ilupre]]. The solution
   !! is carried out by [[lusol]] from SPARSKIT.

      use dsparskit, only: lusol, dvperm

      integer, intent(in) :: neq
         !! Problem size.
      real(rk), intent(in) :: t
         !! Current independent variable (not used).
      real(rk), intent(in) :: y(*)
         !! Current dependent variables (not used).
      real(rk), intent(in) :: ydot(*)
         !! Current derivatives of dependent variables (not used).
      real(rk), intent(in) :: r0(*)
         !! Current residual evaluated at `(t, y, ydot)` (not used).
      real(rk), intent(in) :: wk(*)
         !! Real work space available to this subroutine.
      real(rk), intent(in) :: cj
         !! Scalar used in forming the system Jacobian.
      real(rk), intent(in) :: wght(*)
         !! Error weights for computing norms.
      real(rk), intent(inout) :: rwp(*)
         !! Matrix elements of ILU.
      integer, intent(inout) :: iwp(*)
         !! Array indices for elements of ILU
      real(rk), intent(inout) :: b(*)
         !! Right-hand side vector on input; solution on output.
      real(rk), intent(in) :: epslin
         !! Tolerance for linear system (not used).
      integer, intent(out) :: ierr
         !! Error flag (not used).
      real(rk), intent(inout) :: rpar(*)
         !! Real array used for communication between the calling program and user routines
         !! (not used).
      integer, intent(inout) :: ipar(*)
         !! Integer array used for communication between the calling program and user routines.

      integer :: i, lplu, lju, ljlu, lrownrms, lperm, lqperm, ireorder, isrnorm, ipremeth, &
                 jscalcol

      ! Load IPREMETH, IREORDER and ISRNORM values from IPAR.
      ipremeth = ipar(5)
      ireorder = ipar(7)
      isrnorm = ipar(8)
      jscalcol = ipar(11)

      ! Load pointers into RWP and IWP arrays.
      lplu = ipar(21)
      lju = ipar(22)
      ljlu = ipar(23)
      lrownrms = ipar(24)
      lperm = ipar(25)
      lqperm = ipar(26)

      ! Scale c by multiplying by row-normalization factors, if used.
      if (isrnorm == 1) then
         do i = 1, neq
            b(i) = b(i)*rwp(lrownrms + i - 1)
         end do
      end if

      ! Solve P*x=b for a preconditioner stored as a sparse matrix in compressed sparse row
      ! format. If rows and columns of P were reordered (permuted), permute b, then use inverse
      ! permutation on x.
      if (ipremeth == 1 .or. ipremeth == 2) then
         if (ireorder == 1) call dvperm(neq, b, iwp(lperm))
         call lusol(neq, b, wk, rwp(lplu), iwp(ljlu), iwp(lju))
         if (ireorder == 1) call dvperm(neq, wk, iwp(lqperm))
      end if

      ! Unscale x by dividing by column scaling vector WGHT.
      if (jscalcol == 1) then
         b(1:neq) = wk(1:neq)/wght(1:neq)
      else
         b(1:neq) = wk(1:neq)
      end if

      ierr = 0

   end subroutine psol_ilupre

   subroutine jcalc( &
      neq, t, y, ydot, r0, ml, mu, r1, rewt, res, h, cj, nnzmx, jac, ja, ia, rcoo, jcoo, &
      icoo, ipar, rpar, ires, nre, ierr)
   !! This subroutine calculates the Jacobian matrix by one-sided finite-differences. Lower and
   !! upper bandwidths are used to select the elements to be computed. The Jacobian is stored
   !! in compressed sparse row format.

      use dsparskit, only: coocsr

      integer, intent(in) :: neq
         !! Problem size.
      real(rk), intent(in) :: t
         !! Current independent variable.
      real(rk), intent(inout) :: y(neq)
         !! Current dependent variables.
      real(rk), intent(inout) :: ydot(neq)
         !! Current derivatives of dependent variables.
      real(rk), intent(in) :: r0(neq)
         !! Current residual evaluated at `(t, y, ydot)`.
      integer, intent(in) :: ml
         !! Lower bandwidth.
      integer, intent(in) :: mu
         !! Upper bandwidth.
      real(rk), intent(inout) :: r1(neq)
         !! Real work space available to this subroutine.
      real(rk), intent(in) :: rewt(neq)
         !! Reciprocal error weights for scaling `y` and `ydot`.
      procedure(res_t) :: res
         !! User-defined residuals routine.
      real(rk), intent(in) ::  h
         !! Current step size.
      real(rk), intent(in) :: cj
         !! Scalar used in forming the system Jacobian.
      integer, intent(in) :: nnzmx
         !! Maximum number of nonzeros in Jacobian.
      real(rk), intent(out) :: jac(nnzmx)
         !! Nonzero Jacobian elements.
      integer, intent(out) :: ja(nnzmx)
         !! Column indices of nonzero Jacobian elements.
      integer, intent(out) :: ia(neq + 1)
         !! Pointers to beginning of each row in `jac` and `ja`.
      real(rk), intent(inout) :: rcoo(nnzmx)
         !! Nonzero Jacobian elements.
      integer, intent(out) :: jcoo(nnzmx)
         !! Column indices of nonzero Jacobian elements.
      integer, intent(out) :: icoo(nnzmx)
         !! Row indices of nonzero Jacobian elements.
      integer, intent(inout) :: ipar(*)
         !! User integer workspace.
      real(rk), intent(inout) :: rpar(*)
         !! User real workspace.
      integer, intent(out) :: ires
         !! Error flag for `res` routine.
      integer, intent(inout) :: nre
         !! Number of calls to `res`.
      integer, intent(out) :: ierr
         !! Error flag.

      external :: xerrwd

      integer :: nnz, i, i1, i2, j, jj, mba, meband, meb1, mband
      real(rk) :: jacelem, squround, del, delinv
      character(len=80) :: msg

      ! Set band parameters.
      nnz = 1
      mband = ml + mu + 1
      mba = min(mband, neq)
      meband = mband + ml
      meb1 = meband - 1

      ! Set the machine unit roundoff UROUND and SQRT(UROUND), used to set increments in the
      ! difference quotient procedure.
      squround = sqrt(epsilon(one))

      ! Initial error flags.
      ierr = 0
      ires = 0

      ! Make MBA calls to RES to approximate the Jacobian.
      ! Here, R0(1),...,R0(neq) contains the base RES value, and
      ! R1(1),...,R1(NEQ) contains the perturbed values of RES.
      do j = 1, mba

         do jj = j, neq, mband
            jac(jj) = y(jj)
            jac(jj + neq) = ydot(jj)
            del = squround*max(abs(y(jj)), abs(h*ydot(jj)), abs(one/rewt(jj)))
            del = sign(del, h*ydot(jj))
            del = (y(jj) + del) - y(jj)
            y(jj) = y(jj) + del
            ydot(jj) = ydot(jj) + cj*del
         end do

         call res(t, y, ydot, cj, r1, ires, rpar, ipar)
         if (ires < 0) return

         nre = nre + 1
         do jj = j, neq, mband
            y(jj) = jac(jj)
            ydot(jj) = jac(jj + neq)
            del = squround*max(abs(y(jj)), abs(h*ydot(jj)), abs(one/rewt(jj)))
            del = sign(del, h*ydot(jj))
            del = (y(jj) + del) - y(jj)
            delinv = one/del
            i1 = max(1, (jj - mu))
            i2 = min(neq, (jj + ml))
            do i = i1, i2
               ! Calculate possibly nonzero Jacobian elements for this variable,
               ! and store nonzero elements in coordinate format.
               jacelem = (r1(i) - r0(i))*delinv
               if (jacelem /= zero) then
                  if (nnz > nnzmx) then
                     msg = 'DJCALC -- More storage needed for Jacobian.'
                     call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
                     msg = 'DJCALC -- Increase LENPFAC.'
                     call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
                     msg = 'DJCALC -- Storage exceeded at (I,J) = (I1,I2)'
                     call xerrwd(msg, 80, 0, 0, 2, i, jj, 0, 0.0, 0.0)
                     ierr = 1
                     ires = -2
                     return
                  end if
                  rcoo(nnz) = jacelem
                  jcoo(nnz) = jj
                  icoo(nnz) = i
                  nnz = nnz + 1
               end if
            end do
         end do
      end do
      nnz = nnz - 1

      ! Convert Jacobian from coordinate to compressed sparse row format.
      call coocsr(neq, nnz, rcoo, icoo, jcoo, jac, ja, ia)

   end subroutine jcalc

   subroutine jilu( &
      neq, nnzmx, jac, ja, ia, ju, plu, jlu, rwk1, iwk1, lenplumx, tolilut, &
      lfililut, permtol, premeth, iperm, ierr)
   !! This subroutine computes the incomplete LU decomposition of the Jacobian matrix and returns
   !! it in any given storage format.

      use dsparskit, only: ilut, ilutp
      
      integer, intent(in) :: neq
         !! Problem size.
      integer, intent(in) :: nnzmx
         !! Maximum number of nonzeros in Jacobian.
      real(rk), intent(in) :: jac(nnzmx)
         !! Nonzero Jacobian elements.
      integer, intent(in) :: ja(nnzmx)
         !! Column indices of nonzero Jacobian elements.
      integer, intent(in) :: ia(neq + 1)
         !! Pointers to beginning of each row in `jac` and `ja`.
      integer, intent(out) :: ju(neq)
         !! Pointer to beginning of each row of U in matrix `plu` and `jlu`.
      real(rk), intent(out) :: plu(lenplumx)
         !! Matrix elements of ILU.
      integer, intent(out) :: jlu(lenplumx)
         !! Sizes and array indices for elements of ILU.
      real(rk), intent(inout) :: rwk1(neq)
         !! Real work space.
      integer, intent(inout) :: iwk1(2*neq)
         !! Integer work space.
      integer, intent(in) :: lenplumx
         !! Length of `plu` and `jlu` arrays.
      real(rk), intent(in) :: tolilut
         !! Tolerance for ILU decomposition.
      integer, intent(in) :: lfililut
         !! @todo: find description
      real(rk), intent(in) :: permtol
         !! Tolerance for permutation.
      character(len=8), intent(in) :: premeth
         !! Preconditioner method.
      integer, intent(inout) :: iperm(2*neq)
         !! Integer work space.
      integer, intent(out) :: ierr
         !! Error flag (0 means success, else failure).

      external :: xerrwd ! @todo: replace by module

      character(len=80) :: msg
      logical :: err

      err = .false.

      if (premeth == 'ILUT') then
         ! Use incomplete factorization routine ILUT from SparsKit.
         call ilut(neq, jac, ja, ia, lfililut, tolilut, plu, jlu, ju, lenplumx, rwk1, iwk1, ierr)
         if (ierr /= 0) then
            msg = 'DJILU -- Error return from ILUT: IERR = (I1)'
            call xerrwd(msg, 80, 0, 0, 1, ierr, 0, 0, 0.0, 0.0)
            err = .true.
         end if
      elseif (premeth == 'ILUTP') then
         ! Use incomplete factorization routine ILUTP from SparsKit.
         call ilutp(neq, jac, ja, ia, lfililut, tolilut, permtol, neq, plu, jlu, ju, lenplumx, &
                    rwk1, iwk1, iperm, ierr)
         if (ierr /= 0) then
            msg = 'DJILU -- Error return from ILUTP: IERR = (I1)'
            call xerrwd(msg, 80, 0, 0, 1, ierr, 0, 0, 0.0, 0.0)
            err = .true.
         end if
      end if

      ! Put in other options here for incomplete factorizations.
      ! @todo: refactor this
      if (err) then
         msg = 'DJILU -- IERR /= 0 means one of the following has occurred:'
         call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
         msg = '    IERR >  0   --> Zero pivot encountered at step number IERR.'
         call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
         msg = '    IERR = -1   --> Error. input matrix may be wrong.'
         call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
         msg = '                     (The elimination process has generated a'
         call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
         msg = '                     row in L or U with length > NEQ.)'
         call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
         msg = '    IERR = -2   --> Matrix L overflows.'
         call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
         msg = '    IERR = -3   --> Matrix U overflows.'
         call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
         msg = '    IERR = -4   --> Illegal value for LFILILUT.'
         call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
         msg = '    IERR = -5   --> Zero row encountered.'
         call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
         msg = '    '
         call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
         msg = '    For IERR = -2 or -3, increase the value of LENPLUFAC or'
         call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
         msg = '    decrease the value of LFILILUT if LENPLUFAC cannot be'
         call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
         msg = '    increased.'
         call xerrwd(msg, 80, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
      end if

   end subroutine jilu

   subroutine jreord(neq, nnzmx, jac, ja, ia, awk, jwk, iwk, perm, qperm, levels, mask)
   !! This subroutine reorders the Jacobian matrix.

      use dsparskit, only: atob, bfs, rversp, dperm

      integer, intent(in) :: neq
         !! Problem size.
      integer, intent(in) :: nnzmx
         !! Maximum number of nonzeroes in Jacobian.
      real(rk), intent(in) :: jac(nnzmx)
         !! Nonzero Jacobian elements.
      integer, intent(in) :: ja(nnzmx)
         !! Column indices of nonzero Jacobian elements.
      integer, intent(in) :: ia(neq + 1)
         !! Indices of 1st nonzero element in each row.
      real(rk), intent(inout) :: awk(nnzmx)
         !! Work array.
      integer, intent(inout) :: jwk(nnzmx)
         !! Work array.
      integer, intent(inout) :: iwk(neq + 1)
         !! Work array.
      integer, intent(inout) :: perm(neq)
         !! Array containing the permutation used in reordering the rows and columns of the
         !! Jacobian matrix.
      integer, intent(inout) :: qperm(neq)
         !! Integer array holding the inverse of the permutation in array `perm`.
      integer, intent(inout) :: levels(neq)
         !! Work array used by the [[bfs]] reordering subroutine.
      integer, intent(inout) :: mask(neq)
         !! Work array used by the [[bfs]] reordering subroutine.

      integer :: i, nfirst, nlev, maskval

      ! Copy JAC, JA, and IA to AWK, JWK, and IWK.
      call atob(neq, jac, ja, ia, awk, jwk, iwk)

      ! Perform a Cuthill-McKee reordering of the Jacobian.
      nfirst = 1
      perm(1) = 0
      mask = 1
      maskval = 1
      qperm(1) = 1
      call bfs(neq, jwk, iwk, nfirst, perm, mask, maskval, qperm, levels, nlev)

      ! Reverse the permutation to obtain the reverse Cuthill-McKee reordering.
      call rversp(neq, qperm)

      ! Calculate the inverse of QPERM and put it in PERM.
      do i = 1, neq
         perm(qperm(i)) = i
      end do

      ! Permute rows and columns of Jacobian using PERM.
      call dperm(neq, awk, jwk, iwk, jac, ja, ia, perm, perm, 1)

   end subroutine jreord

end module daskr_ilupre
