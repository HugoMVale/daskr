!----------------------------------------------------------------------------------------------
! Adapted from original Fortran code in `original/preconds/drbgpre.f`
!----------------------------------------------------------------------------------------------

module daskr_rbgpre
!! # Preconditioner Tools for Reaction-Transport Problems. Part II: Block-Grouping in Block-Diagonal Reaction-Based Factor
!!
!! This module is provided to assist in the generation and solution of preconditioner matrices
!! for problems arising from reaction-transport systems. More specifically, the routines
!! [[jac_rbgpre]], and [[psol_rbgpre]] are intended as *auxiliary* routines for the user-supplied
!! routines `jac` and `psol` called by [[daskr]] when the Krylov method is selected.
!! 
!! These routines are meant for a DAE system obtained from a system of reaction-transport PDEs,
!! in which some of the PDE variables obey differential equations, and the rest obey algebraic
!! equations [2, section 4]. It is assumed that the right-hand sides of all the equations have
!! the form of a sum of a reaction term \(R\) and a transport term \(S\), that the transport
!! term is discretized by finite differences, and that in the spatial discretization the PDE
!! variables at each spatial point are kept together. Thus, the DAE system function, in terms
!! of a dependent variable vector \(u\), has the form:
!!   
!!   $$  G(t,u,\dot{u}) = I_d \dot{u} - R(t,u) - S(t,u) $$
!!
!! where \(I_d\) is the identity matrix with zeros in the positions corresponding to the
!! algebraic components and ones in those for the differential components, \(R(t,u)\) denotes
!! the reaction terms (spatial coupling absent), and \(S(t,u)\) denotes the spatial transport
!! terms.     
!! 
!! As shown in Brown et al. [2], two possible preconditioners for such a system are:
!!   
!! * \(P_R = c_J I_d - \partial R/ \partial u\), based on the reaction term \(R\) alone.
!! * \(P_{SR} = (I - c_J^{-1} \partial S/\partial u) (c_J I_d - \partial R / \partial u)\), the
!!   product of two factors (in either order), one being \(P_R\) and the other being based on
!!   the spatial term \(S\) alone.
!!
!! The routines given here can be used to compute the reaction-based factor \(P_R\). More
!! precisely, they provide an approximation \(A_R\) to \(P_R\). The matrix \(P_R\) is block-diagonal,
!! with each block corresponding to one spatial point. In \(A_R\), we compute each block by
!! difference quotient approximations, by way of calls to a user-supplied subroutine `rblock`,
!! that evaluates the reaction terms at a single spatial point. In addition, rather than evaluating
!! a block at every spatial point in the mesh, we use a block-grouping scheme, described in
!! Brown et al. [1]. In this scheme, the mesh points are grouped, as in domain decomposition,
!! and only one block of \(\partial R / \partial u)\) is computed for each group; then in solving
!! \(A_R x = b\), the inverse of the representative block is applied to all the blocks of unknowns
!! in the group. Block-grouping greatly reduces the storage required for the preconditioner.
!!
!! The routines given here are specialized to the case of a 2D problem on a rectangular mesh in
!! the x-y plane, and for a block-grouping arrangement that is rectangular (i.e. the Cartesian
!! product of two 1D groupings). However, they can be easily modified for a different problem
!! geometry or a different grouping arrangement.  It is also assumed that the PDE variables are
!! ordered so that the differential variables appear first, followed by the algebraic variables.
!!
!! ## Usage
!!
!! To use these routines in conjunction with [[daskr]], the user's calling program should include
!! the following, in addition to setting the other [[daskr]] input parameters:
!!
!! * Set `info(12) = 1` to select the Krylov iterative method and `info(15) = 1` to indicate
!!   that a `jac` routine exists.
!!         
!! * Call [[setup_rbgpre]] to set mesh-related and block-grouping inputs.
!!   
!! * A `jac` routine, as prescribed by the [[daskr]] instructions, which calls [[jac_rbgpre]],
!!   and does any other Jacobian-related preprocessing needed for preconditioning. The latter
!!   routine takes as argument `r0` the current value of the reaction vector \(R\). This can be
!!   done, for example, by taking `r0` to be `rpar`, and loading `rpar` with the vector \(R\)
!!   in the last call to the `res` routine; in that case, the calling program must declare `rpar`
!!   to have length at least `neq`. Alternatively, insert a call to `rblock` within the loop 
!!   over mesh points inside [[jac_rbgpre]].
!!   
!! * A `psol` routine, as prescribed by the [[daskr]] instructions, which calls [[psol_rbgpre]]
!!   for the solution of systems \(A_R x = b\), and does any other linear system solving required
!!   by the preconditioner.
!!
!! ## Example
!!
!! The program [[example_web]] demonstrates the use of this preconditioner.
!!   
!! ## References
!!   
!! 1. Peter N. Brown and Alan C. Hindmarsh,
!!    Reduced Storage Matrix Methods in Stiff ODE Systems,
!!    J. Appl. Math. & Comp., 31 (1989), pp. 40-91.
!! 2. Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
!!    Using Krylov Methods in the Solution of Large-Scale Differential-Algebraic Systems,
!!    SIAM J. Sci. Comput., 15 (1994), pp. 1467-1488.
!! 3. Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
!!    Consistent Initial Condition Calculation for Differential-Algebraic Systems,
!!    SIAM Journal on Scientific Computing 19.5 (1998): 1495-1512.

   use daskr_kinds, only: rk, one, half
   implicit none
   private

   ! Module variables >>Thread Unsafe<<
   integer :: mp, mpd, mpsq, meshx, meshy, mxmp, ngx, ngy, ngrp
   integer, allocatable :: jgx(:), jgy(:), jigx(:), jigy(:), jxr(:), jyr(:)
   real(rk) :: srur

   public :: setup_rbgpre, jac_rbgpre, psol_rbgpre

contains

   subroutine setup_rbgpre(mx, my, ns, nsd, nxg, nyg, lid, iwork)
   !! This routine sets the mesh and block-grouping parameters needed to use the routines
   !! [[jac_rbgpre]] and [[psol_rbgpre]], assuming a 2D rectangular problem with uniform
   !! rectangular block-grouping. Additionally, it loads the lengths `lrwp` and `liwp` in array
   !! `iwork`.
      integer, intent(in) :: mx
         !! Number of mesh points in x-direction.
      integer, intent(in) :: my
         !! Number of mesh points in y-direction.
      integer, intent(in) :: ns
         !! Number of PDE variables, the size of each block in the block-diagonal preconditioner
         !! matrix \(P_R\).
      integer, intent(in) :: nsd
         !! Number of differential PDE variables. In the DAE system, the first `nsd` variables
         !! at each spatial point have time derivatives, and the remaining `(ns - nsd)` do not.
      integer, intent(in) :: nxg
         !! Number of groups in x-direction.
      integer, intent(in) :: nyg
         !! Number of groups in y-direction.
      integer, intent(in) :: lid
         !! Flag indicating whether to load the ID array in `iwork`. If `lid > 0`, set the ID
         !! array in `iwork`, indicating which components are differential and which are algebraic.
         !! This value is required if either `info(11) = 1` or `info(16) = 1`, in which case
         !! set `lid = 40` or `40 + neq`, depending on the value of the constraint option
         !! `info(10)`. Otherwise, set `lid = 0`.
      integer, intent(inout) :: iwork(*)
         !! Integer work array.

      integer :: i, i0, jx, jy, maxm

      ! Load the modules variables.
      srur = sqrt(epsilon(one))
      mp = ns
      mpd = nsd
      mpsq = ns*ns
      meshx = mx
      meshy = my
      mxmp = meshx*mp
      ngx = nxg
      ngy = nyg
      ngrp = ngx*ngy

      maxm = max(meshx, meshy)
      allocate(jgx(maxm + 1), jgy(maxm + 1), jigx(maxm), jigy(maxm), jxr(maxm), jyr(maxm))
   
      ! Call set_grouping for each mesh direction to load grouping arrays.
      call set_grouping(meshx, ngx, jgx, jigx, jxr)
      call set_grouping(meshy, ngy, jgy, jigy, jyr)

      ! Set the sizes of the preconditioning storage space segments in RWORK and IWORK.
      iwork(27) = mpsq*ngrp
      iwork(28) = mp*ngrp

      ! If LID > 0, set the ID array in IWORK.
      if (lid == 0) return
      i0 = lid
      do jy = 1, my
         do jx = 1, mx
            do i = 1, mpd
               iwork(i0 + i) = 1
            end do
            do i = mpd + 1, mp
               iwork(i0 + i) = -1
               i0 = i0 + mp
            end do
         end do
      end do

   end subroutine setup_rbgpre

   pure subroutine set_grouping(m, ng, jg, jig, jr)
   !! This routine sets arrays `jg`, `jig`, and `jr` describing a uniform (or nearly uniform)
   !! partition of (1, 2,...,`m`) into `ng` groups.
      integer, intent(in) :: m
         !! Number of mesh points in the direction being grouped.
      integer, intent(in) :: ng
         !! Number of groups in the direction being grouped.
      integer, intent(inout) :: jg(:)
         !! Array of length `ng + 1` of group boundaries in the direction being grouped.
         !! Group `ig` has indices `j = jg(ig),...,jg(ig+1)-1`.
      integer, intent(inout) :: jig(:)
         !! Array of length `m` of group indices vs mesh point index.
         !! Mesh point index `j` is in group `jig(j)`.
      integer, intent(inout) :: jr(:)
         !! Array of length `ng` of representative indices for the groups.
         !! The index for group `ig` is `jr(ig)`.

      integer :: ig, j, len1, mper, ngm1

      mper = m/ng
      do ig = 1, ng
         jg(ig) = 1 + (ig - 1)*mper
      end do
      jg(ng + 1) = m + 1

      ngm1 = ng - 1
      len1 = ngm1*mper
      do j = 1, len1
         jig(j) = 1 + (j - 1)/mper
      end do

      len1 = len1 + 1
      do j = len1, m
         jig(j) = ng
      end do

      ! I couldn't tell for sure if the original intention was to roundoff or truncate...
      ! I assumed the latter, while fixing the mixed precision issues.
      do ig = 1, ngm1
         !jr(ig) = 0.5d0 + (real(ig) - 0.5d0)*real(mper)
         jr(ig) = int(half + (real(ig, rk) - half)*mper)
      end do
      !jr(ng) = 0.5d0*real(1 + ngm1*mper + m)
      jr(ng) = int(half*(1 + ngm1*mper + m))

   end subroutine set_grouping
 
   subroutine jac_rbgpre(t, u, r0, rblock, r1, rewt, cj, bd, ipbd, ierr)
   !! This routine generates and preprocesses a block-diagonal preconditioner matrix, based on
   !! the part of the Jacobian corresponding to the reaction terms \(R\) of the problem, using
   !! block-grouping. It generates a matrix of the form \(c_J I_d - \partial R/ \partial u\).
   !! It calls [[dgefa]] from LINPACK to do the LU decomposition of each diagonal block. The
   !! computation of the diagonal blocks uses the mesh information in the module variables. One
   !! block per group is computed. The Jacobian elements are generated by difference quotients. 
   !! This routine calls a user routine `rblock` to compute the block `(jx,jy)` of \(R\).
      real(rk), intent(in) :: t
         !! Independent variable.
      real(rk), intent(inout) :: u(*)
         !! Current dependent variables.
      real(rk), intent(in) :: r0(*)
         !! Current value of the vector \(R(t,u)\).
      interface
         subroutine rblock(t, jx, jy, uxy, rxy)
         !! User routine that computes a single block of \(R\).
            import :: rk
            real(rk), intent(in) :: t
               !! Current independent variable.
            integer, intent(in) :: jx
               !! Spatial index in x-direction.
            integer, intent(in) :: jy
               !! Spatial index in y-direction.
            real(rk), intent(in) :: uxy(*)
               !! Block of `ns` dependent variables at spatial point `(jx,jy)`.
            real(rk), intent(out) :: rxy(*)
              !! Block `(jx,jy)` of \(R\).
         end subroutine rblock
      end interface
      real(rk), intent(inout) :: r1(*)
        !! Real work space available to this subroutine.
      real(rk), intent(in) :: rewt(*)
         !! Reciprocal error weights.
      real(rk), intent(in) :: cj
         !! Scalar used in forming the system Jacobian.
      real(rk), intent(out) :: bd(*)
         !! LU factors of the diagonal blocks.
      integer, intent(out) :: ipbd(*)
         !! Pivots for the LU factorizations.
      integer, intent(out) :: ierr
         !! Error flag. If no error occurred, `ierr = 0`; if a zero pivot was found at the k-th
         !! stage in one of the LU factorizations, this routine returns `ierr = k > 0`.

      external :: dgefa
      integer :: i, ibd, idiag, ig, igx, igy, iip, j, j0, j00, js, jx, jy
      real(rk) :: del, dfac, fac, uj

      ! Make MP calls to RBLOCK to approximate each diagonal block of dR/du.
      dfac = 1e-2_rk
      ibd = 0
      do igy = 1, ngy
         jy = jyr(igy)
         j00 = (jy - 1)*mxmp
         do igx = 1, ngx
            jx = jxr(igx)
            j0 = j00 + (jx - 1)*mp
            ! If R0 has not been set previously as an array of length NEQ, it can
            ! be set here, as an array of length MP, with the call
            !    CALL RBLOCK (T, JX, JY, U(J0+1), R0)
            ! In this case, change R0(J0+I) below to R0(I).
            do js = 1, mp
               j = j0 + js
               uj = u(j)
               del = max(srur*abs(uj), dfac/rewt(j))
               u(j) = u(j) + del
               fac = -one/del
               call rblock(t, jx, jy, u(j0 + 1), r1)
               do i = 1, mp
                  bd(ibd + i) = (r1(i) - r0(j0 + i))*fac
               end do
               u(j) = uj
               ibd = ibd + mp
            end do
         end do
      end do

      ! Add matrix CJ*Id, and do LU decomposition on blocks.
      ibd = 1
      iip = 1
      do ig = 1, ngrp
         idiag = ibd
         do i = 1, mp
            if (i <= mpd) bd(idiag) = bd(idiag) + cj
            idiag = idiag + (mp + 1)
         end do
         call dgefa(bd(ibd), mp, mp, ipbd(iip), ierr)
         if (ierr /= 0) exit
         ibd = ibd + mpsq
         iip = iip + mp
      end do

   end subroutine jac_rbgpre

   subroutine psol_rbgpre(b, bd, ipbd)
   !! This routine solves a linear system \(A_R x = b\), using the LU factors of the diagonal 
   !! blocks computed in [[jac_rbgpre]] and mesh parameters passed as module variables. The
   !! solution is carried out by [[dgesl]] from LINPACK.
      real(rk), intent(inout) :: b(*)
         !! Right-hand side vector on entry and solution vector on return.
      real(rk), intent(in) :: bd(*)
         !! LU factors of the diagonal blocks. `bd` is the segment `rwp` of `rwork`.
      integer, intent(in) :: ipbd(*)
         !! Pivots for the LU factorizations. `ipbd` is the segment `iwp` of `iwork`.

      external :: dgesl
      integer :: ib, ibd, ig0, igm1, igx, igy, iip, jx, jy

      ib = 1
      do jy = 1, meshy
         igy = jigy(jy)
         ig0 = (igy - 1)*ngx
         do jx = 1, meshx
            igx = jigx(jx)
            igm1 = igx - 1 + ig0
            ibd = 1 + igm1*mpsq
            iip = 1 + igm1*mp
            call dgesl(bd(ibd), mp, mp, ipbd(iip), b(ib), 0)
            ib = ib + mp
         end do
      end do

   end subroutine psol_rbgpre

end module daskr_rbgpre
