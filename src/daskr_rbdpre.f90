!----------------------------------------------------------------------------------------------
! Adapted from original Fortran code in `original/preconds/drbdpre.f`
!----------------------------------------------------------------------------------------------

module daskr_rbdpre
!! # Preconditioner Tools for Reaction-Transport Problems. Part I: Block-Diagonal Reaction-Based Factor without Grouping
!!
!! This module is provided to assist in the generation and solution of preconditioner matrices
!! for problems arising from reaction-transport systems. More specifically, the routines
!! [[jac_rbdpre]], and [[psol_rbdpre]] are intended as *auxiliary* routines for the user-supplied
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
!! that evaluates the reaction terms at a single spatial point. \(A_R\) has one such block for
!! each spatial point in the mesh. For a more economical approximation, see [[daskr_rbgpre]]
!! which uses block-grouping in \(A_R\).
!!
!! The routines given here are specialized to the case of a 2D problem on a rectangular mesh in 
!! the x-y plane. However, they can be easily modified for a different problem geometry.  It is
!! also assumed that the PDE variables are ordered so that the differential variables appear
!! first, followed by the algebraic variables.
!!
!! ## Usage
!!
!! @warning
!! This module uses module variables to share the mesh and block-group parameters among routines
!! and is, thus, **thread unsafe**.  
!!           
!! To use these routines in conjunction with [[daskr]], the user's calling program should include
!! the following, in addition to setting the other [[daskr]] input parameters:
!!
!! * Set `info(12) = 1` to select the Krylov iterative method and `info(15) = 1` to indicate
!!   that a `jac` routine exists.
!!         
!! * Call [[setup_rbdpre]] to set mesh-related inputs. 
!!      
!! * A `jac` routine, as prescribed by the [[daskr]] instructions, which calls [[jac_rbdpre]],
!!   and does any other Jacobian-related preprocessing needed for preconditioning. The latter
!!   routine takes as argument `r0` the current value of the reaction vector \(R\). This can be
!!   done, for example, by taking `r0` to be `rpar`, and loading `rpar` with the vector \(R\)
!!   in the last call to the `res` routine; in that case, the calling program must declare `rpar`
!!   to have length at least `neq`. Alternatively, insert a call to `rblock` within the loop 
!!   over mesh points inside [[jac_rbdpre]].
!!   
!! * A `psol` routine, as prescribed by the [[daskr]] instructions, which calls [[psol_rbdpre]]
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

   use daskr_kinds, only: rk, one
   implicit none
   private

   ! Module variables >>Thread Unsafe<<
   integer :: mp, mpd, mpsq, meshx, meshy, mxmp
   real(rk) :: srur

   public :: setup_rbdpre, jac_rbdpre, psol_rbdpre

contains

   subroutine setup_rbdpre(mx, my, ns, nsd, lid, iwork)
   !! This routine sets the mesh parameters needed to use the routines [[jac_rbdpre]] and 
   !! [[psol_rbdpre]], assuming a 2D rectangular problem. Additionally, it loads the lengths
   !! `lrwp` and `liwp` in array `iwork`.

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
      integer, intent(in) :: lid
         !! Flag indicating whether to load the ID array in `iwork`. If `lid > 0`, set the ID 
         !! array in `iwork`, indicating which components are differential and which are algebraic.
         !! This value is required if either `info(11) = 1` or `info(16) = 1`, in which case
         !! set `lid = 40` or `lid = 40 + neq`, depending on the value of the constraint option
         !! `info(10)`. Otherwise, set `lid = 0`.
      integer, intent(inout) :: iwork(*)
         !! Integer work array.

      integer :: i, i0, jx, jy

      ! Load the modules variables.
      srur = sqrt(epsilon(one))
      mp = ns
      mpd = nsd
      mpsq = ns*ns
      meshx = mx
      meshy = my
      mxmp = meshx*mp

      ! Set the sizes of the preconditioning storage space segments in RWORK and IWORK.
      iwork(27) = mpsq*meshx*meshy
      iwork(28) = mp*meshx*meshy

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
         end do
         i0 = i0 + mp
         end do
      end do

   end subroutine setup_rbdpre

   subroutine jac_rbdpre(t, u, r0, rblock, r1, rewt, cj, bd, ipbd, ierr)
   !! This routine generates and preprocesses a block-diagonal preconditioner matrix, based on
   !! the part of the Jacobian corresponding to the reaction terms \(R\) of the problem. It 
   !! generates a matrix of the form \(c_J I_d - \partial R/ \partial u\). It calls [[dgefa]]
   !! from LINPACK to do the LU decomposition of each diagonal block. The computation of the
   !! diagonal blocks uses the mesh information in the module variables. One block per spatial
   !! point is computed. The Jacobian elements are generated by difference quotients. 
   !! This routine calls a user routine `rblock` to compute the block `(jx,jy)` of \(R\).

      use dlinpack, only: dgefa

      real(rk), intent(in) :: t
         !! Current independent variable.
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
         !! Error flag.
         !! `0`: no error occurred;
         !! `k > 0`: a zero pivot was found at the `k`-th stage in one of the LU factorizations.
      
      real(rk) :: del, dfac, fac, uj
      integer :: i, ibd, idiag, iip, j, j0, js, jx, jy

      ! Make MP calls to RBLOCK to approximate each diagonal block of dR/du.
      dfac = 1e-2_rk
      ibd = 0
      j0 = 0
      do jy = 1, meshy
         do jx = 1, meshx
            ! If R0 has not been set previously as an array of length NEQ, it can
            ! be set here, as an array of length MP, with the call:
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
            j0 = j0 + mp
         end do
      end do

      ! Add matrix CJ*Id, and do LU decomposition on blocks.
      ibd = 1
      iip = 1
      do j = 1, meshx*meshy
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

   end subroutine jac_rbdpre

   subroutine psol_rbdpre(b, bd, ipbd)
   !! This routine solves a linear system \(A_R x = b\), using the LU factors of the diagonal 
   !! blocks computed in [[jac_rbdpre]] and mesh parameters passed as module variables. The
   !! solution is carried out by [[dgesl]] from LINPACK.

      use dlinpack, only: dgesl

      real(rk), intent(inout) :: b(*)
         !! Right-hand side vector on entry and solution vector on return.
      real(rk), intent(in) :: bd(*)
         !! LU factors of the diagonal blocks. `bd` corresponds to the segment `rwp` of `rwork`.
      integer, intent(in) :: ipbd(*)
         !! Pivots for the LU factorizations. `ipbd` corresponds to the segment `iwp` of `iwork`.

      integer :: ib, ibd, jx, jy

      ib = 1
      ibd = 1
      do jy = 1, meshy
         do jx = 1, meshx
            call dgesl(bd(ibd), mp, mp, ipbd(ib), b(ib), 0)
            ib = ib + mp
            ibd = ibd + mpsq
         end do
      end do

   end subroutine psol_rbdpre

end module daskr_rbdpre
