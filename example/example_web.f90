!----------------------------------------------------------------------------------------------
! Adapted from original Fortran code in `original/examples/dweb.f`
!----------------------------------------------------------------------------------------------

module web_par
   !! Common set of parameters for [[example_web]].
   
   use daskr_kinds, only: rk, one
   implicit none

   integer, parameter :: np = 1, ns = 2*np, mx = 40, my = mx, mxns = mx*ns
   real(rk), parameter :: pi = 4*atan(one)

   real(rk) :: aa, ee, gg, bb, dprey, dpred
   real(rk) :: acoef(ns, ns), alpha, ax, ay, bcoef(ns), beta, cox(ns), coy(ns), diff(ns), dx, dy

contains

   impure subroutine setpar()
   !! This routine sets the basic problem parameters which are passed to the various routines
   !! via the module [[web_par]].
   
      integer :: i, j

      ax = one
      ay = one

      aa = one
      ee = 1e4_rk
      gg = 0.5e-6_rk
      bb = one
      dprey = one
      dpred = 0.05_rk
      alpha = 5e1_rk
      beta = 3e2_rk

      dx = ax/(mx - 1)
      dy = ay/(my - 1)

      do j = 1, np
         do i = 1, np
            acoef(np + i, j) = ee
            acoef(i, np + j) = -gg
         end do
         acoef(j, j) = -aa
         acoef(np + j, np + j) = -aa
         bcoef(j) = bb
         bcoef(np + j) = -bb
         diff(j) = dprey
         diff(np + j) = dpred
      end do

      do i = 1, ns
         cox(i) = diff(i)/dx**2
         coy(i) = diff(i)/dy**2
      end do

   end subroutine setpar

end module web_par

module web_m
   !! Procedures for [[example_web]].

   use daskr_kinds, only: rk, zero, one
   implicit none
   private

   public :: setid, cinit, output, res, jac, psol, rt, c1_average

contains

   pure subroutine setid(mx, my, ns, nsd, lid, iwork)
   !! This routine sets the ID array in `iwork`, indicating which components are differential
   !! and which are algebraic.

      integer, intent(in) :: mx
      integer, intent(in) :: my
      integer, intent(in) :: ns
      integer, intent(in) :: nsd
      integer, intent(in) :: lid
      integer, intent(inout) :: iwork(*)

      integer :: i, i0, i00, jx, jy, nsdp1

      nsdp1 = nsd + 1
      do jy = 1, my
         i00 = mx*ns*(jy - 1) + lid
         do jx = 1, mx
            i0 = i00 + ns*(jx - 1)
            do i = 1, nsd
               iwork(i0 + i) = 1
            end do
            do i = nsdp1, ns
               iwork(i0 + i) = -1
            end do
         end do
      end do

   end subroutine setid

   pure subroutine cinit(c, cdot, pred_ic, rpar)
   !! This routine computes and loads the vectors of initial values.

      use web_par, only: ax, ay, alpha, beta, dx, dy, ns, np, mx, my, mxns, pi
   
      real(rk), intent(out) :: c(:)
      real(rk), intent(out) :: cdot(:)
      real(rk), intent(in)  :: pred_ic
      real(rk), intent(inout) :: rpar(:)

      integer :: i, ioff, iyoff, jx, jy, npp1
      real(rk) :: argx, argy, fac, t, x, y

      ! Load C
      npp1 = np + 1
      do jy = 1, my
         y = (jy - 1)*dy
         argy = 16*y*y*(ay - y)*(ay - y)
         iyoff = mxns*(jy - 1)
         do jx = 1, mx
            x = (jx - 1)*dx
            argx = 16*x*x*(ax - x)*(ax - x)
            ioff = iyoff + ns*(jx - 1)
            fac = one + alpha*x*y + beta*sin(4*pi*x)*sin(4*pi*y)
            do i = 1, np
               c(ioff + i) = 10.0_rk + i*argx*argy
            end do
            do i = npp1, ns
               c(ioff + i) = pred_ic
            end do
         end do
      end do

      ! Load CDOT
      t = zero
      call f(t, c, cdot, rpar)
      do jy = 1, my
         iyoff = mxns*(jy - 1)
         do jx = 1, mx
            ioff = iyoff + ns*(jx - 1)
            do i = npp1, ns
               cdot(ioff + i) = zero
            end do
         end do
      end do

   end subroutine cinit

   impure subroutine output(t, c, ns, mx, my, lout)
   !! This routine prints the values of the individual species densities at the current time
   !! `t`, to logical unit `lout`.
   
      real(rk), intent(in) :: t
      real(rk), intent(in) :: c(ns, mx, my)
      integer, intent(in) :: ns
      integer, intent(in) :: mx
      integer, intent(in) :: my
      integer, intent(in) :: lout

      logical, save :: first = .true.
      integer :: i, jx, jy

      if (first) then
         write (lout, '(a16)', advance='no') 't'
         do i = 1, ns
            do jy = 1, my
               do jx = 1, mx
                  write (lout, '(1x, a16)', advance='no') &
                     "c("//to_string(i)//","//to_string(jx)//","//to_string(jy)//")"
               end do
            end do
         end do
         write (lout, *) ''
         first = .false.
      end if

      write (lout, '(es16.5e3)', advance="no") t
      do i = 1, ns
         do jy = 1, my
            do jx = 1, mx
               write (lout, '(1x, es16.5e3)', advance="no") c(i, jx, jy)
            end do
         end do
      end do
      write (lout, *) ''

   end subroutine output

   pure subroutine res(t, c, cdot, cj, delta, ires, rpar, ipar)
   !! This routine computes the residuals vector, using subroutine `f` for the right-hand
   !! sides.

      use web_par, only: ns, np, mx, my, mxns
   
      real(rk), intent(in) :: t
      real(rk), intent(in) :: c(*)
      real(rk), intent(in) :: cdot(*)
      real(rk), intent(in) :: cj
      real(rk), intent(out) :: delta(*)
      integer, intent(inout) :: ires
      real(rk), intent(inout) :: rpar(*)
      integer, intent(inout) :: ipar(*)

      integer :: i, ic0, ici, iyoff, jx, jy

      call f(t, c, delta, rpar)

      do jy = 1, my
         iyoff = mxns*(jy - 1)
         do jx = 1, mx
            ic0 = iyoff + ns*(jx - 1)
            do i = 1, ns
               ici = ic0 + i
               if (i > np) then
                  delta(ici) = -delta(ici)
               else
                  delta(ici) = cdot(ici) - delta(ici)
               end if
            end do
         end do
      end do

   end subroutine res

   pure subroutine f(t, c, cdot, rpar)
   !! This routine computes the right-hand sides of all the equations and returns them in the
   !! array `cdot`. The interaction rates are computed by calls to `rates`, and these are saved
   !! in `rpar` for later use in preconditioning.

      use web_par, only: cox, coy, ns, mx, my, mxns
   
      real(rk), intent(in) :: t
      real(rk), intent(in) :: c(*)
      real(rk), intent(out) :: cdot(*)
      real(rk), intent(inout) :: rpar(*)

      integer :: i, ic, ici, idxl, idxu, idyl, idyu, iyoff, jx, jy
      real(rk) :: dcxli, dcxui, dcyli, dcyui

      do jy = 1, my
         iyoff = mxns*(jy - 1)
         idyu = mxns
         if (jy == my) idyu = -mxns
         idyl = mxns
         if (jy == 1) idyl = -mxns
         do jx = 1, mx
            ic = iyoff + ns*(jx - 1) + 1
            ! Get interaction rates at one point (X,Y)
            call rates(t, jx, jy, c(ic), rpar(ic))
            idxu = ns
            if (jx == mx) idxu = -ns
            idxl = ns
            if (jx == 1) idxl = -ns
            do i = 1, ns
               ici = ic + i - 1
               ! Do differencing in Y
               dcyli = c(ici) - c(ici - idyl)
               dcyui = c(ici + idyu) - c(ici)
               ! Do differencing in X
               dcxli = c(ici) - c(ici - idxl)
               dcxui = c(ici + idxu) - c(ici)
               ! Collect terms and load CDOT elements
               cdot(ici) = coy(i)*(dcyui - dcyli) + cox(i)*(dcxui - dcxli) + rpar(ici)
            end do
         end do
      end do

   end subroutine f

   pure subroutine rates(t, jx, jy, c, rxy)
   !! This routine computes one block of the rate term of the system, namely block `(jx, jy)`,
   !! for use in preconditioning.

      use web_par, only: alpha, beta, acoef, bcoef, dx, dy, pi, ns
   
      real(rk), intent(in) :: t
      integer, intent(in) :: jx
      integer, intent(in) :: jy
      real(rk), intent(in) :: c(*)
      real(rk), intent(out) :: rxy(*)

      integer :: i
      real(rk) :: fac, x, y

      y = (jy - 1)*dy
      x = (jx - 1)*dx
      fac = one + alpha*x*y + beta*sin(4*pi*x)*sin(4*pi*y)

      rxy(1:ns) = zero
      do i = 1, ns
         rxy(1:ns) = rxy(1:ns) + c(i)*acoef(1:ns, i)
      end do

      do i = 1, ns
         rxy(i) = c(i)*(bcoef(i)*fac + rxy(i))
      end do

   end subroutine rates

   subroutine jac(res, ires, neq, t, c, cdot, rewt, savr, wk, h, cj, rwp, iwp, ierr, rpar, ipar)
   !! This routine interfaces to subroutines [[jac_rbdpre]] or [[jac_rbgpre]], depending on the 
   !! flag `jbg = ipar(2)`, to generate and preprocess the block-diagonal Jacobian corresponding
   !! to the reaction term.
   !!
   !! * If `jbg == 0`, we call [[jac_rbdpre]], without block-grouping.
   !! * If `jbg == 1`, we call [[jac_rbgpre]], with block-grouping.
   !!
   !! Array `rpar`, containing the current \(R\) vector, is passed to [[jac_rbdpre]] and
   !! [[jac_rbgpre]] as argument `r0`, consistent with the loading of `rpar` in  subroutine `f`.
   !! The procedure name `rates` is passed as the name of the routine which computes the individual
   !! blocks of \(R\).

      use daskr_rbdpre, only: jac_rbdpre
      use daskr_rbgpre, only: jac_rbgpre
      use daskr, only: res_t
      
      procedure(res_t) :: res
      integer, intent(out) :: ires
      integer, intent(in) :: neq
      real(rk), intent(in) :: t
      real(rk), intent(inout) :: c(*)
      real(rk), intent(in) :: cdot(*)
      real(rk), intent(in) :: rewt(*)
      real(rk), intent(inout) :: savr(*)
      real(rk), intent(inout) :: wk(*)
      real(rk), intent(in) :: h
      real(rk), intent(in) :: cj
      real(rk), intent(inout) :: rwp(*)
      integer, intent(inout) :: iwp(*)
      integer, intent(inout) :: ierr
      real(rk), intent(inout) :: rpar(*)
      integer, intent(inout) :: ipar(*)

      integer :: jbg

      jbg = ipar(2)
      if (jbg == 0) then
         call jac_rbdpre(t, c, rpar, rates, wk, rewt, cj, rwp, iwp, ierr)
      else
         call jac_rbgpre(t, c, rpar, rates, wk, rewt, cj, rwp, iwp, ierr)
      end if

   end subroutine jac

   subroutine psol(neq, t, c, cdot, savr, wk, cj, wght, rwp, iwp, b, epslin, ierr, rpar, ipar)
   !! This routine applies the inverse of a product preconditioner matrix to the vector in the
   !! array `b`. Depending on the flag `jpre`, this involves a call to `gauss_seidel`, for the
   !! inverse of the spatial factor, and/or a call to [[psol_rbdpre]] or [[psol_rbgpre]] for the
   !! inverse of the reaction-based factor \(c_J I_d - \partial R / \partial y \).
   !! If `jbg == 0`, the latter factor does not use block-grouping and [[jac_rbdpre]] is called.
   !! Otherwise, if `jbg == 1`, block-grouping is used and [[jac_rbgpre]] is called. The array
   !! `b` is overwritten with the solution.

      use daskr_rbdpre, only: psol_rbdpre
      use daskr_rbgpre, only: psol_rbgpre
      
      integer, intent(in) :: neq
      real(rk), intent(in) :: t
      real(rk), intent(in) :: c(*)
      real(rk), intent(in) :: cdot(*)
      real(rk), intent(in) :: savr(*)
      real(rk), intent(inout) :: wk(*)
      real(rk), intent(in) :: cj
      real(rk), intent(in) :: wght(*)
      real(rk), intent(inout) :: rwp(*)
      integer, intent(inout) :: iwp(*)
      real(rk), intent(inout) :: b(*)
      real(rk), intent(in) :: epslin
      integer, intent(inout) :: ierr
      real(rk), intent(inout) :: rpar(*)
      integer, intent(inout) :: ipar(*)

      integer :: jbg, jpre
      real(rk) :: hl0

      jpre = ipar(1)
      jbg = ipar(2)
      ierr = 0
      hl0 = one/cj     

      if (jpre == 2 .or. jpre == 3) call gauss_seidel(neq, hl0, b, wk)

      if (jpre /= 2) then
         if (jbg == 0) call psol_rbdpre(b, rwp, iwp)
         if (jbg == 1) call psol_rbgpre(b, rwp, iwp)
      end if

      if (jpre == 4) call gauss_seidel(neq, hl0, b, wk)

   end subroutine psol

   pure subroutine gauss_seidel(n, hl0, z, x)
   !! This routine provides the inverse of the spatial factor for a product preconditoner in an
   !! ns-species reaction-diffusion problem. It performs `itmax` Gauss-Seidel iterations to
   !! compute an approximation to \(A_S^{-1} z\), where \(A_S = I - h_{l0} J_d \), and \(J_d\)
   !! represents the diffusion contributions to the Jacobian. The solution vector is returned
   !! in `z`.

      use web_par, only: cox, coy, ns, mx, my, mxns
   
      integer, intent(in) :: n
      real(rk), intent(in) :: hl0
      real(rk), intent(inout) :: z(*)
      real(rk), intent(inout) :: x(*)

      integer, parameter :: itmax = 5
      integer :: i, ic, ici, iter, iyoff, jx, jy
      real(rk) :: elamda, beta1(ns), gamma1(ns), beta2(ns), gamma2(ns), dinv(ns)

      ! Write matrix as A = D - L - U.
      ! Load local arrays BETA, BETA2, GAMMA1, GAMMA2, and DINV.
      do i = 1, ns
         elamda = one/(one + 2*hl0*(cox(i) + coy(i)))
         beta1(i) = hl0*cox(i)*elamda
         beta2(i) = 2*beta1(i)
         gamma1(i) = hl0*coy(i)*elamda
         gamma2(i) = 2*gamma1(i)
         dinv(i) = elamda
      end do

      ! Zero X in all its components, since X is added to Z at the end.
      x(1:n) = zero

      ! Load array X with (D-inverse)*Z for first iteration.
      do jy = 1, my
         iyoff = mxns*(jy - 1)
         do jx = 1, mx
            ic = iyoff + ns*(jx - 1)
            do i = 1, ns
               ici = ic + i
               x(ici) = dinv(i)*z(ici)
               z(ici) = zero
            end do
         end do
      end do

      do iter = 1, itmax

         ! Calculate (D-inverse)*U*X
         if (iter > 1) then

            jy = 1
            jx = 1
            ic = ns*(jx - 1)
            do i = 1, ns
               ici = ic + i
               x(ici) = beta2(i)*x(ici + ns) + gamma2(i)*x(ici + mxns)
            end do
            do jx = 2, mx - 1
               ic = ns*(jx - 1)
               do i = 1, ns
                  ici = ic + i
                  x(ici) = beta1(i)*x(ici + ns) + gamma2(i)*x(ici + mxns)
               end do
            end do
            jx = mx
            ic = ns*(jx - 1)
            do i = 1, ns
               ici = ic + i
               x(ici) = gamma2(i)*x(ici + mxns)
            end do
            do jy = 2, my - 1
               iyoff = mxns*(jy - 1)
               jx = 1
               ic = iyoff
               do i = 1, ns
                  ici = ic + i
                  x(ici) = beta2(i)*x(ici + ns) + gamma1(i)*x(ici + mxns)
               end do
               do jx = 2, mx - 1
                  ic = iyoff + ns*(jx - 1)
                  do i = 1, ns
                     ici = ic + i
                     x(ici) = beta1(i)*x(ici + ns) + gamma1(i)*x(ici + mxns)
                  end do
               end do
               jx = mx
               ic = iyoff + ns*(jx - 1)
               do i = 1, ns
                  ici = ic + i
                  x(ici) = gamma1(i)*x(ici + mxns)
               end do
            end do
            jy = my
            iyoff = mxns*(jy - 1)
            jx = 1
            ic = iyoff
            do i = 1, ns
               ici = ic + i
               x(ici) = beta2(i)*x(ici + ns)
            end do
            do jx = 2, mx - 1
               ic = iyoff + ns*(jx - 1)
               do i = 1, ns
                  ici = ic + i
                  x(ici) = beta1(i)*x(ici + ns)
               end do
            end do
            jx = mx
            ic = iyoff + ns*(jx - 1)
            do i = 1, ns
               ici = ic + i
               x(ici) = zero
            end do

         end if

         ! Calculate [(I - (D-inverse)*L)]-inverse * X
         jy = 1
         do jx = 2, mx - 1
            ic = ns*(jx - 1)
            do i = 1, ns
               ici = ic + i
               x(ici) = x(ici) + beta1(i)*x(ici - ns)
            end do
         end do
         jx = mx
         ic = ns*(jx - 1)
         do i = 1, ns
            ici = ic + i
            x(ici) = x(ici) + beta2(i)*x(ici - ns)
         end do
         do jy = 2, my - 1
            iyoff = mxns*(jy - 1)
            jx = 1
            ic = iyoff
            do i = 1, ns
               ici = ic + i
               x(ici) = x(ici) + gamma1(i)*x(ici - mxns)
            end do
            do jx = 2, mx - 1
               ic = iyoff + ns*(jx - 1)
               do i = 1, ns
                  ici = ic + i
                  x(ici) = (x(ici) + beta1(i)*x(ici - ns)) + gamma1(i)*x(ici - mxns)
               end do
            end do
            jx = mx
            ic = iyoff + ns*(jx - 1)
            do i = 1, ns
               ici = ic + i
               x(ici) = (x(ici) + beta2(i)*x(ici - ns)) + gamma1(i)*x(ici - mxns)
            end do
         end do
         jy = my
         iyoff = mxns*(jy - 1)
         jx = 1
         ic = iyoff
         do i = 1, ns
            ici = ic + i
            x(ici) = x(ici) + gamma2(i)*x(ici - mxns)
         end do
         do jx = 2, mx - 1
            ic = iyoff + ns*(jx - 1)
            do i = 1, ns
               ici = ic + i
               x(ici) = (x(ici) + beta1(i)*x(ici - ns)) + gamma2(i)*x(ici - mxns)
            end do
         end do
         jx = mx
         ic = iyoff + ns*(jx - 1)
         do i = 1, ns
            ici = ic + i
            x(ici) = (x(ici) + beta2(i)*x(ici - ns)) + gamma2(i)*x(ici - mxns)
         end do

         ! Add increment X to Z
         do i = 1, n
            z(i) = z(i) + x(i)
         end do

      end do

   end subroutine gauss_seidel

   pure subroutine c1_average(c, c1ave)
   !! This routine computes the spatial average value of \(c_1\).

      use web_par, only: mx, my, ns, mxns
      real(rk), intent(in) :: c(*)
      real(rk), intent(out) :: c1ave

      integer :: ioff, iyoff, jx, jy
      real(rk) :: total

      total = zero
      do jy = 1, my
         iyoff = mxns*(jy - 1)
         do jx = 1, mx
            ioff = iyoff + ns*(jx - 1)
            total = total + c(ioff + 1)
         end do
      end do

      c1ave = total/(mx*my)

   end subroutine c1_average

   pure subroutine rt(neq, t, c, cdot, nrt, rval, rpar, ipar)
   !! Roots routine.

      integer, intent(in) :: neq
      real(rk), intent(in) :: t
      real(rk), intent(in) :: c(neq)
      real(rk), intent(in) :: cdot(neq)
      integer, intent(in) :: nrt
      real(rk), intent(out) :: rval(nrt)
      real(rk), intent(in) :: rpar(*)
      integer, intent(in) :: ipar(*)

      real(rk) :: c1ave

      call c1_average(c, c1ave)
      rval(1) = c1ave - 20.0_rk

   end subroutine rt

   function to_string(value) result(string)
      integer, intent(in) :: value
      character(len=:), allocatable :: string
      character(len=128) :: buffer

      write (buffer, *) value
      string = trim(adjustl(buffer))

   end function

end module web_m

program example_web
!! Example program for [[daskr]]:
!! DAE system derived from \(s\)-species interaction PDE in 2 dimensions.
!!
!! This program solves a DAE system that arises from a system of partial differential equations.
!! The PDE system is a food web population model, with predator-prey interaction and diffusion
!! on the unit square in two dimensions. The dependent variable vector is
!!
!! $$ c = [c_1, c_2 , ..., c_s] $$
!!
!! and the PDEs are as follows:
!!
!! $$\begin{aligned}
!! \frac{\partial c_i}{\partial t} &= d_i \left( \frac{\partial^2 c_i}{\partial x^2}
!!           + \frac{\partial^2 c_i}{\partial y^2} \right) + R_i \quad i=1,...,s/2 \\
!!                               0 &= d_i \left( \frac{\partial^2 c_i}{\partial x^2}
!!           + \frac{\partial^2 c_i}{\partial y^2} \right) + R_i \quad i=s/2+1,...,s
!! \end{aligned}$$
!!
!! where the rate of formation of species \(i\) is given by:
!!
!! $$ R_i(x,y,c) = c_i \left( b_i + \sum_{j=1}^s a_{ij} c_j \right) $$
!!
!! The number of species is \(s = 2 p\), with the first \(p\) being prey and the last \(p\)
!! being predators. The coefficients \(a_{ij}\), \(b_i\), \(d_i\) are:
!!
!! $$\begin{aligned}
!!   a_{ii} &= -a  \quad (\mathrm{all}\; i) \\
!!   a_{ij} &= -g  \quad (i \le p,\; j > p) \\
!!   a_{ij} &=  e  \quad (i > p,\; j \le p)
!! \end{aligned}$$
!!
!! $$\begin{aligned}
!!   b_i    &=  b (1 + \alpha x y + \beta \sin(4 \pi x) \sin(4 \pi y)) \quad (i \le p) \\
!!   b_i    &= -b (1 + \alpha x y + \beta \sin(4 \pi x) \sin(4 \pi y)) \quad (i > p)
!! \end{aligned}$$
!!
!! $$\begin{aligned}
!!   d_i    &= d_{prey} \quad (i \le p) \\
!!   d_i    &= d_{pred} \quad (i > p)
!! \end{aligned}$$
!!
!! The various scalar parameters are set in subroutine [[setpar]].
!!
!! The boundary conditions are of Neumann type (zero normal derivative) everywhere. A polynomial
!! in \(x\) and \(y\) is used to set the initial conditions. The PDEs are discretized by central
!! differencing on a \( M_x \times M_y\) mesh.
!!
!! The root function is:
!!
!! $$ r(t,c,\dot{c}) = -20 + \iint c_1 dx dy $$
!!
!! The DAE system is solved by [[daskr]] with three different method options:
!!
!! 1. Direct band method for the linear systems (internal Jacobian),
!! 2. Preconditioned Krylov method for the linear systems, without block-grouping in the
!!    reaction-based factor, and
!! 3. Preconditioned Krylov method for the linear systems, with block-grouping in the
!!    reaction-based factor.
!!
!! In the Krylov cases, the preconditioner is the product of two factors:
!!
!! * The spatial factor uses a fixed number of Gauss-Seidel iterations based on the diffusion
!!   terms only.
!! * The reaction-based factor is a block-diagonal matrix based on the partial derivatives of
!!   the interaction terms \(R\) only.
!!
!! With block-grouping, only a subset of the \((s \times s)\) blocks are computed. An integer
!! flag, `jpre`, is set in the main program to specify whether the preconditioner is to use only
!! one of the two factors or both, and in which order.
!!
!! The reaction-based preconditioner factor is set up and solved in independent modules:
!!
!! * [[daskr_rbdpre]] in the case of no block-grouping, and
!! * [[daskr_rbgpre]] in the case of block-grouping.
!!
!! These modules are provided separately for general use on reaction-transport problems.
!!
!! Two output files are written: one with the problem description and performance statistics on
!! and one with solution profiles at selected output times. The solution file is written only
!! in the case of the direct method.
!!
!! **References**
!!
!! 1. Peter N. Brown and Alan C. Hindmarsh, "Reduced Storage Matrix Methods in Stiff ODE
!!    Systems", J. Appl. Math. & Comp., 31 (1989), pp. 40-91.
!! 2. Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold, "Using Krylov Methods in the
!!    Solution of Large-Scale Differential-Algebraic Systems", SIAM J. Sci. Comput., 15 (1994),
!!    pp. 1467-1488.
!! 3. Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold, "Consistent Initial Condition
!!    Calculation for Differential-Algebraic Systems", LLNL Report UCRL-JC-122175, August 1995;
!!    SIAM J. Sci. Comput., 19 (1998), pp. 1495 - 1512.

   use daskr_kinds, only: rk, zero, one
   use daskr_rbdpre, only: setup_rbdpre, jac_rbdpre, psol_rbdpre
   use daskr_rbgpre, only: setup_rbgpre, jac_rbgpre, psol_rbgpre
   use web_par, only: aa, alpha, bb, beta, dpred, dprey, ee, gg, mx, my, mxns, np, ns, setpar
   use web_m, only: setid, cinit, output, res, jac, psol, rt, c1_average
   implicit none

   integer, parameter :: neq = ns*mx*my, maxm = max(mx, my)
   integer, parameter :: lrw = 63 + (3*ns*maxm + 11)*neq + maxm, liw = 40 + 2*neq
   integer :: idid, info(20), iout, ipar(2), iwork(liw), jbg, jpre, jroot, leniw, &
              lenrw, ldout, lcout, meth, ncfl, ncfn, ng, nli, nlidif, nni, nnidif, nout, &
              npe, nps, nqu, nre, nrt, nrte, nst, nxg, nyg
   real(rk) :: atol, avlin, c1ave, c(neq), cdot(neq), hu, pred_ic, rpar(neq), rtol, &
               rwork(lrw), t, tout

   ! Dimension solution arrays and work arrays.
   !
   ! When INFO(12) = 0, with INFO(5) = 0, INFO(6) = 1:
   !       The length required for RWORK is
   !       60 + 3*NRT + (2*ML+MU+11)*NEQ + 2*(NEQ/(ML+MU+1) + 1) .
   !       For MX = MY = (even number) and ML = MU = NS*MX, this length is
   !       60 + 3*NRT + (3*NS*MX + 11)*NEQ + MY .
   !       The length required for IWORK is  40 + 2*NEQ .
   !
   ! When INFO(12) = 1:
   !       The length required for RWORK is
   !       101 + 3*NRT + 19*NEQ + LENWP = 104 + 19*NEQ + NS*NS*NGRP .
   !       The length required for IWORK is
   !       40 + NEQ + LENIWP = 40 + NEQ + NS*NGRP .
   !
   ! The dimensions for the various arrays are set using parameters
   !       MAXN    which must be >= NEQ = NS*MX*MY,
   !       MAXM    which must be >= MAX(MX,MY).

   ! Open output files.
   open (newunit=ldout, file='./example/example_web_d.out', action="write", position="rewind")
   open (newunit=lcout, file='./example/example_web_c.out', action="write", position="rewind")

   ! Call SETPAR to set basic problem parameters.
   call setpar

   ! Set NRT = number of root functions.
   nrt = 1

   write (ldout, '(a, /)') 'DWEB: Example program for DASKR package'
   write (ldout, '(a, i4)') 'Food web problem with NS species, NS =', ns
   write (ldout, '(a, /)') 'Predator-prey interaction and diffusion on a 2-D square'
   write (ldout, '(a, e12.4, a, e12.4, a, e12.4)') &
      'Matrix parameters..  a =', aa, '  e =', ee, '  g =', gg
   write (ldout, '(21x, a, e12.4)') &
      'b =', bb
   write (ldout, '(a, e12.4, a, e12.4)') &
      'Diffusion coefficients: dprey =', dprey, '  dpred =', dpred
   write (ldout, '(a, e12.4, a, e12.4)') &
      'Rate parameters alphaa =', alpha, ' and beta =', beta
   write (ldout, '(a, 2i4, 5x, a, i7)') &
      'Mesh dimensions (MX,MY) =', mx, my, ' Total system size is NEQ =', neq
   write (ldout, '(a)') 'Root function is R(Y) = average(c1) - 20'

   ! Set the flat initial guess for the predators.
   pred_ic = 1e5_rk

   ! Set remaining method parameters for DASKR.
   ! These include the INFO array and tolerances.
   info = 0

   ! Set INFO(11) = 1, indicating I.C. calculation requested.
   info(11) = 1

   ! Set INFO(14) = 1 to get the computed initial values.
   info(14) = 1

   ! Set INFO(15) = 1 to signal that a preconditioner setup routine is to be called in the
   ! Krylov case.
   info(15) = 1

   ! Set INFO(16) = 1 to get alternative error test (on the differential variables only).
   info(16) = 1

   ! Set the tolerances.
   rtol = 1e-5_rk
   atol = rtol

   write (ldout, '(/, a, e10.2, a, e10.2)') &
      'Tolerance parameters: RTOL =', rtol, '  ATOL =', atol
   write (ldout, '(a, i2, a)') &
      'Internal I.C. calculation flag INFO(11) =', info(11), '  (0 = off, 1 = on)'
   write (ldout, '(a, e10.2)') &
      'Predator I.C. guess =', pred_ic
   write (ldout, '(a, i2, a)') &
      'Alternate error test flag INFO(16) =', info(16), '  (0 = off, 1 = on)'

   ! Set NOUT = number of output times.
   nout = 18

   ! Loop over method options:
   ! METH = 0 means use INFO(12) = 0 (direct)
   ! METH = 1 means use INFO(12) = 1 (Krylov) without block-grouping in the reaction-based
   !          factor in the preconditioner.
   ! METH = 2 means use INFO(12) = 1 (Krylov) with block-grouping in the reaction-based factor
   !          in the preconditioner.
   ! A block-grouping flag JBG, communicated through IPAR, is set to 0 (no block-grouping) or
   ! 1 (use block-grouping) with METH = 1 or 2.
   ! Reset INFO(1) = 0 and INFO(11) = 1.

   do meth = 0, 2

      info(12) = min(meth, 1)
      info(1) = 0
      info(11) = 1
      jbg = meth - 1
      ipar(2) = jbg

      write (ldout, '(/, 80("."), //, a, i2, a)') &
         'Linear solver method flag INFO(12) =', info(12), '  (0 = direct, 1 = Krylov)'

      ! In the case of the direct method, set INFO(6) = 1 to signal a banded Jacobian, set
      ! IWORK(1) = IWORK(2) = MX*NS, the half-bandwidth, and call SETID to set the IWORK
      ! segment containing the block ID indicating the differential and algebraic components.
      if (info(12) == 0) then
         info(6) = 1
         iwork(1) = mxns
         iwork(2) = mxns
         call setid(mx, my, ns, np, 40, iwork)
         write (ldout, '(a, i4)') 'Difference-quotient banded Jacobian, half-bandwidths =', mxns
      end if

      ! In the case of the Krylov method, set and print various preconditioner parameters.
      if (info(12) == 1) then

         ! First set the preconditioner choice JPRE.
         !  JPRE = 1 means reaction-only (block-diagonal) factor A_R
         !  JPRE = 2 means spatial factor (Gauss-Seidel) A_S
         !  JPRE = 3 means A_S * A_R
         !  JPRE = 4 means A_R * A_S
         ! Use IPAR to communicate JPRE to the preconditioner solve routine.
         jpre = 3
         ipar(1) = jpre
         write (ldout, '(a, i3)') 'Preconditioner flag is JPRE =', jpre
         write (ldout, '(a)') &
            '(1 = reaction factor A_R, 2 = spatial factor A_S, 3 = A_S*A_R, 4 = A_R*A_S)'

         ! Call the preconditioner setup routines to set the 2D mesh parameters and 
         ! block-grouping data, and the IWORK segment ID indicating the differential and 
         ! algebraic components.
         if (jbg == 0) then
            call setup_rbdpre(mx, my, ns, np, 40, iwork)
            write (ldout, '(a)') 'No block-grouping in reaction factor'
         end if
         if (jbg == 1) then
            nxg = 5
            nyg = 5
            ng = nxg*nyg
            call setup_rbgpre(mx, my, ns, np, nxg, nyg, 40, iwork)
            write (ldout, '(a)') 'Block-grouping in reaction factor'
            write (ldout, '(a, i5, a, i3, a, i3, a)') &
               'Number of groups =', ng, ' (NGX by NGY, NGX =', nxg, ',  NGY =', nyg, ')'
         end if
      end if

      ! Set the initial T and TOUT, and call CINIT to set initial values.
      t = zero
      tout = 1.0e-8_rk
      call cinit(c, cdot, pred_ic, rpar)

      nli = 0
      nni = 0

      ! Header for "ldout"
      write (ldout, '(/10x, a)') &
         't      <c1>  NSTEP   NRE   NNI   NLI   NPE    NQ          H    AVLIN'

      ! Loop over output times and call DASKR. At each output time, print average c1 value and
      ! performance data.
      ! The first call, with IOUT = 0, is to calculate initial values only.
      ! After the first call, reset INFO(11) = 0 and the initial TOUT.
      ! If a root was found, we flag this, and return to the DASKR call.
      do iout = 0, nout

         do
            call daskr(res, neq, t, c, cdot, tout, info, rtol, atol, &
                       idid, rwork, lrw, iwork, liw, rpar, ipar, jac, psol, &
                       rt, nrt, jroot)

            nst = iwork(11)
            nre = iwork(12)
            npe = iwork(13)
            nnidif = iwork(19) - nni
            nni = iwork(19)
            nlidif = iwork(20) - nli
            nli = iwork(20)
            nqu = iwork(8)
            hu = rwork(7)
            avlin = zero
            if (nnidif > 0) avlin = one*nlidif/nnidif

            if (meth == 0) then
               call output(t, c, ns, mx, my, lcout)
            end if

            call c1_average(c, c1ave)
            write (ldout, '(e11.5, f10.5, i7, i6, i6, i6, i6, i6, e11.2, f9.4)') &
               t, c1ave, nst, nre, nni, nli, npe, nqu, hu, avlin

            if (idid == 5) then
               write (ldout, '(15x, a, i3)') '*****   Root found, JROOT =', jroot
            else
               exit
            end if

         end do

         if (idid < 0) then
            write (ldout, '(//a, e12.4//)') 'Final time reached =', t
            exit
         end if

         if (tout > 0.9_rk) then
            tout = tout + one
         else
            tout = 10*tout
         end if

         if (iout == 0) then
            info(11) = 0
            tout = 1.0e-8_rk
            nli = 0
            nni = 0
         end if

      end do

      lenrw = iwork(18)
      leniw = iwork(17)
      nst = iwork(11)
      nre = iwork(12)
      npe = iwork(13)
      nni = iwork(19)
      nli = iwork(20)
      nps = iwork(21)
      if (nni > 0) avlin = real(nli)/real(nni)
      ncfn = iwork(15)
      ncfl = iwork(16)
      nrte = iwork(36)

      write (ldout, '(//a)') 'Final statistics for this run..'
      write (ldout, '(a, i8, a, i6)') 'RWORK size =', lenrw, '  IWORK size =', leniw
      write (ldout, '(a, i5)') 'Number of time steps              =', nst
      write (ldout, '(a, i5)') 'Number of residual evaluations    =', nre
      write (ldout, '(a, i5)') 'Number of root fn. evaluations    =', nrte
      write (ldout, '(a, i5)') 'Number of Jac. or prec. evals.    =', npe
      write (ldout, '(a, i5)') 'Number of preconditioner solves   =', nps
      write (ldout, '(a, i5)') 'Number of nonlinear iterations    =', nni
      write (ldout, '(a, i5)') 'Number of linear iterations       =', nli
      write (ldout, '(a, f8.4)') 'Average Krylov subspace dimension =', avlin
      write (ldout, '(i3, a, i5, a)') ncfn, ' nonlinear conv. failures,', ncfl, ' linear conv. failures'

   end do

   close (unit=ldout)
   close (unit=lcout)

end program example_web
