!----------------------------------------------------------------------------------------------
! Adapted from original Fortran code in `original/solver/ddaskr.f`
!----------------------------------------------------------------------------------------------
module daskr

   use daskr_kinds, only: rk
   implicit none
   
   abstract interface
      subroutine res_t(t, y, ydot, cj, delta, ires, rpar, ipar)
         import :: rk
         real(rk), intent(in) :: t
         real(rk), intent(in) :: y(*)
         real(rk), intent(in) :: ydot(*)
         real(rk), intent(in) :: cj
         real(rk), intent(out) :: delta(*)
         integer, intent(inout) :: ires
         real(rk), intent(inout) :: rpar(*)
         integer, intent(inout) :: ipar(*)
      end subroutine res_t

      subroutine jacd_t(t, y, ydot, pd, cj, rpar, ipar)
         import :: rk
         procedure(res_t) :: res
         real(rk), intent(in) :: t
         real(rk), intent(in) :: y(*)
         real(rk), intent(in) :: ydot(*)
         real(rk), intent(out) :: pd(:,:)
         real(rk), intent(in) :: cj
         real(rk), intent(inout) :: rpar(*)
         integer, intent(inout) :: ipar(*)
      end subroutine jacd_t

      subroutine jack_t(res, ires, neq, t, y, ydot, rewt, savr, wk, h, cj, rwp, iwp, ierr, rpar, ipar)
         import :: rk
         procedure(res_t) :: res
         integer, intent(out) :: ires
         integer, intent(in) :: neq
         real(rk), intent(in) :: t
         real(rk), intent(inout) :: y(*)
         real(rk), intent(inout) :: ydot(*)
         real(rk), intent(in) :: rewt(*)
         real(rk), intent(inout) :: savr(*)
         real(rk), intent(inout) :: wk(*) !(neq)
         real(rk), intent(in) :: h
         real(rk), intent(in) :: cj
         real(rk), intent(inout) :: rwp(*)
         integer, intent(inout) :: iwp(*)
         integer, intent(inout) :: ierr
         real(rk), intent(inout) :: rpar(*)
         integer, intent(inout) :: ipar(*)
      end subroutine jack_t

      subroutine psol_t(neq, t, y, ydot, savr, wk, cj, wght, rwp, iwp, b, epslin, ierr, rpar, ipar)
         import :: rk
         integer, intent(in) :: neq
         real(rk), intent(in) :: t
         real(rk), intent(in) :: y(*)
         real(rk), intent(in) :: ydot(*)
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
      end subroutine psol_t
   end interface

end module daskr

subroutine dlinsk( &
   neq, y, t, ydot, savr, cj, tscale, p, pnorm, &
   wght, sqrtn, rsqrtn, lsoff, stptol, iret, res, ires, psol, &
   rwm, iwm, rhok, fnorm, icopt, id, rwp, iwp, r, epslin, ynew, ydotnew, &
   pwk, icnflg, icnstr, rlx, rpar, ipar)
!! This routine uses a linesearch algorithm to calculate a new \( (y,\dot{y}) \) pair,
!! denoted \( (y_{new},\dot{y}_{new}) \), such that:
!!
!!   $$ f(y_{new},\dot{y}_{new}) \leq (1 - 2 \alpha \lambda) f(y,\dot{y})$$
!!
!! where \(0 < \lambda \le 1\), and \( f(y,\dot{y}) \) is defined as:
!!
!! $$ f(y,\dot{y}) = (1/2) \| (P^{-1} G(t,y,\dot{y}) \|^2 $$
!!
!! where the norm is the weighted root mean square vector norm, \(G\) is the DAE system
!! residual function, and \(P\) is the preconditioner used in the Krylov iteration.

   use daskr_kinds, only: rk, one, two
   use daskr, only: res_t, psol_t
   external :: xerrwd

   integer, intent(in) :: neq
      !! Problem size.
   real(rk), intent(inout) :: y(neq)
      !! On entry, the current solution vector. On return, the new solution vector.
   real(rk), intent(in) :: t
      !! Independent variable.
   real(rk), intent(inout) :: ydot(neq)
      !! On entry, the current derivative vector. On return, the new derivative vector.
   real(rk), intent(out) :: savr(neq)
      !! New residual vector.
   real(rk), intent(in) :: cj
      !! Scalar used in forming the system Jacobian.
   real(rk), intent(in) :: tscale ! @todo: what is "t"?
      !! Scale factor in `t`; used for stopping tests if nonzero.
   real(rk), intent(inout) :: p(neq)
      !! Approximate Newton step used in backtracking.
   real(rk), intent(inout) :: pnorm
      !! Weighted root mean square norm of `p`.
   real(rk), intent(inout) :: wght(neq)
      !! Scaling factors.
   real(rk), intent(in) :: sqrtn
      !! Square root of `neq`.
   real(rk), intent(in) :: rsqrtn
      !! Reciprocal of square root of `neq`.
   integer, intent(in) :: lsoff ! @todo: convert to logical
      !! Flag showing whether the linesearch algorithm is to be invoked.
      !! `0`: do the linesearch.
      !! `1`: turn off linesearch.
   real(rk), intent(in) :: stptol
      !! Tolerance used in calculating the minimum lambda (`rl`) value allowed.
   integer, intent(out) :: iret
      !! Return flag.
      !! `0`: a satisfactory (y,y') was found.
      !! `1`: the routine failed to find a new (y,y') pair that was sufficiently distinct
      !! from the current one.
      !! `2`: a failure in `res` or `psol`.
   procedure(res_t) :: res
      !! User-defined residuals routine.
   integer, intent(out) :: ires
      !! Error flag from `res`.
   procedure(psol_t) :: psol
      !! User-defined preconditioner routine.
   real(rk), intent(inout) :: rwm(*)
      !! Real workspace for the linear system solver.
   integer, intent(inout) :: iwm(*)
      !! Integer workspace for the linear system solver.
   real(rk), intent(in) :: rhok
      !! Weighted norm of preconditioned Krylov residual (not used).
   real(rk), intent(inout) :: fnorm
      !! Value of \( \sqrt{2f} \) for the current \( (y,\dot{y}) \) on input and output.
   integer, intent(in) :: icopt
      !! Initial condition flag.
   integer, intent(in) :: id(neq)
      !! Array indicating which variables are differential and which are algebraic.
   real(rk), intent(inout) :: rwp(*)
      !! Real workspace for the linear system solver.
   integer, intent(inout) :: iwp(*)
      !! Integer workspace for the linear system solver.
   real(rk), intent(out) :: r(neq)
      !! Work space for  the residual vector.
   real(rk), intent(in) :: epslin
      !! Tolerance for linear system.
   real(rk), intent(out) :: ynew(neq)
      !! New solution vector.
   real(rk), intent(out) :: ydotnew(neq)
      !! New derivative of solution vector.
   real(rk), intent(inout) :: pwk(*)
      !! Real work array used by `psol`.
   integer, intent(in) :: icnflg
      !! Constraint flag. If nonzero, then constraint violations in the proposed new
      !! approximate solution will be checked for, and the maximum step length will be
      !! adjusted accordingly.
   integer, intent(out) :: icnstr(neq)
      !! Flags for checking constraints.
   real(rk), intent(in) :: rlx
      !! Factor for restricting update size in [[dcnstr]].
   real(rk), intent(inout) :: rpar(*)
      !! User real workspace.
   integer, intent(inout) :: ipar(*)
      !! User integer workspace.

   integer, parameter :: lnres = 12, lnpsol = 21, lkprint = 31
   real(rk), parameter :: alpha = 1e-4_rk
   integer :: ierr, ivar, kprint
   real(rk) :: f1norm, f1normp, fnormp, ratio, ratio1, rl, rlmin, slpi, tau
   character(len=80) :: msg

   kprint = iwm(lkprint)
   f1norm = (fnorm**2)/2
   ratio = one

   if (kprint >= 2) then
      msg = '------ IN ROUTINE DLINSK-- PNRM = (R1)'
      call xerrwd(msg, 38, 921, 0, 0, 0, 0, 1, pnorm, zero)
   end if
   tau = pnorm
   rl = one

   ! Check for violations of the constraints, if any are imposed.
   ! If any violations are found, the step vector P is rescaled, and the
   ! constraint check is repeated, until no violations are found.
   if (icnflg /= 0) then
      do
         call dyypnw(neq, y, ydot, cj, rl, p, icopt, id, ynew, ydotnew)
         call dcnstr(neq, y, ynew, icnstr, tau, rlx, iret, ivar)
         if (iret /= 1) exit
         ratio1 = tau/pnorm
         ratio = ratio*ratio1
         p = p*ratio1
         pnorm = tau
         if (kprint >= 2) then
            msg = '------ CONSTRAINT VIOL., PNRM = (R1), INDEX = (I1)'
            call xerrwd(msg, 50, 922, 0, 1, ivar, 0, 1, pnorm, zero)
         end if
         if (pnorm <= stptol) then
            iret = 1
            return
         end if
      end do
   end if

   slpi = -two*f1norm*ratio
   rlmin = stptol/pnorm
   if ((lsoff == 0) .and. (kprint >= 2)) then
      msg = '------ MIN. LAMBDA = (R1)'
      call xerrwd(msg, 25, 923, 0, 0, 0, 0, 1, rlmin, zero)
   end if

   ! Begin iteration to find RL value satisfying alpha-condition.
   ! Update YNEW and YDOTNEW, then compute norm of new scaled residual and
   ! perform alpha condition test.
   do
      call dyypnw(neq, y, ydot, cj, rl, p, icopt, id, ynew, ydotnew)
      call dfnrmk(neq, ynew, t, ydotnew, savr, r, cj, tscale, wght, &
                  sqrtn, rsqrtn, res, ires, psol, 0, ierr, fnormp, epslin, &
                  rwp, iwp, pwk, rpar, ipar)
      iwm(lnres) = iwm(lnres) + 1
      if (ires >= 0) iwm(lnpsol) = iwm(lnpsol) + 1
      if (ires /= 0 .or. ierr /= 0) then
         iret = 2
         return
      end if

      if (lsoff == 1) goto 150

      f1normp = (fnormp**2)/2
      if (kprint >= 2) then
         msg = '------ LAMBDA = (R1)'
         call xerrwd(msg, 20, 924, 0, 0, 0, 0, 1, rl, zero)
         msg = '------ NORM(F1) = (R1),  NORM(F1NEW) = (R2)'
         call xerrwd(msg, 43, 925, 0, 0, 0, 0, 2, f1norm, f1normp)
      end if

      if (f1normp > f1norm + alpha*slpi*rl) goto 200

      ! Alpha-condition is satisfied, or linesearch is turned off.
      ! Copy YNEW,YPNEW to Y,YPRIME and return.
150   continue
      iret = 0
      y = ynew
      ydot = ydotnew
      fnorm = fnormp
      if (kprint >= 1) then
         msg = '------ LEAVING ROUTINE DLINSK, FNRM = (R1)'
         call xerrwd(msg, 42, 926, 0, 0, 0, 0, 1, fnorm, zero)
      end if
      return

      ! Alpha-condition not satisfied.  Perform backtrack to compute new RL
      ! value.  If RL is less than RLMIN, i.e. no satisfactory YNEW,YPNEW can
      ! be found sufficiently distinct from Y,YPRIME, then return IRET = 1.
200   continue
      if (rl < rlmin) then
         iret = 1
         return
      end if

      rl = rl/2
   end do

end subroutine dlinsk


subroutine dfnrmk( &
   neq, y, t, ydot, savr, r, cj, tscale, wght, &
   sqrtn, rsqrtn, res, ires, psol, irin, ierr, &
   rnorm, epslin, rwp, iwp, wk, rpar, ipar)
!! This routine calculates the scaled preconditioned norm of the nonlinear function used
!! in the nonlinear iteration for obtaining consistent initial conditions. Specifically,
!! it calculates the weighted root-mean-square norm of the vector:
!!
!! $$  r = P^{-1} G(t,y,\dot{y}) $$
!!
!! where \(P\) is the preconditioner matrix and \(G\) is the DAE equation vector.
   
   use daskr_kinds, only: rk, zero
   use daskr, only: res_t, psol_t
   implicit none

   integer, intent(in) :: neq
      !! Problem size.
   real(rk), intent(in) :: y(neq)
      !! Solution vector.
   real(rk), intent(in) :: t
      !! Independent variable.
   real(rk), intent(in) :: ydot(neq)
      !! Derivative of solution vector after successful step.
   real(rk), intent(inout) :: savr(neq)
      !! Saved residual vector.
   real(rk), intent(out) :: r(neq)
      !! Result vector.
   real(rk), intent(in) :: cj
      !! Scalar used in forming the system Jacobian.
   real(rk), intent(in) :: tscale ! @todo: what is "t?
      !! Scale factor in `t`; used for stopping tests if nonzero. 
   real(rk), intent(inout) :: wght(neq)
      !! Scaling factors.
   real(rk), intent(in) :: sqrtn
      !! Square root of `neq`.
   real(rk), intent(in) :: rsqrtn
      !! Reciprocal of square root of `neq`.   
   procedure(res_t) :: res
      !! User-defined residuals routine.
   integer, intent(out) :: ires
      !! Error flag from `res`.
   procedure(psol_t) :: psol
      !! User-defined preconditioner routine.
   integer, intent(in) :: irin
      !! Flag indicating whether the current residual vector is input in `savr`.
      !! `0`: it is not.
      !! `1`: it is.
   integer, intent(out) :: ierr
      !! Error flag from `psol`.
   real(rk), intent(out) :: rnorm
      !! Weighted root-mean-square norm of `r`.
   real(rk), intent(in) :: epslin
      !! Tolerance for linear system.
   real(rk), intent(inout) :: rwp(*)
      !! Real work array used by preconditioner `psol`.
   integer, intent(inout) :: iwp(*)
      !! Integer work array used by preconditioner `psol`.
   real(rk), intent(inout) :: wk(neq)
      !! Real work array used by `psol`.
   real(rk), intent(inout) :: rpar(*)
      !! User real workspace.
   integer, intent(inout) :: ipar(*)
      !! User integer workspace.

   real(rk) :: ddwnrm ! @todo: remove this once inside module

   ! Call RES routine if IRIN = 0.
   if (irin == 0) then
      ires = 0
      call res(t, y, ydot, cj, savr, ires, rpar, ipar)
      if (ires < 0) return
   end if

   ! Apply inverse of left preconditioner to vector R.
   ! First scale WT array by 1/sqrt(N), and undo scaling afterward.
   call dcopy(neq, savr, 1, r, 1)
   call dscal(neq, rsqrtn, wght, 1)
   ierr = 0
   call psol(neq, t, y, ydot, savr, wk, cj, wght, rwp, iwp, r, epslin, ierr, rpar, ipar)
   call dscal(neq, sqrtn, wght, 1)
   if (ierr /= 0) return

   ! Calculate norm of R.
   rnorm = ddwnrm(neq, r, wght, rpar, ipar)
   if (tscale > zero) rnorm = rnorm*tscale*abs(cj)

end subroutine dfnrmk

subroutine dnedk( &
   t, y, ydot, neq, res, jac, psol, &
   h, wt, jstart, idid, rpar, ipar, phi, gama, savr, delta, e, &
   rwm, iwm, cj, cjold, cjlast, s, uround, epli, sqrtn, rsqrtn, &
   epscon, jcalc, jflag, kp1, nonneg, ntype, iernls)
!! This routine solves a nonlinear system of algebraic equations of the form:
!!
!!  $$ G(t, y, \dot{y}) = 0 $$
!!
!! for the unknown \(y\). The method used is a matrix-free Newton scheme.

   use daskr_kinds, only: rk, zero, one
   use daskr, only: res_t, jack_t, psol_t
   implicit none

   real(rk), intent(in) :: t
      !! Independent variable.
   real(rk), intent(inout) :: y(neq)
      !! Solution vector.
   real(rk), intent(inout) :: ydot(neq)
      !! Derivative of solution vector after successful step.
   integer, intent(in) :: neq
      !! Problem size.
   procedure(res_t) :: res
      !! User-defined residuals routine.
   procedure(jack_t) :: jac
      !! User-defined Jacobian routine.
   procedure(psol_t) :: psol
      !! User-defined preconditioner routine.
   real(rk), intent(in) :: h
      !! Step size.
   real(rk), intent(inout) :: wt(neq) ! @todo: confusing notation: wt ?= ewt ?= whgt
      !! Weights for error control.
   integer, intent(in) :: jstart ! @todo: ? convert to logical
      !! Flag indicating whether this is the first call to this routine.
      !! If `jstart = 0`, then this is the first call, otherwise it is not.
   integer, intent(out) :: idid
      !! Completion flag.
   real(rk), intent(inout) :: rpar(*)
      !! User real workspace.
   integer, intent(inout) :: ipar(*)
      !! User integer workspace.
   real(rk), intent(in) :: phi(neq, kp1)
      !! Array of divided differences used by the Newton iteration.
   real(rk), intent(in) :: gama(kp1)
      !! Array used to predict `y` and `ydot`.
   real(rk), intent(out) :: savr(neq)
      !! Saved residual vector.
   real(rk), intent(inout) :: delta(neq)
      !! Real workspace.
   real(rk), intent(out) :: e(neq)
      !! Error accumulation vector.
   real(rk), intent(inout) :: rwm(*)
      !! Real workspace for the linear system solver.
   integer, intent(inout) :: iwm(*)
      !! Integer workspace for the linear system solver.
   real(rk), intent(in) :: cj
      !! Scalar used in forming the system Jacobian.
   real(rk), intent(out) :: cjold
      !! Value of `cj` as of the last call to [[ditmd]]. Accounts for changes in `cj`
      !! needed to decide whether to call [[ditmd]].
   real(rk), intent(in) :: cjlast
      !! Previous value of `cj`.
   real(rk), intent(out) :: s
      !! Convergence factor for the Newton iteration. See [[dnsk]] for details. On the
      !! first Newton iteration with an updated preconditioner, `s = 100`. The value of
      !! `s` is preserved from call to call so that the rate estimate from a previous
      !! step can be applied to the current step.
   real(rk), intent(in) :: uround ! @todo: ? remove
      !! Unit roundoff (not used).
   real(rk), intent(in) :: epli ! @note: not the same as 'epslin'!
      !! Convergence test constant for the iterative solution of linear systems. See
      !! [[daskr]] for details.
   real(rk), intent(in) :: sqrtn
      !! Square root of `neq`.
   real(rk), intent(in) :: rsqrtn
      !! Reciprocal of square root of `neq`.
   real(rk), intent(in) :: epscon
      !! Tolerance for convergence of the Newton iteration.
   integer, intent(out) :: jcalc
      !! Flag indicating whether the Jacobian matrix needs to be updated.
      !! `-1`: call [[ditmd]] to update the Jacobian matrix.
      !!  `0`: Jacobian matrix is up-to-date.
      !!  `1`: Jacobian matrix is outdated, but [[ditmd]] will not be called unless `jcalc` is set to `-1`.
   integer, intent(in) :: jflag
      !! Flag indicating whether a `jac` routine is supplied by the user.
   integer, intent(in) :: kp1
      !! The current order + 1. Updated across calls.
   integer, intent(in) :: nonneg
      !! Flag indicating whether nonnegativity constraints are imposed on the solution.
      !! `0`: no nonnegativity constraints.
      !! `1`: nonnegativity constraints are imposed.
   integer, intent(in) :: ntype
      !! Type of the nonlinear solver.
      !! `1`: modified Newton with iterative linear solver.
      !! `2`: modified Newton with user-supplied linear solver.
   integer, intent(out) :: iernls
      !! Error flag for the nonlinear solver.
      !!  `0`: nonlinear solver converged.
      !!  `1`: recoverable error inside non-linear solver.
      !! `-1`: unrecoverable error inside non-linear solver.

   integer, parameter :: lnre = 12, lnje = 13, llocwp = 29, llciwp = 30, &
                         muldel = 0, maxit = 4
   real(rk), parameter :: xrate = 0.25_rk

   integer :: iernew, ierpj, iersl, iertyp, ires, j, liwp, lwp
   real(rk) :: delnrm, epslin, temp1, temp2, tolnew
   real(rk) :: ddwnrm ! @todo: remove this once inside module

   ! Verify that this is the correct subroutine.
   iertyp = 0
   if (ntype /= 1) then
      iertyp = 1
      goto 380
   end if

   ! If this is the first step, perform initializations.
   if (jstart == 0) then
      cjold = cj
      jcalc = -1
      s = 1e2_rk
   end if

   ! Perform all other initializations.
   iernls = 0
   lwp = iwm(llocwp)
   liwp = iwm(llciwp)

   ! Decide whether to update the preconditioner.
   if (jflag /= 0) then
      temp1 = (one - xrate)/(one + xrate)
      temp2 = one/temp1
      if ((cj/cjold < temp1) .or. (cj/cjold > temp2)) jcalc = -1
      if (cj /= cjlast) s = 1e2_rk
   else
      jcalc = 0
   end if

   ! Looping point for updating preconditioner with current stepsize.
300 continue

   ! Initialize all error flags to zero.
   ierpj = 0
   ires = 0
   iersl = 0
   iernew = 0

   ! Predict the solution and derivative and compute the tolerance for the Newton iteration.
   y = phi(:, 1)
   ydot = zero
   do j = 2, kp1
      y = y + phi(:, j)
      ydot = ydot + gama(j)*phi(:, j)
   end do

330 continue
   epslin = epli*epscon
   tolnew = epslin

   ! Call RES to initialize DELTA.
   iwm(lnre) = iwm(lnre) + 1
   call res(t, y, ydot, cj, delta, ires, rpar, ipar)
   if (ires < 0) goto 380

   ! If indicated, update the preconditioner.
   ! Set JCALC to 0 as an indicator that this has been done.
   if (jcalc == -1) then
      iwm(lnje) = iwm(lnje) + 1
      jcalc = 0
      call jac(res, ires, neq, t, y, ydot, wt, delta, e, h, cj, rwm(lwp), iwm(liwp), ierpj, rpar, ipar)
      cjold = cj
      s = 1e2_rk
      if (ires < 0) goto 380
      if (ierpj /= 0) goto 380
   end if

   ! Call the nonlinear Newton solver.

   call dnsk(t, y, ydot, neq, res, psol, wt, rpar, ipar, savr, &
             delta, e, rwm, iwm, cj, sqrtn, rsqrtn, epslin, epscon, &
             s, temp1, tolnew, muldel, maxit, ires, iersl, iernew)

   if (iernew > 0 .and. jcalc /= 0) then
      ! The Newton iteration had a recoverable failure with an old
      ! preconditioner. Retry the step with a new preconditioner.
      jcalc = -1
      goto 300
   end if

   if (iernew /= 0) goto 380

   ! The Newton iteration has converged. If nonnegativity of solution is required, set
   ! the solution nonnegative, if the perturbation to do it is small enough. If the
   ! change is too large, then consider the corrector iteration to have failed.
   if (nonneg == 0) goto 390
   delta = min(y, zero)
   delnrm = ddwnrm(neq, delta, wt, rpar, ipar)
   if (delnrm > epscon) goto 380
   e = e - delta
   goto 390

   ! Exits from nonlinear solver.
   ! No convergence with current preconditioner.
   ! Compute IERNLS and IDID accordingly.
380 continue
   if ((ires <= -2) .or. (iersl < 0) .or. (iertyp /= 0)) then
      iernls = -1
      if (ires <= -2) idid = -11
      if (iersl < 0) idid = -13
      if (iertyp /= 0) idid = -15
   else
      iernls = 1
      if (ires == -1) idid = -10
      if (ierpj /= 0) idid = -5
      if (iersl > 0) idid = -14
   end if

390 continue
   jcalc = 1

end subroutine dnedk

subroutine dnsk( &
   t, y, ydot, neq, res, psol, wt, rpar, ipar, &
   savr, delta, e, rwm, iwm, cj, sqrtn, rsqrtn, epslin, epscon, &
   s, confac, tolnew, muldel, maxit, ires, iersl, iernew)
!! This routines solves a nonlinear system of algebraic equations of the form:
!!
!!  $$ G(t, y, \dot{y}) = 0 $$
!!
!! for the unknown \(y\). The method used is a modified Newton scheme.

   use daskr_kinds, only: rk, zero, one
   use daskr, only: res_t, psol_t
   implicit none

   real(rk), intent(in) :: t
      !! Independent variable.
   real(rk), intent(inout) :: y(neq)
      !! Solution vector.
   real(rk), intent(inout) :: ydot(neq)
      !! Derivative of solution vector.
   integer, intent(in) :: neq
      !! Problem size.
   procedure(res_t) :: res
      !! User-defined residuals routine.
   procedure(psol_t) :: psol
      !! User-defined preconditioner routine.
   real(rk), intent(inout) :: wt(neq) ! @todo: confusing notation: wt ?= ewt ?= whgt
      !! Weights for error control.
   real(rk), intent(inout) :: rpar(*)
      !! User real workspace.
   integer, intent(inout) :: ipar(*)
      !! User integer workspace.
   real(rk), intent(out) :: savr(neq)
      !! Saved residual vector.
   real(rk), intent(inout) :: delta(neq)
      !! Real workspace.
   real(rk), intent(out) :: e(neq)
      !! Error accumulation vector.
   real(rk), intent(inout) :: rwm(*)
      !! Real workspace for the linear system solver.
   integer, intent(inout) :: iwm(*)
      !! Integer workspace for the linear system solver.
   real(rk), intent(in) :: cj
      !! Scalar used in forming the system Jacobian.
   real(rk), intent(in) :: sqrtn
      !! Square root of `neq`.
   real(rk), intent(in) :: rsqrtn
      !! Reciprocal of square root of `neq`.
   real(rk), intent(in) :: epslin
      !! Tolerance for linear system solver.
   real(rk), intent(in) :: epscon
      !! Tolerance for convergence of the Newton iteration.
   real(rk), intent(out) :: s
      !! Convergence factor for the Newton iteration. `s=rate/(1 - rate)`, where
      !! `rate` is the estimated rate of convergence of the Newton iteration.
      !! The closer `rate` is to 0, the faster the Newton iteration is converging.
      !! the closer `rate` is to 1, the slower the Newton iteration is converging.
   real(rk), intent(in) :: confac
      !! Residual scale factor to improve convergence.
   real(rk), intent(in) :: tolnew
      !! Tolerance on the norm of the Newton correction in the alternative Newton
      !! convergence test.
   integer, intent(in) :: muldel ! @todo: convert to logical
      !! Flag indicating whether or not to multiply `delta` by `confac`.
   integer, intent(in) :: maxit
      !! Maximum number of Newton iterations allowed.
   integer, intent(out) :: ires
      !! Error flag from `res`. If `ires < -1`, then `iernew = -1`.
   integer, intent(out) :: iersl
      !! Error flag from the linear system solver. See [[dslvk]] for details.
      !! If `iersl < 0`, then `iernew = -1`.
      !! If `iersl = 1`, then `iernew = 1`.
   integer, intent(out) :: iernew
      !! Error flag for the Newton iteration.
      !!  `0`: the Newton iteration converged.
      !!  `1`: recoverable error inside Newton iteration.
      !! `-1`: unrecoverable error inside Newton iteration.

   integer, parameter :: lnni = 19, lnre = 12

   integer :: m
   real(rk) :: delnrm, oldnrm, rate, rhok
   real(rk) :: ddwnrm ! @todo: remove this once inside module
   logical :: converged

   ! Initialize Newton counter M and accumulation vector E.
   m = 0
   e = zero

   ! Corrector loop.
   converged = .false.
   do
      iwm(lnni) = iwm(lnni) + 1

      ! If necessary, multiply residual by convergence factor.
      if (muldel == 1) then
         delta = delta*confac
      end if

      ! Save residual in SAVR.
      savr = delta

      ! Compute a new iterate. Store the correction in DELTA.
      call dslvk(neq, y, t, ydot, savr, delta, wt, rwm, iwm, &
                 res, ires, psol, iersl, cj, epslin, sqrtn, rsqrtn, rhok, &
                 rpar, ipar)
      if ((ires /= 0) .or. (iersl /= 0)) exit

      ! Update Y, E, and YDOT.
      y = y - delta
      e = e - delta
      ydot = ydot - cj*delta

      ! Test for convergence of the iteration.
      delnrm = ddwnrm(neq, delta, wt, rpar, ipar)
      if (m == 0) then
         oldnrm = delnrm
         if (delnrm <= tolnew) then
            converged = .true.
            exit
         end if
      else
         rate = (delnrm/oldnrm)**(one/m)
         if (rate > 0.9_rk) exit
         s = rate/(one - rate)
      end if

      if (s*delnrm <= epscon) then
         converged = .true.
         exit
      end if
      
      ! The corrector has not yet converged. Update M and test whether
      ! the maximum number of iterations have been tried.
      m = m + 1
      if (m >= maxit) exit

      ! Evaluate the residual, and go back to do another iteration.
      iwm(lnre) = iwm(lnre) + 1
      call res(t, y, ydot, cj, delta, ires, rpar, ipar)
      if (ires < 0) exit

   end do

   if (.not. converged) then
      if ((ires <= -2) .or. (iersl < 0)) then
         iernew = -1
      else
         iernew = 1
      end if
   end if

end subroutine dnsk

subroutine dslvk( &
   neq, y, t, ydot, savr, x, ewt, rwm, iwm, res, ires, psol, &
   iersl, cj, epslin, sqrtn, rsqrtn, rhok, rpar, ipar)
!! This routine uses a restart algorithm and interfaces to [[dspigm]] for the solution of
!! the linear system arising from a Newton iteration.

   use daskr_kinds, only: rk, zero
   use daskr, only: res_t, psol_t
   use blas_interfaces, only: dscal, dcopy
   implicit none

   integer, intent(in) :: neq
      !! Problem size.
   real(rk), intent(in) :: y(neq)
      !! Current dependent variables.
   real(rk), intent(in) :: t
      !! Current independent variable.
   real(rk), intent(in) :: ydot(neq)
      !! Current derivatives of dependent variables.
   real(rk), intent(in) :: savr(neq)
      !! Current residual evaluated at `(t, y, ydot)`.
   real(rk), intent(inout) :: x(neq)
      !! On entry, the right-hand side vector of the linear system to be solved,
      !! and, on exit, the solution vector.
   real(rk), intent(inout) :: ewt(neq)
      !! Nonzero elements of the diagonal scaling matrix.
   real(rk), intent(inout) :: rwm(*)
      !! Real work space containing data for the algorithm (Krylov basis vectors,
      !! Hessenberg matrix, etc.).
   integer, intent(inout) :: iwm(*)
      !! Integer work space containing data for the algorithm.
   procedure(res_t) :: res
      !! User-defined residuals routine.
   integer, intent(out) :: ires
      !! Error flag from `res`.
   procedure(psol_t) :: psol
      !! User-defined preconditioner routine.
   integer, intent(out) :: iersl
      !! Error flag.
      !!  `0`: no trouble occurred (or user `res` routine returned `ires < 0`).
      !!  `1`: iterative method failed to converge ([[dspigm]] returned `iflag > 0`).
      !! `-1`: nonrecoverable error in the iterative solver, and an error exit will occur.
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
   do while ((iflag == 1) .and. (nrsts < nrmax) .and. (ires == 0))
      nrsts = nrsts + 1
      if (nrsts > 0) call dcopy(neq, rwm(ldl), 1, rwm(lr), 1)
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
   if (ires < 0) then
      ncfl = ncfl + 1
   elseif (iflag /= 0) then
      ncfl = ncfl + 1
      if (iflag > 0) iersl = 1
      if (iflag < 0) iersl = -1
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
   use daskr, only: res_t, psol_t
   use blas_interfaces, only: dnrm2, dcopy, dscal, daxpy
   implicit none

   integer, intent(in) :: neq
      !! Problem size.
   real(rk), intent(in) :: t
      !! Current independent variable.
   real(rk), intent(in) :: y(neq)
      !! Current dependent variables.
   real(rk), intent(in) :: ydot(neq)
      !! Current derivatives of dependent variables.
   real(rk), intent(in) :: savr(neq)
      !! Current residual evaluated at `(t, y, ydot)`.
   real(rk), intent(inout) :: r(neq)
      !! On entry, the right hand side vector. Also used as workspace when computing
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
   procedure(res_t) :: res
      !! User-defined residuals routine.
   integer, intent(out) :: ires
      !! Error flag from `res`.
   integer, intent(inout) :: nres
      !! Number of calls to `res`.
   procedure(psol_t) :: psol
      !! User-defined preconditioner routine.
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
      !! `0`: convergence in `lgmr` iterations, `lgmr <= maxl`.
      !! `1`: convergence test did not pass in `maxl` iterations, but the new
      !! residual norm (`rho`) is less than the old residual norm (`rnrm`), and so `z`
      !! is computed.
      !! `2`: convergence test did not pass in `maxl` iterations, new residual
      !! norm (`rho`) is greater than or equal to old residual norm (`rnrm`), and the
      !! initial guess, `z = 0`, is returned.
      !! `3`: recoverable error in `psol` caused by the preconditioner being out of date.
      !! `-1`: unrecoverable error in `psol`.
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

   ! Apply inverse of left preconditioner to vector R if NRSTS == 0.
   ! Form V(:,1), the scaled preconditioned right hand side.
   if (nrsts == 0) then
      call psol(neq, t, y, ydot, savr, wk, cj, wght, rwp, iwp, r, epslin, ierr, rpar, ipar)
      npsol = 1
      if (ierr /= 0) goto 300
      v(:, 1) = r*wght
   else
      v(:, 1) = r
   end if

   ! Calculate norm of scaled vector V(:,1) and normalize it
   ! If, however, the norm of V(:,1) (i.e. the norm of the preconditioned
   ! residual) is <= EPLIN, then return with Z=0.
   rnrm = dnrm2(neq, v, 1)
   if (rnrm <= epslin) then
      rhok = rnrm
      return
   end if
   tem = one/rnrm
   call dscal(neq, tem, v(1, 1), 1)

   ! Zero out the HES array.
   hes = zero

   ! Main loop to compute the vectors V(:,2) to V(:,MAXL).
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
      if (ires < 0) return
      if (ierr /= 0) goto 300

      ! Call routine DORTH to orthogonalize the new vector VNEW = V(:,LL+1).
      call dorth(v(1, ll + 1), v, hes, neq, ll, maxlp1, kmp, snormw)
      hes(ll + 1, ll) = snormw

      ! Call routine DHEQR to update the factors of HES.
      call dheqr(hes, maxlp1, ll, q, info, ll)
      if (info == ll) goto 120

      ! Update RHO, the estimate of the norm of the residual R - A*ZL.
      ! If KMP < MAXL, then the vectors V(:,1),...,V(:,LL+1) are not
      ! necessarily orthogonal for LL > KMP.  The vector DL must then
      ! be computed, and its norm used in the calculation of RHO.
      prod = prod*q(2*ll)
      rho = abs(prod*rnrm)
      if ((ll > kmp) .and. (kmp < maxl)) then
         if (ll == kmp + 1) then
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
      ! If failed and LL < MAXL, then continue iterating.
      if (rho <= epslin) goto 200
      if (ll == maxl) goto 100

      ! Rescale so that the norm of V(1,LL+1) is one.
      tem = one/snormw
      call dscal(neq, tem, v(1, ll + 1), 1)

   end do

100 continue
   if (rho < rnrm) goto 150

120 continue
   iflag = 2
   z = zero
   return

150 continue
   ! The tolerance was not met, but the residual norm was reduced.
   ! If performing restarting (IRST > 0) calculate the residual vector
   ! RL and store it in the DL array.  If the incomplete version is
   ! being used (KMP < MAXL) then DL has already been calculated.
   iflag = 1
   if (irst > 0) then

      if (kmp == maxl) then
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
   if (ierr < 0) iflag = -1
   if (ierr > 0) iflag = 3

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
   use daskr, only: res_t, psol_t

   integer, intent(in) :: neq
      !! Problem size.
   real(rk), intent(in) :: y(neq)
      !! Current dependent variables.
   real(rk), intent(in) :: t
      !! Current independent variable.
   real(rk), intent(in) :: ydot(neq)
      !! Current derivatives of dependent variables.
   real(rk), intent(in) :: savr(neq)
      !! Current residual evaluated at `(t, y, ydot)`.
   real(rk), intent(in) :: v(neq)
      !! Orthonormal vector (can be the same as `z`).  
   real(rk), intent(in) :: wght(neq)
      !! Scaling factors. `1/wght(i)` are the diagonal elements of the matrix \(D\).
   real(rk), intent(out) :: ydottemp(neq)
      !! Work array used to store the incremented value of `ydot`.
   procedure(res_t) :: res
      !! User-defined residuals routine.
   integer, intent(out) :: ires
      !! Error flag from `res`.
   procedure(psol_t) :: psol
      !! User-defined preconditioner routine.
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
   if (ires < 0) return

   ! Set Z = (dF/dY) * VBAR using difference quotient.
   ! (VBAR is old value of VTEMP before calling RES)
   z = vtemp - savr

   ! Apply inverse of left preconditioner to Z.
   call psol(neq, t, y, ydot, savr, ydottemp, cj, wght, rwp, iwp, z, epslin, ierr, rpar, ipar)
   npsol = npsol + 1
   if (ierr /= 0) return

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
      !! On entry, vector containing a scaled product of the Jacobian and the vector `v(:,ll)`.
      !! On return, the new vector orthogonal to `v(:,i0)`, where `i0 = max(1, ll - kmp + 1)`.
   real(rk), intent(in) :: v(n, ll)
      !! Matrix containing the previous `ll` orthogonal vectors `v(:,1)` to `v(:,ll)`.
   real(rk), intent(inout) :: hes(ldhes, ll)
      !! On entry, an upper Hessenberg matrix of shape `(ll, ll)` containing in `hes(i,k)`,
      !! for `k < ll`, the scaled inner products of `a*v(:,k)` and `v(:,i)`.
      !! On return, an upper Hessenberg matrix with column `ll` filled in with the scaled
      !! inner products of `a*v(:,ll)` and `v(:,i)`.
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
   if (vnrm + snormw/1000 /= vnrm) return ! @todo: fix this comparison

   sumdsq = zero
   do i = i0, ll
      tem = -ddot(n, v(1, i), 1, vnew, 1)
      if (hes(i, ll) + tem/1000 == hes(i, ll)) cycle ! @todo: fix this comparison
      hes(i, ll) = hes(i, ll) - tem
      call daxpy(n, tem, v(1, i), 1, vnew, 1)
      sumdsq = sumdsq + tem**2
   end do
   if (sumdsq == zero) return

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
         if (km1 < 1) goto 20
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
         if (t2 /= zero) goto 30
         c = one
         s = zero
         goto 50

30       continue
         if (abs(t2) < abs(t1)) goto 40
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
         if (a(k, k) == zero) info = k
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
      if (t2 /= zero) goto 110
      c = one
      s = zero
      goto 130

110   continue
      if (abs(t2) < abs(t1)) goto 120
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
      if (a(n, n) == zero) info = n

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
