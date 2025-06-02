module linpack
!! Selected routines from LINPACK.

   use daskr_kinds, only: rk, zero, one
   implicit none

contains

   pure subroutine dgbfa(abd, lda, n, ml, mu, ipvt, info)
   !! DGBFA factors a real band matrix by elimination.
   !
   !  Discussion:
   !
   !    DGBFA is usually called by DGBCO, but it can be called
   !    directly with a saving in time if RCOND is not needed.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    17 May 2005
   !
   !  Author:
   !
   !    FORTRAN90 version by John Burkardt.
   !
   !  Reference:
   !
   !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
   !    LINPACK User's Guide,
   !    SIAM, 1979,
   !    ISBN13: 978-0-898711-72-1,
   !    LC: QA214.L56.
   !
   !  Parameters:
   !
   !    Input/output, real ( kind = 8 ) ABD(LDA,N).  On input, the matrix in band
   !    storage.  The columns of the matrix are stored in the columns of ABD
   !    and the diagonals of the matrix are stored in rows ML+1 through
   !    2*ML+MU+1 of ABD.  On output, an upper triangular matrix in band storage
   !    and the multipliers which were used to obtain it.  The factorization
   !    can be written A = L*U where L is a product of permutation and unit lower
   !    triangular matrices and U is upper triangular.
   !
   !    Input, integer ( kind = 4 ) LDA, the leading dimension of the array ABD.
   !    2*ML + MU + 1 <= LDA is required.
   !
   !    Input, integer ( kind = 4 ) N, the order of the matrix.
   !
   !    Input, integer ( kind = 4 ) ML, MU, the number of diagonals below and above
   !    the main diagonal.  0 <= ML < N, 0 <= MU < N.
   !
   !    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
   !
   !    Output, integer ( kind = 4 ) INFO, error flag.
   !    0, normal value.
   !    K, if U(K,K) == 0.0.  This is not an error condition for this
   !      subroutine, but it does indicate that DGBSL will divide by zero if
   !      called.  Use RCOND in DGBCO for a reliable indication of singularity.

      use blas_interfaces, only: idamax, daxpy, dscal

      real(rk), intent(inout) :: abd(lda, n)
      integer, intent(in) :: lda
      integer, intent(in) :: n
      integer, intent(in) :: ml
      integer, intent(in) :: mu
      integer, intent(out) :: ipvt(n)
      integer, intent(out) :: info

      integer :: i
      integer :: i0
      integer :: j
      integer :: j0
      integer :: j1
      integer :: ju
      integer :: jz
      integer :: k
      integer :: l
      integer :: lm
      integer :: m
      integer :: mm
      real(rk) :: t

      m = ml + mu + 1
      info = 0
   !
   !  Zero initial fill-in columns.
   !
      j0 = mu + 2
      j1 = min(n, m) - 1

      do jz = j0, j1
         i0 = m + 1 - jz
         do i = i0, ml
            abd(i, jz) = zero
         end do
      end do

      jz = j1
      ju = 0
   !
   !  Gaussian elimination with partial pivoting.
   !
      do k = 1, n - 1
   !
   !  Zero out the next fill-in column.
   !
         jz = jz + 1
         if (jz <= n) then
            abd(1:ml, jz) = zero
         end if
   !
   !  Find L = pivot index.
   !
         lm = min(ml, n - k)
         l = idamax(lm + 1, abd(m, k), 1) + m - 1
         ipvt(k) = l + k - m
   !
   !  Zero pivot implies this column already triangularized.
   !
         if (abd(l, k) == zero) then

            info = k
   !
   !  Interchange if necessary.
   !
         else

            if (l /= m) then
               t = abd(l, k)
               abd(l, k) = abd(m, k)
               abd(m, k) = t
            end if
   !
   !  Compute multipliers.
   !
            t = -one/abd(m, k)
            call dscal(lm, t, abd(m + 1, k), 1)
   !
   !  Row elimination with column indexing.
   !
            ju = min(max(ju, mu + ipvt(k)), n)
            mm = m

            do j = k + 1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l, j)
               if (l /= mm) then
                  abd(l, j) = abd(mm, j)
                  abd(mm, j) = t
               end if
               call daxpy(lm, t, abd(m + 1, k), 1, abd(mm + 1, j), 1)
            end do

         end if

      end do

      ipvt(n) = n

      if (abd(m, n) == zero) then
         info = n
      end if

   end subroutine dgbfa

   pure subroutine dgbsl(abd, lda, n, ml, mu, ipvt, b, job)
   !! DGBSL solves a real banded system factored by DGBCO or DGBFA.
   !
   !  Discussion:
   !
   !    DGBSL can solve either A * X = B  or  A' * X = B.
   !
   !    A division by zero will occur if the input factor contains a
   !    zero on the diagonal.  Technically this indicates singularity
   !    but it is often caused by improper arguments or improper
   !    setting of LDA.  It will not occur if the subroutines are
   !    called correctly and if DGBCO has set 0.0 < RCOND
   !    or DGBFA has set INFO == 0.
   !
   !    To compute inverse(A) * C  where C is a matrix with P columns:
   !
   !      call dgbco ( abd, lda, n, ml, mu, ipvt, rcond, z )
   !
   !      if ( rcond is too small ) then
   !        exit
   !      end if
   !
   !      do j = 1, p
   !        call dgbsl ( abd, lda, n, ml, mu, ipvt, c(1,j), 0 )
   !      end do
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    17 May 2005
   !
   !  Author:
   !
   !    FORTRAN90 version by John Burkardt.
   !
   !  Reference:
   !
   !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
   !    LINPACK User's Guide,
   !    SIAM, 1979,
   !    ISBN13: 978-0-898711-72-1,
   !    LC: QA214.L56.
   !
   !  Parameters:
   !
   !    Input, real ( kind = 8 ) ABD(LDA,N), the output from DGBCO or DGBFA.
   !
   !    Input, integer ( kind = 4 ) LDA, the leading dimension of the array ABD.
   !
   !    Input, integer ( kind = 4 ) N, the order of the matrix.
   !
   !    Input, integer ( kind = 4 ) ML, MU, the number of diagonals below and above
   !    the main diagonal.  0 <= ML < N, 0 <= MU < N.
   !
   !    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from DGBCO or DGBFA.
   !
   !    Input/output, real ( kind = 8 ) B(N).  On input, the right hand side.
   !    On output, the solution.
   !
   !    Input, integer ( kind = 4 ) JOB, job choice.
   !    0, solve A*X=B.
   !    nonzero, solve A'*X=B.

      use blas_interfaces, only: daxpy, ddot

      real(rk), intent(in) :: abd(lda, n)
      integer, intent(in) :: lda
      integer, intent(in) :: n
      integer, intent(in) :: ml
      integer, intent(in) :: mu
      integer, intent(in) :: ipvt(n)
      real(rk), intent(inout) :: b(n)
      integer, intent(in) :: job

      integer :: k
      integer :: l
      integer :: la
      integer :: lb
      integer :: lm
      integer :: m
      real(rk) :: t

      m = mu + ml + 1
   !
   !  JOB = 0, Solve A * x = b.
   !
   !  First solve L * y = b.
   !
      if (job == 0) then

         if (0 < ml) then

            do k = 1, n - 1
               lm = min(ml, n - k)
               l = ipvt(k)
               t = b(l)
               if (l /= k) then
                  b(l) = b(k)
                  b(k) = t
               end if
               call daxpy(lm, t, abd(m + 1, k), 1, b(k + 1), 1)
            end do

         end if
   !
   !  Now solve U * x = y.
   !
         do k = n, 1, -1
            b(k) = b(k)/abd(m, k)
            lm = min(k, m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call daxpy(lm, t, abd(la, k), 1, b(lb), 1)
         end do
   !
   !  JOB nonzero, solve A' * x = b.
   !
   !  First solve U' * y = b.
   !
      else

         do k = 1, n
            lm = min(k, m) - 1
            la = m - lm
            lb = k - lm
            t = ddot(lm, abd(la, k), 1, b(lb), 1)
            b(k) = (b(k) - t)/abd(m, k)
         end do
   !
   !  Now solve L' * x = y.
   !
         if (0 < ml) then

            do k = n - 1, 1, -1
               lm = min(ml, n - k)
               b(k) = b(k) + ddot(lm, abd(m + 1, k), 1, b(k + 1), 1)
               l = ipvt(k)
               if (l /= k) then
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
               end if
            end do

         end if

      end if

   end subroutine dgbsl

   pure subroutine dgefa(a, lda, n, ipvt, info)
   !! DGEFA factors a real general matrix.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    07 March 2001
   !
   !  Author:
   !
   !    FORTRAN90 version by John Burkardt.
   !
   !  Reference:
   !
   !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
   !    LINPACK User's Guide,
   !    SIAM, 1979,
   !    ISBN13: 978-0-898711-72-1,
   !    LC: QA214.L56.
   !
   !  Parameters:
   !
   !    Input/output, real ( kind = 8 ) A(LDA,N).
   !    On intput, the matrix to be factored.
   !    On output, an upper triangular matrix and the multipliers used to obtain
   !    it.  The factorization can be written A=L*U, where L is a product of
   !    permutation and unit lower triangular matrices, and U is upper triangular.
   !
   !    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
   !
   !    Input, integer ( kind = 4 ) N, the order of the matrix A.
   !
   !    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
   !
   !    Output, integer ( kind = 4 ) INFO, singularity indicator.
   !    0, normal value.
   !    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
   !    but it does indicate that DGESL or DGEDI will divide by zero if called.
   !    Use RCOND in DGECO for a reliable indication of singularity.

      use blas_interfaces, only: idamax

      real(rk), intent(inout) :: a(lda, n)
      integer, intent(in) :: lda
      integer, intent(in) :: n
      integer, intent(out) :: ipvt(n)
      integer, intent(out) :: info

      integer :: j
      integer :: k
      integer :: l
      real(rk) :: t
   !
   !  Gaussian elimination with partial pivoting.
   !
      info = 0

      do k = 1, n - 1
   !
   !  Find L = pivot index.
   !
         l = idamax(n - k + 1, a(k, k), 1) + k - 1
         ipvt(k) = l
   !
   !  Zero pivot implies this column already triangularized.
   !
         if (a(l, k) == zero) then
            info = k
            cycle
         end if
   !
   !  Interchange if necessary.
   !
         if (l /= k) then
            t = a(l, k)
            a(l, k) = a(k, k)
            a(k, k) = t
         end if
   !
   !  Compute multipliers.
   !
         t = -one/a(k, k)
         a(k + 1:n, k) = a(k + 1:n, k)*t
   !
   !  Row elimination with column indexing.
   !
         do j = k + 1, n
            t = a(l, j)
            if (l /= k) then
               a(l, j) = a(k, j)
               a(k, j) = t
            end if
            a(k + 1:n, j) = a(k + 1:n, j) + t*a(k + 1:n, k)
         end do

      end do

      ipvt(n) = n

      if (a(n, n) == zero) then
         info = n
      end if

   end subroutine dgefa

   pure subroutine dgesl(a, lda, n, ipvt, b, job)
   !! DGESL solves a real general linear system A * X = B.
   !
   !  Discussion:
   !
   !    DGESL can solve either of the systems A * X = B or A' * X = B.
   !
   !    The system matrix must have been factored by DGECO or DGEFA.
   !
   !    A division by zero will occur if the input factor contains a
   !    zero on the diagonal.  Technically this indicates singularity
   !    but it is often caused by improper arguments or improper
   !    setting of LDA.  It will not occur if the subroutines are
   !    called correctly and if DGECO has set 0.0 < RCOND
   !    or DGEFA has set INFO == 0.
   !
   !  Licensing:
   !
   !    This code is distributed under the GNU LGPL license.
   !
   !  Modified:
   !
   !    07 March 2001
   !
   !  Author:
   !
   !    FORTRAN90 version by John Burkardt.
   !
   !  Reference:
   !
   !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
   !    LINPACK User's Guide,
   !    SIAM, 1979,
   !    ISBN13: 978-0-898711-72-1,
   !    LC: QA214.L56.
   !
   !  Parameters:
   !
   !    Input, real ( kind = 8 ) A(LDA,N), the output from DGECO or DGEFA.
   !
   !    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
   !
   !    Input, integer ( kind = 4 ) N, the order of the matrix A.
   !
   !    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from DGECO or DGEFA.
   !
   !    Input/output, real ( kind = 8 ) B(N).
   !    On input, the right hand side vector.
   !    On output, the solution vector.
   !
   !    Input, integer ( kind = 4 ) JOB.
   !    0, solve A * X = B;
   !    nonzero, solve A' * X = B.

      real(rk), intent(in) :: a(lda, n)
      integer, intent(in) :: lda
      integer, intent(in) :: n
      integer, intent(in) :: ipvt(n)
      real(rk), intent(inout) :: b(n)
      integer, intent(in) :: job

      integer :: k
      integer :: l
      real(rk) :: t
   !
   !  Solve A * X = B.
   !
      if (job == 0) then

         do k = 1, n - 1

            l = ipvt(k)
            t = b(l)

            if (l /= k) then
               b(l) = b(k)
               b(k) = t
            end if

            b(k + 1:n) = b(k + 1:n) + t*a(k + 1:n, k)

         end do

         do k = n, 1, -1
            b(k) = b(k)/a(k, k)
            t = -b(k)
            b(1:k - 1) = b(1:k - 1) + t*a(1:k - 1, k)
         end do

      else
   !
   !  Solve A' * X = B.
   !
         do k = 1, n
            t = dot_product(a(1:k - 1, k), b(1:k - 1))
            b(k) = (b(k) - t)/a(k, k)
         end do

         do k = n - 1, 1, -1

            b(k) = b(k) + dot_product(a(k + 1:n, k), b(k + 1:n))
            l = ipvt(k)

            if (l /= k) then
               t = b(l)
               b(l) = b(k)
               b(k) = t
            end if

         end do

      end if

   end subroutine dgesl

end module linpack
