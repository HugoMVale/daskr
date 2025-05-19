module blas_interfaces
  !! Single and double precision interfaces for the BLAS procedures used by odrpack.
   use daskr_kinds, only: sp, dp
   implicit none

   interface

      ! pure subroutine saxpy(n, a, x, incx, y, incy)
      ! !! Computation `Y = A*X + Y`.
      !    import :: sp
      !    integer, intent(in) :: n
      !    real(sp), intent(in) :: a
      !    real(sp), intent(in) :: x(*)
      !    integer, intent(in) :: incx
      !    real(sp), intent(inout) :: y(*)
      !    integer, intent(in) :: incy
      ! end subroutine saxpy

      pure subroutine daxpy(n, a, x, incx, y, incy)
      !! Computation `Y = A*X + Y`.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: a
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine daxpy

      ! pure subroutine scopy(n, x, incx, y, incy)
      ! !! Vector copy `Y = X`.
      !    import :: sp
      !    integer, intent(in) :: n
      !    real(sp), intent(in) :: x(*)
      !    integer, intent(in) :: incx
      !    real(sp), intent(inout) :: y(*)
      !    integer, intent(in) :: incy
      ! end subroutine scopy

      pure subroutine dcopy(n, x, incx, y, incy)
      !! Vector copy `Y = X`.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine dcopy

      ! pure real(sp) function sdot(n, x, incx, y, incy)
      ! !! Inner product of vectors.
      !    import :: sp
      !    integer, intent(in) :: n
      !    real(sp), intent(in) :: x(*)
      !    integer, intent(in) :: incx
      !    real(sp), intent(in) :: y(*)
      !    integer, intent(in) :: incy
      ! end function sdot

      pure real(dp) function ddot(n, x, incx, y, incy)
      !! Inner product of vectors.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(in) :: y(*)
         integer, intent(in) :: incy
      end function ddot

      ! pure real(sp) function snrm2(n, x, incx)
      ! !! Euclidean length (L2 Norm) of vector.
      !    import :: sp
      !    integer, intent(in) :: n
      !    real(sp), intent(in) :: x(*)
      !    integer, intent(in) :: incx
      ! end function snrm2

      pure real(dp) function dnrm2(n, x, incx)
      !! Euclidean length (L2 Norm) of vector.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
      end function dnrm2

      ! pure subroutine sscal(n, a, x, incx)
      ! !! Vector scale `X = A*X`.
      !    import :: sp
      !    integer, intent(in) :: n
      !    real(sp), intent(in) :: a
      !    real(sp), intent(inout) :: x(*)
      !    integer, intent(in) :: incx
      ! end subroutine sscal

      pure subroutine dscal(n, a, x, incx)
      !! Vector scale `X = A*X`.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: a
         real(dp), intent(inout) :: x(*)
         integer, intent(in) :: incx
      end subroutine dscal

      ! pure integer function isamax(n, x, incx)
      ! !! Find largest component of vector.
      !    import :: sp
      !    integer, intent(in) :: n
      !    real(sp), intent(in) :: x(*)
      !    integer, intent(in) :: incx
      ! end function isamax

      pure integer function idamax(n, x, incx)
      !! Find largest component of vector.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
      end function idamax

   end interface

end module blas_interfaces