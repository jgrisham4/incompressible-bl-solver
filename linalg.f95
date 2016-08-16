      module linalg
        private
        public :: thomas

      contains

        !======================================================
        ! Subroutine for solving a tridiagonal matrix using
        ! Thomas algorithm.
        !======================================================
        subroutine thomas(a,b,c,d,x)
          implicit none
          double precision, dimension(:), intent(inout) :: a,b,c,d
          double precision, dimension(:), intent(out) :: x

          ! Local variables
          integer :: j,nj

          ! Finding size of array
          nj = size(d)

          ! Converting tridiagonal system to upper triangular
          do j=2,nj
            d(j) = d(j) - b(j)/d(j-1)*a(j-1)
            c(j) = c(j) - b(j)/d(j-1)*c(j-1)
          end do

          ! Using back substitution to find the solution
          x(nj) = c(nj)/d(nj) 
          do j=nj-1,1,-1
            x(j) = (c(j) - a(j)*x(j+1))/d(j)
          end do

        end subroutine thomas

      end module linalg
