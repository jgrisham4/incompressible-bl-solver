      module incomp
        use linalg
        private
        public :: generate_grid, explicit_soln, crank_nicolson

      contains

        !=========================================================
        ! Simple grid generation subroutine
        !=========================================================

        subroutine generate_grid(dx,dy,imax,jmax,x0,x,y)
          implicit none
          double precision, intent(in) :: dx,dy,x0
          integer, intent(in) :: imax,jmax
          double precision, dimension(:,:), intent(out) :: x,y

          ! Local variables
          integer :: i,j

          ! Generating grid
          do j=1,jmax
            do i=1,imax
              x(i,j) = (dble(i)-1.0)*dx
              y(i,j) = (dble(j)-1.0)*dy
            end do
          end do

          ! Writing grid to file
          open(1, file="grid.dat")
          write(1,*) 'VARIABLES="x","y"'
          write(1,*) 'ZONE I = ',imax,' J=',jmax
          do j=1,jmax
            do i=1,imax
              write(1,*) x(i,j), y(i,j)
            end do
          end do
          close(1)

        end subroutine generate_grid

        !=========================================================
        ! Explicit solution subroutine following Wu's approach.
        !
        ! Arguments: 
        ! dx spacing in the x-direction.
        ! dy spacing in the y-direction.
        ! imax number of points in the x-direction.
        ! jmax number of points in the y-direction.
        ! x array which contains x-coordinates of grid points.
        ! y array which contains y-coordinates of grid points.
        ! nu kinematic viscosity (m^2/s).
        ! Ue edge velocity (m/s).
        ! rho density (kg/m^3).
        ! of name of output file.
        ! p inviscid pressure distribution which must be supplied.
        ! u resulting x-component of the velocity field.
        ! v resulting y-component of the velocity field.
        !=========================================================

        subroutine explicit_soln(dx,dy,imax,jmax,x,y,nu,Ue,rho,of,p,u,v)
          implicit none
          double precision, intent(in) :: dx,dy,nu,Ue,rho
          integer, intent(in) :: imax,jmax
          character (len=20), intent(in) :: of
          double precision, dimension(:,:), intent(in) :: x,y
          double precision, dimension(:), intent(in) :: p
          double precision, dimension(:,:), intent(out) :: u,v

          ! Local variables
          integer :: i,j
          double precision :: Q

          ! Must assign initial velocity profile
          do j=1,jmax
              u(1,j) = Ue
              v(1,j) = 0.0
          end do

          ! Enforcing no slip condition on the wall and Ue on the top
          ! This assumes that the domain is large enough to contain the
          ! whole boundary layer
          do i=1,imax
            u(i,1) = 0.0
            v(i,1) = 0.0
            u(i,jmax) = Ue
            !v(i,jmax) = 0.0
          end do

          ! Marching downstream
          do i=1,imax-1

            ! Finding u along the current slice
            do j=2,jmax-1
              Q = nu*dx/(u(i,j)*dy*dy)
              u(i+1,j) = Q*(u(i,j+1) + u(i,j-1)) &
                - (2.0*Q - 1.0)*u(i,j) &
                - 1.0/(rho*u(i,j))*(p(i+1) - p(i)) &
                - v(i,j)/u(i,j)*dx/dy*0.5*(u(i,j+1) - u(i,j-1))
            end do

            ! Finding v along the current slice
            do j=2,jmax
              v(i+1,j) = v(i+1,j-1) & 
                - dy/(2.0*dx)*(u(i+1,j) - u(i,j) + u(i+1,j-1) - u(i,j-1))
            end do

          end do

          ! Finding the edge of the boundary layer
          do j=1,jmax
            if (u(imax,j)/Ue.gt.0.99) then
              print *, "delta(", x(imax,1), ") = ", y(1,j), " m"
              exit
            end if
          end do

          ! Writing solution to file
          open(2,file=of)
          write(2,*) 'VARIABLES="x","y","u","v"'
          write(2,*) 'ZONE I = ',imax,' J=',jmax
          do j=1,jmax
            do i=1,imax
              write(2,*) x(i,j), y(i,j), u(i,j), v(i,j)
            end do
          end do
          close(2)

          ! Writing final profile to file
          open(3,file="explicit_profile.dat")
          write(3,*) 'VARIABLES="y","u","v"'
          write(3,*) 'ZONE I = ',jmax
          do j=1,jmax
            write (3,*) y(imax,j), u(imax,j), v(imax,j)
          end do
          close(3)

        end subroutine explicit_soln

        !======================================================
        ! Crank-Nicolson method for boundary layer equations.
        ! The approach used below follows the approach outlined
        ! in chapter 7 of the book by Pletcher, Tannehill and
        ! Anderson.
        !======================================================

        subroutine crank_nicolson(dx,dy,imax,jmax,x,y,nu,Ue,rho,of,p,u,v)
          implicit none
          double precision, intent(in) :: dx,dy,nu,Ue,rho
          integer, intent(in) :: imax,jmax
          character (len=20), intent(in) :: of
          double precision, dimension(:,:), intent(in) :: x,y
          double precision, dimension(:), intent(in) :: p
          double precision, dimension(:,:), intent(out) :: u,v

          ! Local variables
          integer :: i,j
          integer :: errA,errB,errC,errD,erru
          double precision, allocatable :: A(:),B(:),C(:),D(:)
          double precision, allocatable :: uip1(:)
          double precision :: dpdx

          ! Allocating memory for diagonal entries
          allocate(A(jmax),stat=errA)
          allocate(B(jmax),stat=errB)
          allocate(C(jmax),stat=errC)
          allocate(D(jmax),stat=errD)
          allocate(uip1(jmax),stat=erru)
          if ((errA.ne.0).or.(errB.ne.0).or.(errC.ne.0).or.(errD.ne.0)) then 
            print *,"ERROR: Couldn't allocate memory in crank_nicolson."
            stop
          end if
          if (erru.ne.0) then
            print *,"ERROR: Couldn't allocate memory in crank_nicolson."
            stop
          end if

          ! Assigning incoming profile
          do j=1,jmax
            u(1,j) = Ue
            v(1,j) = 0.0
          end do

          ! Assigning no slip conditions and top condition
          do i=1,imax
            u(i,1) = 0.0
            v(i,1) = 0.0
            u(i,jmax) = Ue
          end do

          ! Marching in the x-direction 
          do i=1,imax-1
            
            ! Finding pressure gradient for external flow
            dpdx = (p(i+1) - p(i))/dx

            ! Assembling tridiagonal system
            j = 2
            A(j) = -nu/dy**2 + v(i,j)/(2.0*dy)
            B(j) = 0.0
            C(j) = u(i,j)**2/dx - (-nu/dy**2 - v(i,j)/(2.0*dy))*u(i+1,j-1) &
              - 1.0/rho*dpdx
            D(j) = 2.0*nu/dy**2 + u(i,j)/dx
            j = jmax-1
            A(j) = 0.0
            B(j) = -nu/dy**2 - v(i,j)/(2.0*dy)
            C(j) = u(i,j)**2/dx - (-nu/dy**2 + v(i,j)/(2.0*dy))*u(i+1,j+1) &
              - 1.0/rho*dpdx
            D(j) = 2.0*nu/dy**2 + u(i,j)/dx
            do j=3,jmax-2
              A(j) = -nu/dy**2 + v(i,j)/(2.0*dy)
              B(j) = -nu/dy**2 - v(i,j)/(2.0*dy)
              C(j) = u(i,j)**2/dx - 1.0/rho*dpdx
              D(j) = 2.0*nu/dy**2 + u(i,j)/dx
            end do

            ! Solving tridiagonal system for u along current slice
            call thomas(A(2:jmax-1),B(2:jmax-1),C(2:jmax-1),D(2:jmax-1),uip1(2:jmax-1))
            u(i+1,2:jmax-1) = uip1(2:jmax-1)

            ! Finding v
            do j=2,jmax-1
              v(i+1,j) = v(i+1,j-1) &
                - dy/(2.0*dx)*(u(i+1,j)-u(i,j)+u(i+1,j-1)-u(i,j-1))
            end do

          end do

          ! Finding the edge of the boundary layer
          do j=1,jmax
            if (u(imax,j)/Ue.gt.0.99) then
              print *, "delta(", x(imax,1), ") = ", y(1,j), " m"
              exit
            end if
          end do

          ! Writing solution to file
          open(2,file=of)
          write(2,*) 'VARIABLES="x","y","u","v"'
          write(2,*) 'ZONE I = ',imax,' J=',jmax
          do j=1,jmax
            do i=1,imax
              write(2,*) x(i,j), y(i,j), u(i,j), v(i,j)
            end do
          end do
          close(2)

          ! Writing final profile to file
          open(3,file="crank-nicolson-profile.dat")
          write(3,*) 'VARIABLES="y","u","v"'
          write(3,*) 'ZONE I = ',jmax
          do j=1,jmax
            write (3,*) y(imax,j), u(imax,j), v(imax,j)
          end do
          close(3)

        end subroutine crank_nicolson

      end module incomp
