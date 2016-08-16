       !=======================================================
       ! This code solves the incompressible boundary layer
       ! equations using one of two options:
       !
       ! (1) an explicit finite difference formulation from 
       !     Wu (1961).  
       ! (2) a Crank-Nicolson scheme.
       !
       ! The options can be changed by commenting out the 
       ! calls to the subroutines on lines 91 and 92. 
       !
       ! Notes: 
       ! - This code currently assumes a constant value of Ue.
       ! - A utility script for comparing with the Blasius 
       !   solution is provided.  It was written in Python 
       !   and is named comparison-blasius.py.
       !
       ! Author: James Grisham
       ! Date: 02/06/2016
       !=======================================================

       program solver
         use incomp
         implicit none
         double precision, allocatable :: x(:,:),y(:,:),u(:,:),v(:,:)
         double precision, allocatable :: y0(:), u0(:)
         double precision, allocatable :: p(:)
         double precision :: Ue,mu,nu,rho,p_e
         double precision :: dx,dy,xmax,ymax,x0
         double precision :: cf,q,dudy,cf_exact
         integer :: imax,jmax,npts,i,j,aerr1,aerr2,mx,my
         character (len=20) :: ex_file, im_file, output_file

         ! Defining inputs
         nu = 15.61434409e-6  ! m^2/s
         rho = 1.22           ! kg/m^3
         mu = nu*rho          ! kg/(m-s)
         Ue = 5.0             ! m/s
         p_e = 101.325e3      ! Pa
         ex_file = "explicit.tec"
         im_file = "crank-nicolson.tec"
         !im_file = "mesh2.tec"
         output_file = "mesh2.out"

         ! Grid inputs (assuming the grid starts at x=x0, y=0)
         mx = 50*4
         my = 5*4
         imax = 4*mx + 1
         jmax = 4*my + 1
         xmax = 0.3
         ymax = 0.01
         x0 = 0.0

         ! Printing some information
         print *, "Reynolds number: ", rho*Ue*xmax/mu

         ! Setting inviscid pressure
         allocate(p(imax),stat=aerr1)
         if (aerr1.ne.0) then
           print *, "ERROR: Couldn't allocate memory for p."
           stop
         end if
         do i=1,imax
           p(i) = p_e
         end do

         ! Finding grid spacings
         dx = xmax/dble(imax-1)
         dy = ymax/dble(jmax-1)

         ! Allocating memory for the grid
         allocate(x(imax,jmax),stat=aerr1)
         allocate(y(imax,jmax),stat=aerr2)
         if ((aerr1.ne.0) .or. (aerr2.ne.0)) then
           print *, "ERROR: Couldn't allocate memory."
           stop
         end if

         ! Allocating memory for the solution
         allocate(u(imax,jmax),stat=aerr1)
         allocate(v(imax,jmax),stat=aerr2)
         if ((aerr1.ne.0) .or. (aerr2.ne.0)) then
           print *, "ERROR: Couldn't allocate memory."
           stop
         end if

         ! Generating the grid
         call generate_grid(dx,dy,imax,jmax,x0,x,y)

         ! Solving the problem using FDM
         !call explicit_soln(dx,dy,imax,jmax,x,y,nu,Ue,rho,ex_file,p,u,v)
         call crank_nicolson(dx,dy,imax,jmax,x,y,nu,Ue,rho,im_file,p,u,v)

         ! Computing cf at the last slice
         cf_exact = 0.664/sqrt(Ue*xmax/nu)
         q = 0.5*rho*Ue**2
         dudy = (-3.0*u(imax,1) + 4.0*u(imax,2) - u(imax,3))/(2.0*dy)
         cf = mu*dudy/q
         open(5,file=output_file)
         write(5,*) (x(2,1) - x(1,1))*(y(1,2)-y(1,1)), cf
         close(5)


       end program solver
