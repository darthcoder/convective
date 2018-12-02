program heat

implicit none

! define variables
integer, parameter :: I = 80,nmax = I**2    ! The number of divisions
integer:: l,m,n,info,lda,ldb                ! indexing variables
real (kind=8) :: rho = 1.0                  ! density of water
real (kind=8) :: gam = 0.1                  ! diffusion coefficient 
real (kind=8) :: u = 0.1                    ! Velocity field value
real (kind=8) :: Length = 1.0               ! Length of the flow domain
real (kind=8) :: Breadth = 0.5              ! Breadth of the flow domain
real (kind=8) :: del_x, del_y               ! Grid spacing
real (kind=8) :: T_e0, T_w0, T_n0, T_s0     ! Temperature on the boundaries
real (kind=8), dimension(nmax) :: T         ! The temperature vector
real (kind=8), dimension(nmax) :: b         ! The rhs
integer, dimension(nmax) :: ipiv            ! for lapack
real (kind=8), dimension(nmax,nmax) :: A    ! The coefficient matrix     
real (kind=8) :: dxe, dxw, dyn, dys         ! Node distances       
real (kind=8) :: fe,fw,fn,fs                ! Convection strength
real (kind=8) :: de,dw,dn,ds                ! Diffusion strength
real (kind=8), external :: f,d              ! From functions 
character(len=20) :: text_string            ! for file name

! Boundary Values
T_e0 = 50.0!+273.16
T_w0 = 20.0!+273.16
T_s0 = 80.0!+273.16
T_n0 = 40.0!+273.16

! Initialize T, b and A
do m = 1,nmax
    T(m) = 0.0
    b(m) = 0.0
    do l = 1,nmax
        A(l,m) = 0.0
    end do
end do

! Grid Spacings
del_x = Length/real(I)
del_y = Breadth/real(I)

! Discretize the equation on the domain
! We will use the formulation forwarded by Patankar in the following
! We will consider the solution for central, upwind and hybrid methods


n = 1  ! n is the equation index
do l = 1,nmax ! l is from i,j 
    ! now that all the coefficients are defined we need to put them
    ! into the respective matrices. 
    ! we will try to put them in such a manner that the block diagonal
    ! structure of the problem is preserved. 
    ! For that it is necessary that the equations are numbered correctly
    ! we begin with the point l = 1 and keep moving upwards toward the
    ! point  l = I**2
    
    if (mod(l,I).eq.1.and.l.le.I) then  ! ie l=1
        dxw = del_x/2.0
        dys = del_y/2.0
        dxe = del_x
        dyn = del_y
        fe = f(rho,u,del_y)
        fw = f(rho,u,del_y)
        fn = 0.0
        fs = 0.0

        de = d(gam,del_y,dxe)
        dw = d(gam,del_y,dxw)
        dn = d(gam,del_x,dyn)
        ds = d(gam,del_x,dys)
        A(n,l+1) = -de
        A(n,l+I) = -dn
        A(n,l) = fe+de+dw+dn+ds
        b(n) = dw*T_w0+ds*T_s0
        n = n + 1

    else if (((mod(l,I)/=0).or.(mod(l,I)/=1)).and.l.le.I) then ! ie l<i
        dxe = del_x
        dxw = del_x
        dyn = del_y
        dys = del_y/2.0

        fe = f(rho,u,del_y) 
        fw = f(rho,u,del_y)
        fn = 0.0
        fs = 0.0
        
        de = d(gam,del_y,dxe)
        dw = d(gam,del_y,dxw)
        dn = d(gam,del_x,dyn)
        ds = d(gam,del_x,dys)
        
        A(n,l+1) = -de
        A(n,l-1) = -dw
        A(n,l+I) = -dn
        A(n,l) = de+dw+dn+ds+fe
        b(n) =  ds*T_s0
        n = n + 1

    else if (mod(l,I).eq.0.and.l.le.I) then !ie l=i
        dxe = del_x/2.0
        dys = del_y/2.0
        dyn = del_y
        dxw = del_x

        fe = f(rho,u,del_y) 
        fw = f(rho,u,del_y)
        fn = 0.0
        fs = 0.0
        
        de = d(gam,del_y,dxe)
        dw = d(gam,del_y,dxw)
        dn = d(gam,del_x,dyn)
        ds = d(gam,del_x,dys)

        A(n,l-1) = -dw
        A(n,l+I) = -dn
        A(n,l) = de+dw+dn+ds+fe
        b(n) = (de)*T_e0+ds*T_s0
        n = n + 1
    
    else if (mod(l,I).eq.1.and.l.gt.I.and.l.lt.(nmax-I)) then !l=ki+1
        dxw = del_x/2
        dyn = del_y
        dxe = del_x
        dys = del_y

        fe = f(rho,u,del_y) 
        fw = f(rho,u,del_y)
        fn = 0.0
        fs = 0.0
        
        de = d(gam,del_y,dxe)
        dw = d(gam,del_y,dxw)
        dn = d(gam,del_x,dyn)
        ds = d(gam,del_x,dys)
        
        A(n,l+1) = -de
        A(n,l+I) = -dn
        A(n,l-I) = -ds
        A(n,l) = fe+de+dw+dn+ds
        b(n) = (dw)*T_w0
        n = n + 1
    
    else if (((mod(l,I)/=1).or.(mod(l,I)/=0)).and.l.gt.I.and.l.lt.(nmax-I)) then
        dxe = del_x        !l in domain
        dxw = del_x
        dyn = del_y
        dys = del_y

        fe = f(rho,u,del_y) 
        fw = f(rho,u,del_y)
        fn = 0.0
        fs = 0.0
        
        de = d(gam,del_y,dxe)
        dw = d(gam,del_y,dxw)
        dn = d(gam,del_x,dyn)
        ds = d(gam,del_x,dys)

        A(n,l+1) = -de
        A(n,l-1) = -dw
        A(n,l+I) = -dn
        A(n,l-I) = -ds
        A(n,l) = fe+de+dw+dn+ds
        b(n) = 0.0
        n  = n + 1
    
    else if (mod(l,I).eq.0.and.l.gt.I.and.l.le.(nmax-I)) then
        dxe = del_x/2.0     !l = ki
        dxw = del_x
        dyn = del_y
        dys = del_y

        fe = f(rho,u,del_y) 
        fw = f(rho,u,del_y)
        fn = 0.0
        fs = 0.0
        
        de = d(gam,del_y,dxe)
        dw = d(gam,del_y,dxw)
        dn = d(gam,del_x,dyn)
        ds = d(gam,del_x,dys)
        
        A(n,l-1) = -dw
        A(n,l-I) = -ds
        A(n,l+I) = -dn
        A(n,l) = fe+de+dw+dn+ds
        b(n) = (de)*T_e0
        n = n + 1

    else if (mod(l,I).eq.1.and.l.gt.nmax-I) then    ! l=nmax-i+1
        dxw = del_x/2.0			
        dxe = del_x
        dyn = del_y/2.0
        dys = del_y

        fe = f(rho,u,del_y) 
        fw = f(rho,u,del_y)
        fn = 0.0
        fs = 0.0
        
        de = d(gam,del_y,dxe)
        dw = d(gam,del_y,dxw)
        dn = d(gam,del_x,dyn)
        ds = d(gam,del_x,dys)

        A(n,l+1) = -de
        A(n,l-I) = -ds
        A(n,l) = fe+de+dw+dn+ds
        b(n) = (dw)*T_w0+dn*T_n0
        n = n + 1

    else if (((mod(l,I)/=0).and.(mod(l,I)/=1)).and.(l.gt.nmax-I)) then
        dxe = del_x					!l>nmax-i
        dxw = del_x
        dyn = del_y/2.0
        dys = del_y

        fe = f(rho,u,del_y) 
        fw = f(rho,u,del_y)
        fn = 0.0
        fs = 0.0
        
        de = d(gam,del_y,dxe)
        dw = d(gam,del_y,dxw)
        dn = d(gam,del_x,dyn)
        ds = d(gam,del_x,dys)

        A(n,l+1) = -de
        A(n,l-1) = -dw
        A(n,l-I) = -ds
        A(n,l) = de+dw+dn+ds+fe
        b(n) = dn*T_n0
        n = n + 1

    else 
        dxe = del_x/2.0
        dxw = del_x
        dyn = del_y/2.0
        dys = del_y

        fe = f(rho,u,del_y) 
        fw = f(rho,u,del_y)
        fn = 0.0
        fs = 0.0
        
        de = d(gam,del_y,dxe)
        dw = d(gam,del_y,dxw)
        dn = d(gam,del_x,dyn)
        ds = d(gam,del_x,dys)
        
        A(n,l-1) = -dw
        A(n,l-I) = -ds
        A(n,l) = fe+de+dw+dn+ds
        b(n) = (de)*T_e0+dn*T_n0
        n = n + 1
    end if

end do

! do l = 1,nmax
!     print *, (A(l,m), m=1,I**2)
!     print *, b(l)
! end do

! solve using lapack
lda = nmax
ldb = nmax
call dgesv(nmax,1,A,lda,ipiv,b,ldb, info)


print *, b

open(unit=20, file='up80.txt')
write(20,*)b
close(20)

end program heat


real (kind=8) function f(r,u,dz)
    implicit none
    real (kind=8), intent(in):: r,u,dz
    
    f = r*u*dz
end function f

real (kind=8) function d(g,de,dw)
    implicit none
    real (kind=8), intent(in):: g,de,dw
    
    d = g*de/dw
end function d 







