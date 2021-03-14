module hamilton
use types, only: dp
use basis, only: basis_t, derivative
implicit none


contains
function fnm(x,L,n,m) result(res)
    !Kronig-Penney helper function. L = width of well, n,m = indices, x = position of barrier edge
    real(dp), intent(in) :: x,L
    integer, intent(in)  :: n,m
    real(dp)             :: res
    real(dp)             :: pi

    pi = 4*atan(1.0_dp)
    if(n==m) then
        res = x/L - sin(2.0_dp*pi*n*x/L)/(2*pi*n)
    else

        res = sin((m-n)*pi*x/L)/(pi*(m-n)) - sin((m+n)*pi*x/L)/(pi*(m+n))
    end if
end function

function hnm(x,b,L,n,m) result(res)
    !Kronig-Penney helper function. b= width of barrier, n,m=indices, x=position of barrier
    real(dp), intent(in) :: x,b,L
    integer, intent(in)  :: n,m
    real(dp)             :: res

    res = fnm(x+b/2.0_dp,L,n,m) - fnm(x-b/2.0_dp,L,n,m)
    
     
end function
subroutine ext_kronig(H,npeak,V0,a,b,L)
    !The Kronig-Penney potential, with parameters V0, a, b as in the original paper.
	integer, intent(in)	   :: npeak
	real(dp), intent(in)   :: a,b,L,V0
	real(dp), intent(inout):: H(:,:)
	real(dp) 			   :: xi(npeak),pi, start, ends
	integer				   :: i,j,k
    pi = 4*atan(1.0_dp)
	do i=1, npeak
		xi(i) = (-0.5_dp)*a + (i) * a
	end do	
	
	do i=1, size(H(:,1))
	    do j=i, size(H(:,1))
	        do k=1, npeak
                H(i,j) = H(i,j) + V0*hnm(xi(k),b,L,i,j)

	        end do
	        if(i==j) then
                H(i,j) = H(i,j) + (pi**2.0_dp)*(i)**2.0_dp/(L**2.0_dp)
	        end if
	        H(j,i) = H(i,j)
	    end do
	end do

end subroutine
function get_matrix_element(bas,pot,i,j) result(mel)
    !Integral for matrx elements
    type(basis_t)           :: bas
    real(dp),intent(in)     :: pot(:)
    integer, intent(in)     :: i,j 
    real(dp)                :: dx,mel

    dx = bas%x(2)-bas%x(1)

    mel = sum(dx * bas%f(:,i)*pot(:)*bas%f(:,j))
    
end function
subroutine get_zero_energy(H,bas)
    !The zero energy of our basis functions.
    type(basis_t)           :: bas
    real(dp), intent(inout) :: H(:,:)
    integer                 :: i

    do i=1, size(bas%f(1,:))
        H(i,i) = H(i,i) + (4.0_dp*atan(1.0_dp))**2.0_dp*(i)**2.0_dp/(2.0_dp*bas%L**2)
    end do

end subroutine
subroutine get_hartree(H,bas,nx)
    type(basis_t)           :: bas
    real(dp), intent(inout) :: H(:,:)
    real(dp), intent(in)    :: nx(:)
    real(dp)                :: hart(size(bas%x)), temp,lambda
    integer                 :: i,j
    lambda = 1.0_dp
    hart(:) = nx(:)*lambda
    do i=1, size(bas%f(1,:))
        do j=i, size(bas%f(1,:))
            temp = get_matrix_element(bas,hart, i,j)
            H(i,j) = H(i,j)+ temp
            H(j,i) = H(i,j)
        end do
    end do
    close(12)
end subroutine
subroutine get_exchange(H,bas,nx)
    !Exchange potential in one dimension, along with the (rather insignificant) correlation.
    type(basis_t)           :: bas
    real(dp), intent(inout) :: H(:,:)
    real(dp), intent(in)    :: nx(:)
    real(dp)                :: ex(size(bas%x)),correl(size(bas%x)),temp,pi
    integer                 :: i,j
    pi = 4.0_dp*atan(1.0_dp)
    !ex(:) = -(3.0_dp/pi)**(1.0_dp/3.0_dp)*nx(:)**(1.0_dp/3.0_dp)
    ex(:) = -0.5_dp * nx(:)
    correl(:) = -1.0_dp/(2*pi**2.0_dp) * (nx(:)**2.0_dp + 2.0_dp/pi**2.0_dp)/(nx(:)+2.0_dp/pi**2.0_dp)**2.0_dp *&
               (0.829_dp + (1-0.829_dp)*exp(-10.0_dp*nx(:)**2.0_dp)) +&
               20.0_dp * nx(:) * nx(:)**2.0_dp/(nx(:)+2.0_dp/pi**2.0_dp)*(1-0.829_dp)*exp(-10.0_dp*nx(:)**2.0_dp)
    correl(:) = 0.0_dp
    do i=1, size(bas%f(1,:))
        do j=i, size(bas%f(1,:))
            temp = get_matrix_element(bas,ex(:)+correl(:), i,j)
            H(i,j) = H(i,j)+ temp
            H(j,i) = H(i,j)
        end do
    end do
    open(12, file="exchange.dt")
    do i=1, size(bas%x)
        write(12, *) bas%x(i),ex(i)
    end do

    close(12)
end subroutine
function beta_derivative_all(bas,rhob,db) result(drhoa)
    !Numerical beta derivative (Nagy 2011, 2008)
    type(basis_t)           :: bas
    real(dp)                :: rhob(:,:), drhoa(size(rhob(:,1)),size(rhob(1,:))  ),db
    real(dp)                :: dermat(size(rhob(:,1)),size(rhob(:,1)))
    integer                 :: i
    dermat = 0.0_dp
    do i=1, 2
        drhoa(i,:) = -11.0_dp*rhob(i,:) + 18.0_dp*rhob(i+1,:) - 9.0_dp * rhob(i+2,:) + 2.0_dp*rhob(i+3,:)
        drhoa(i,:) = drhoa(i,:)/(6.0_dp*db)
    end do
    do i=3, 3
        drhoa(i,:)= 1.0_dp * rhob(i-2,:) - 8.0_dp * rhob(i-1,:) + 8.0_dp * rhob(i+1,:) - 1.0_dp * rhob(i+2,:)
        drhoa(i,:)= drhoa(i,:)/(12.0_dp*db)
    end do
    do i=4, 5
        drhoa(i,:)= -2.0_dp * rhob(i-3,:) + 9.0_dp * rhob(i-2,:)-18.0_dp*rhob(i-1,:)+11.0_dp*rhob(i,:)
        drhoa(i,:)=drhoa(i,:)/(6.0_dp*db)
    end do

end function
function beta_derivative_exact(bas,energies,db,rho) result(betaterm)
    !Calculates the beta/gamma terms exactly by solving a system of linear equations, as per Nagy's paper
    !. J Chem Phys 135 ("Functional derivative..") eq 21 
    type(basis_t):: bas
    real(dp)     :: rho(:), energies(:),betaterm(size(energies)),db
    real(dp)     :: orbs(size(energies)), A(size(energies),size(energies)), Aw(size(energies),size(energies))
    real(dp)     :: betagrid(size(energies)), B(size(energies))
    integer      :: LDA, LDB, N, NRHS,i, IPIV(size(energies)), INFO
    do i=1, size(energies)
        betagrid(i) = (i-1)*db
    end do
    B(:) = rho(:)
    N = size(energies)
    NRHS = 1
    LDA = N
    LDB = N
    do i=1, size(energies)
        A(i,:) = 2.0_dp*exp(betagrid(i)*energies(:))
    end do
    Aw(:,:)=A(:,:)
    call dgesv(N, NRHS, Aw, LDA, IPIV, B, LDB, INFO)
    if (INFO /= 0) then
        print *, "dgesv failed with code ", INFO
        stop INFO
    end if
    orbs(:) = B(:)

    do i=1, size(energies)
        betaterm(i) = sum(energies(:)*A(i,:)*orbs(:))
    end do
    
end function
subroutine get_pauli_appro(H,bas,nx,mu)
    !Approximate Pauli potential developed for thesis
    type(basis_t),intent(in):: bas
    real(dp),intent(inout)  :: H(:,:),nx(:),mu
    real(dp)                :: ex(size(bas%x)), hart(size(bas%x)), pi, dx, con, ef, e0, a(size(bas%x)),b(size(bas%x))
    real(dp)                :: pot(size(bas%x)), vp(size(bas%x)), f(size(bas%x)), dnx(size(bas%x)),beta(size(bas%x)),lambda
    integer                 :: i,j
    lambda = 1.0_dp
    pi = 4.0_dp*atan(1.0_dp)
    dx = bas%x(2)-bas%x(1)
    ex = 0.0_dp
    dnx = derivative(bas,nx)
    hart = 0.0_dp
    con = 2.0_dp*sqrt(2.0_dp)/(3.0_dp*pi)
    !ef = 52.52_dp
    !e0 = 22.835_dp
    e0 = -0.14_dp
    ef = 9.01e-2_dp
    goto 10
    do i=1, size(bas%x)
        hart(i) = sum(dx*nx(:)/sqrt( (bas%x(i)-bas%x(:))**2.0_dp + 0.1_dp  ))
    end do
    10 continue
    hart(:) = lambda * nx(:)
    ex(:) = -(3.0_dp/pi)**(1.0_dp/3.0_dp)*nx(:)**(1.0_dp/3.0_dp)
    pot = ex+hart
    a(:) = ef - pot(:)
    b(:) = e0 - pot(:)
    do i=1, size(bas%x)
        if (a(i) < 0.0_dp) a(i) = 0.0_dp
        if (b(i) < 0.0_dp) b(i) = 0.0_dp
    end do
    beta = con * (sqrt(a(:))*(2.0_dp*pot(:)+ef) - sqrt(b(:))*(2.0_dp*pot(:) + e0))&
         + e0*1.0_dp/2.0_dp * sin(pi/4.0_dp * bas%x(:))**2.0_dp
    vp(:) = ef
    vp(:) = vp(:)-2.0_dp*beta(:)/nx(:)
    f = 0.0_dp
    do i=1, size(bas%x)
        do j=1, size(bas%x)-i+1
            f(i) = f(i)-  dnx(i+j-1)*beta(i+j-1)*dx
        end do
    end do
    vp(:) = vp(:) + 2.0_dp*f(:)/(nx(:)**2.0_dp)
    vp(1:150) = vp(150)
    do i=1,999
        vp(i+1001) = vp(1000-i)
    end do
    vp(1001) = vp(1000)
    open(13, file="approvp.dt")
    do i=1,size(bas%x)
        write(13, *) bas%x(i),vp(i) 
    end do
    close(13)
    do i=1, size(bas%f(1,:))
        do j=i, size(bas%f(1,:))
            H(i,j) = H(i,j)+ get_matrix_element(bas,vp(:), i,j)
            H(j,i) = H(i,j)
        end do
    end do   
end subroutine
subroutine get_pauli_exact(H,bas,psi,energies,nx)
    !Exact pauli potential from Nagy 2008
    type(basis_t),intent(in):: bas
    real(dp),intent(inout)  :: H(:,:)
    real(dp),intent(in)     :: nx(:),energies(:),psi(:,:)
    real(dp)                :: vp(size(bas%x)),vp2(size(bas%x)), f(size(bas%x))
    real(dp)                :: pi,dx, eqenergies(bas%nelectrons/2),beta(size(bas%x))
    real(dp)                :: dbeta(size(bas%x)), dnx(size(bas%x))
    integer                 :: i,j
    pi = 4.0_dp*atan(1.0_dp)
    dx = bas%x(2)-bas%x(1)
    vp(:) = 0.0_dp
    f(:) = 0.0_dp
    beta(:) = 0.0_dp
    dnx = derivative(bas,nx(:))

    do i=1, bas%nelectrons/2
        beta(:) = beta(:) + psi(:,i)**2.0_dp*2.0_dp*energies(i)
    end do
     
    print *, energies(size(energies))
    dbeta = derivative(bas,beta)
    f(:) =  - dbeta(:) + energies(size(energies))*dnx(:)
    !f(:) = !energies(size(energies))*dnx(:)!f(:)
    do i=1, size(bas%x)
        do j=1, size(bas%x)-i+1
            vp(i) = vp(i)- nx(i+j-1)*f(i+j-1)*dx
        end do
    end do
    !vp(:) = 0.0_dp
    !vp(:) = nx(:)*beta(:)*0.5_dp
    open(12, file ="vp3.dt")
    do i=1, size(bas%x)
        write(12, *) bas%x(i), beta(i)
    end do
    
    close(12)
    vp(:) = 2.0_dp/(nx(:)+0.001_dp)**2.0_dp * vp(:)

    do i=1, size(bas%f(1,:))
        do j=i, size(bas%f(1,:))
            H(i,j) = H(i,j)+ get_matrix_element(bas,vp(:), i,j)
            H(j,i) = H(i,j)
        end do
    end do 
    open(12, file ="vp4.dt")
    do i=1, size(bas%x)
        write(12, *) bas%x(i), vp(i)
    end do
    close(12)
    
end subroutine
subroutine get_pauli_exact_symmetry(H,bas,psi,energies,nx)
    !Exact Pauli potential with symmetry enforced
    type(basis_t),intent(in):: bas
    real(dp),intent(inout)  :: H(:,:)
    real(dp),intent(in)     :: nx(:),energies(:),psi(:,:)
    real(dp)                :: vp(size(bas%x)),vp2(size(bas%x)), f(size(bas%x))
    real(dp)                :: pi,dx, eqenergies(bas%nelectrons/2),beta(size(bas%x))
    real(dp)                :: dbeta(size(bas%x)), dnx(size(bas%x))
    integer                 :: i,j
    pi = 4.0_dp*atan(1.0_dp)
    dx = bas%x(2)-bas%x(1)
    vp(:) = 0.0_dp
    f(:) = 0.0_dp
    beta(:) = 0.0_dp
    dnx = derivative(bas,nx(:))

    do i=1, bas%nelectrons/2
        beta(:) = beta(:) + psi(:,i)**2.0_dp*2.0_dp*energies(i)
    end do
     
    print *, energies(size(energies))
    dbeta = derivative(bas,beta)
    f(:) =  - dbeta(:) + energies(size(energies))*dnx(:)
    !f(:) = !energies(size(energies))*dnx(:)!f(:)
    do i=1, size(bas%x)
        do j=1, size(bas%x)-i+1
            vp(i) = vp(i)- nx(i+j-1)*f(i+j-1)*dx
        end do
    end do
    !vp(:) = 0.0_dp
    !vp(:) = nx(:)*beta(:)*0.5_dp
    open(12, file ="vp3.dt")
    do i=1, size(bas%x)
        write(12, *) bas%x(i), beta(i)
    end do
    
    close(12)
    vp(:) = 2.0_dp/(nx(:)+0.001_dp)**2.0_dp * vp(:)
    do i=1,999
        vp(i+1001) = vp(1000-i)
    end do
        vp(1001) = vp(1000)
    do i=1, size(bas%f(1,:))
        do j=i, size(bas%f(1,:))
            H(i,j) = H(i,j)+ get_matrix_element(bas,vp(:), i,j)
            H(j,i) = H(i,j)
        end do
    end do 
    open(12, file ="vp4.dt")
    do i=1, size(bas%x)
        write(12, *) bas%x(i), vp(i)
    end do
    close(12)
    
end subroutine
subroutine get_pauli_beta_exact(H,bas,energies,db,nx)
    !Pauli term calculated with exact Beta term, but with partial integration in place.
    real(dp),intent(inout)  :: H(:,:,:)
    real(dp),intent(in)     :: nx(:,:),db,energies(:)
    type(basis_t),intent(in):: bas
    real(dp)                :: dbeta(size(nx(:,1)),size(nx(1,:))),beta(size(nx(:,1)),size(nx(1,:)))
    real(dp)                :: f(size(nx(:,1)),size(nx(1,:))), vp(size(nx(:,1)),size(nx(1,:))),dx,k1,k2,k3,k4
    real(dp)                :: dnx(size(nx(:,1)),size(nx(1,:)))
    integer                 :: i,j,k
    !db = 0.01_dp
    f(:,:) = 0.0_dp
    vp(:,:) = energies(size(energies))
    dx = bas%x(2)-bas%x(1)
    do i=1, size(nx(:,1))
        dnx(i,:) = derivative(bas,nx(i,:))
    end do
    !beta(:) = 0.0_dp
    dbeta(:,:) = 0.0_dp
    do i=1,size(bas%x)
        beta(:,i)  = beta_derivative_exact(bas,energies,db,nx(:,i))
    end do
    do k=1, size(nx(:,1))
        vp(k,:) = vp(k,:)  -2.0_dp*beta(k,:)/(nx(k,:))
    end do
    f(:,:) = 0.0_dp
    do k=1, size(nx(:,1))
        do i=1, size(bas%x)
            do j=1, size(bas%x)-i+1
                f(k,i) = f(k,i)-  dnx(k,i+j-1)*beta(k,i+j-1)*dx
            end do
        end do
    end do
    do k=1, size(nx(:,1))
        vp(k,:) = vp(k,:) + 2.0_dp*f(k,:)/(nx(k,:))**2.0_dp
    end do
    do k=1,size(nx(:,1))
    do i=1,999
        vp(k,i+1001) = vp(k,1000-i)
    end do
        vp(k,1001) = vp(k,1000)
    end do
    open(12, file="brho-exact.dt")
    do i=1, size(nx(1,:))
        write(12, * ) bas%x(i), vp(1,i),vp(2,i),vp(3,i),vp(4,i),vp(5,i),vp(6,i),vp(7,i)
    end do
    close(12)
    
    open(12, file="beta-exact.dt")
    do i=1, size(nx(1,:))
        write(12, * ) bas%x(i), beta(1,i)
    end do
    close(12)
    do k=1, size(nx(:,1))
        do i=1, size(bas%f(1,:))
            do j=i, size(bas%f(1,:))
                H(k,i,j) = H(k,i,j)+ get_matrix_element(bas,vp(k,:), i,j)
                H(k,j,i) = H(k,i,j)
            end do
        end do
    end do

end subroutine

subroutine get_pauli_beta_rho(H,bas,mu,db,nx)
    !Pauli term calculated with approximate numerical beta derivative.
    real(dp),intent(inout)  :: H(:,:,:)
    real(dp),intent(in)     :: mu,nx(:,:),db
    type(basis_t),intent(in):: bas
    real(dp)                :: dbeta(size(nx(:,1)),size(nx(1,:))),beta(size(nx(:,1)),size(nx(1,:)))
    real(dp)                :: f(size(nx(:,1)),size(nx(1,:))), vp(size(nx(:,1)),size(nx(1,:))),dx,k1,k2,k3,k4
    real(dp)                :: dnx(size(nx(:,1)),size(nx(1,:)))
    integer                 :: i,j,k
    !db = 0.01_dp
    f(:,:) = 0.0_dp
    vp(:,:) = mu
    dx = bas%x(2)-bas%x(1)
    do i=1, size(nx(:,1))
        dnx(i,:) = derivative(bas,nx(i,:))
    end do
    !beta(:) = 0.0_dp
    dbeta(:,:) = 0.0_dp

    beta(:,:)  = beta_derivative_all(bas,nx,db)
    do k=1, size(nx(:,1))
        vp(k,:) = vp(k,:)  -2.0_dp*beta(k,:)/(nx(k,:))
    end do
    f(:,:) = 0.0_dp
    do k=1, size(nx(:,1))
        do i=1, size(bas%x)
            do j=1, size(bas%x)-i+1
                f(k,i) = f(k,i)-  dnx(k,i+j-1)*beta(k,i+j-1)*dx
            end do
        end do
    end do
    do k=1, size(nx(:,1))
        vp(k,:) = vp(k,:) + 2.0_dp*f(k,:)/(nx(k,:))**2.0_dp
    end do
    do k=1,size(nx(:,1))
        vp(k,1:150) = vp(k,150)
        vp(k,1850:2000) = vp(k,2000)
    end do

    do k=1,size(nx(:,1))
    do i=1,999
        vp(k,i+1001) = vp(k,1000-i)
    end do
        vp(k,1001) = vp(k,1000)
    end do
    open(12, file="brho.dt")
    do i=1, size(nx(1,:))
        write(12, * ) bas%x(i), vp(1,i),vp(2,i),vp(3,i),vp(4,i),vp(5,i)
    end do
    close(12)
    
    open(12, file="beta.dt")
    do i=1, size(nx(1,:))
        write(12, * ) bas%x(i), beta(1,i),beta(2,i)
    end do
    close(12)
    do k=1, size(nx(:,1))
        do i=1, size(bas%f(1,:))
            do j=i, size(bas%f(1,:))
                H(k,i,j) = H(k,i,j)+ get_matrix_element(bas,vp(k,:), i,j)
                H(k,j,i) = H(k,i,j)
            end do
        end do
    end do

end subroutine


subroutine get_pauli_ready(H,bas,vp)
    real(dp), intent(inout) :: H(:,:)
    type(basis_t)           :: bas
    real(dp), intent(in)    :: vp(:)
    integer                 :: i,j,k
    do i=1, size(bas%f(1,:))
        do j=i, size(bas%f(1,:))
            H(i,j) = H(i,j)+ get_matrix_element(bas,vp(:), i,j)
            H(j,i) = H(i,j)
        end do
    end do 

end subroutine
subroutine get_pauli_continuous_appro(H,bas,nx)
    !Approximation of the Pauli potential based on the idea beta = g(x)rho(x)
    real(dp), intent(inout) :: H(:,:)
    type(basis_t)           :: bas
    real(dp)                :: nx(:),vp(size(nx)),f(size(nx)), g(size(nx)), nxsq(size(nx))
    real(dp)                :: dnxsq(size(nx)),ddnxsq(size(nx)), dnx(size(nx))
    real(dp)                :: dx
    integer                 :: i,j
    dx = bas%x(2)-bas%x(1)
    nxsq(:) = sqrt(nx(:)) 
    dnxsq(:) = derivative(bas,nxsq(:))
    ddnxsq(:) = derivative(bas,dnxsq(:))
    dnx(:) = derivative(bas,nx(:))
    g(:) = -1/sqrt(nx) * 1.0_dp/2.0_dp *ddnxsq(:) - nx(:)**2.0_dp * (4.0_dp*atan(1.0_dp) )**2.0_dp/6.0_dp + 8.97_dp
    g(:) = 0.5_dp *g(:)
    f(:) = g(:) * derivative(bas,nx)**2.0_dp
    vp(:) = 0.0_dp
    do i=1, size(bas%x)
        do j=1, size(bas%x)-i+1
            vp(i) = vp(i)-  f(i+j-1)*dx
        end do
    end do
    vp(:) = -2.0_dp/nx(:)**2.0_dp * vp(:) !+ 2*g(:) 
    open(12, file="cont-vp.dt")
    do i=1, size(nx(:))
        write(12, * ) bas%x(i), vp(i), 2*g(i)*dnx(i)/nx(i)**2
    end do
    close(12)
end subroutine
subroutine get_thomas_fermi(H,bas,nx)
    !Thomas-Fermi potential for OFDFT calculations with the Euler equation.
    type(basis_t)   :: bas
    real(dp)        :: H(:,:),nx(:), pot(size(nx)), drho(size(nx)), drho2(size(nx)),pi
    integer         :: j,i
    pi = 4*atan(1.0_dp)
    pot(:) = pi**2.0_dp /2.0_dp * nx(:)**2.0_dp
    drho = derivative(bas,nx)
    drho2 = derivative(bas,drho)

    !do i=1,999
    !pot(i+1001) = pot(1000-i)
    !end do
    !pot(1001) = pot(1000)
    !pot(:) = pot(:) - (1.0_dp/8.0_dp * (drho*drho)/((nx(:))**2.0_dp) &
    !        - 1.0_dp/4.0_dp * drho2(:)/(nx(:)))
    do i=1,999
        pot(i+1001) = pot(1000-i)
    end do
    pot(1001) = pot(1000)
    open(11,file="tfvp.dt")
    do i=1, size(nx)
        write(11, *) bas%x(i), pot(i),nx(i),drho2(i),drho(i)
    end do

    close(11)
    do i=1, size(bas%f(1,:))
        do j=i, size(bas%f(1,:))
            H(i,j) = H(i,j)+ get_matrix_element(bas,pot(:), i,j)
            H(j,i) = H(i,j)
        end do
    end do     
    

    
end subroutine
end module hamilton