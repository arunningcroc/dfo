module energies
use types, only: dp
use basis, only: basis_t, derivative, expand_psi_ofdft
implicit none


contains

function get_kinetic_energy(bas,psi) result(ke)
    real(dp), intent(in) :: psi(:,:)
    type(basis_t)        :: bas
    real(dp)             :: ke, dx, step
    real(dp),allocatable :: d2s(:,:),dpsi(:)
    integer              :: i,nel

    nel = bas%nelectrons/2

    allocate(d2s(size(bas%x),nel), dpsi(size(bas%x)))
    !Calculate second derivatives
    dpsi = 0.0_dp
    do i=1, nel
        dpsi(:)  = derivative(bas, psi(:,i))
        d2s(:,i) = derivative(bas, dpsi(:))
    end do
    dx = bas%x(2)-bas%x(1)
    !Calculate the integrals -1/2 * psi * nabla^2 * psi
    ke = 0.0_dp

    do i=1, nel
        step = sum(psi(:,i)*d2s(:,i)*dx)
        ke = ke - 0.5_dp * step
    end do

    deallocate(d2s,dpsi)
end function

function get_delta_energy(bas,nx) result(en)
    type(basis_t)           :: bas
    real(dp), intent(in)    :: nx(:)
    real(dp)                :: en, dx
    dx = bas%x(2)-bas%x(1)
    en = sum(nx(:)**2 * 0.5_dp * dx)
    

end function
function get_exchange_energy(bas,nx) result(ex)
    type(basis_t)       :: bas
    real(dp),intent(in) :: nx(:)
    real(dp)            :: ex, pi, dx
    pi = 4.0_dp*atan(1.0_dp)
    dx = bas%x(2)-bas%x(1)
    ex = -0.25_dp* sum(nx(:)**2.0_dp   *dx )! - 1.0_dp/(2.0_dp*pi**2.0_dp)*sum(dx*(nx(:)**2.0_dp)/(nx(:)+2.0_dp/pi**2.0_dp) *&
        !(0.829_dp + (1.0_dp-0.829_dp))*exp(-10.0_dp*nx(:)**2.0_dp))
        

    
end function

function get_nagy_total_energy(bas,alpha,sn,vec,mu) result(k)
    type(basis_t)           :: bas
    real(dp)                :: k, g(size(bas%x)), beta(size(bas%x)),pi,alpha,mu,rho(size(bas%x)), drho(size(bas%x))
    real(dp)                :: sn(:),vec(:), mupar, psi(size(bas%x)),dx, xint(size(bas%x)), nonloc(size(bas%x))
    integer                 :: i,j
    !Beta term for non-interacting system
    beta = 0.0_dp
    !Rho for non-interacting system
    rho = 0.0_dp
    !x-integral of the non local part
    xint = 0.0_dp
    !Non local part of the potential
    nonloc = 0.0_dp
    mupar = 10.0_dp

    dx = bas%x(2)-bas%x(1)
    pi = 4*atan(1.0_dp)
    psi = expand_psi_ofdft(bas,vec(:)+alpha*sn(:))
    do i=1, bas%nelectrons/2
        beta(:) = beta(:) + 2.0_dp*(i**2.0_dp * pi**2.0_dp/(2*bas%L))*bas%f(:,i)**2.0_dp
    end do
    do i=1, bas%nelectrons/2
        rho(:) = rho(:) + 2.0_dp * bas%f(:,i)**2.0_dp
    end do
    g(:) = beta(:)/rho(:)

    drho(:) = derivative(bas,rho)
    nonloc(:) = sum(derivative(bas,psi(:)**2.0_dp)*psi(:)**2.0_dp * bas%dx)/psi(:)**2.0_dp
    
    
    !goto 10

    goto 20
    open(12, file="gnonloc.dt")
    
    do i=1,size(bas%x)
        write(12,*) bas%x(i), nonloc(i),g(i),drho(i),psi(i)**2.0_dp
    end do
    close(12)
    stop
    20 continue
    k = sum(( abs( derivative(bas,psi(:)**2.0_dp) ) )**2.0_dp/(psi(:)**2.0_dp)*dx)*1.0_dp/(8.0_dp)
    k = k+mupar*sum(psi(:)**2.0_dp*dx) - 2.0_dp * sum(g(:)*psi(:)**2.0_dp*dx) &
        +sum(0.25_dp*psi(:)**4.0_dp*dx)-mu*sum(psi(:)**2.0_dp*dx)-sum(nonloc(:)*bas%dx)
    !print *, "here"
    !print *, k 
    !stop
end function
function get_nagy_gradient(bas,vec,mu) result(grad)
    type(basis_t)           :: bas
    real(dp)                :: vec(:),mu,psi(size(bas%x)),grad(bas%nfuns),xint(size(bas%x)),drho(size(bas%x)),mupar
    real(dp)                :: nonloc(size(bas%x),bas%nfuns),rho(size(bas%x)),beta(size(bas%x)),dx,g(size(bas%x)),pi
    integer                 :: i,j,k,z
    !Beta term for non-interacting system
    beta = 0.0_dp
    !Rho for non-interacting system
    rho = 0.0_dp
    !x-integral of the non local part
    xint = 0.0_dp
    !Non local part of the potential
    nonloc = 0.0_dp
    grad = 0.0_dp
    mupar = 10.0_dp
    dx = bas%x(2)-bas%x(1)
    pi = 4*atan(1.0_dp)
    psi = expand_psi_ofdft(bas,vec(:))
    dx = bas%x(2)-bas%x(1)
    do i=1, bas%nelectrons/2
        beta(:) = beta(:) + 2.0_dp*(i**2.0_dp * pi**2.0_dp/(2*bas%L))*bas%f(:,i)**2.0_dp
    end do
    do i=1, bas%nelectrons/2
        rho(:) = rho(:) + 2.0_dp * bas%f(:,i)**2.0_dp
    end do
    g(:) = beta(:)/rho(:)
    
    drho(:) = derivative(bas,rho)
    nonloc = 0.0_dp
    do i=1, bas%nfuns
        nonloc(:,i) = nonloc(:,i)+ 2.0_dp*sum((bas%fder(:,i)*psi(:)+&
                      bas%f(:,i)*derivative(bas,psi(:)) )*psi(:)**2.0_dp*bas%dx)/psi(:)**2.0_dp
        nonloc(:,i) = nonloc(:,i) + sum(bas%f(:,i)*derivative(bas,psi(:)**2.0_dp)*bas%dx)/psi(:)**2.0_dp
        nonloc(:,i) = nonloc(:,i) - 2.0_dp*bas%f(:,i)*sum(psi(:)**2.0_dp*derivative(bas,psi(:)**2.0_dp)*bas%dx)/psi(:)**3.0_dp
    end do
    !grad(:) = get_vw_gradient
    !goto 10

    !10 continu
    do i=1, bas%nfuns
        !print *, i
        grad(i) = grad(i) + sum(2.0_dp*mupar * bas%f(:,i) *psi(:)* dx)-sum(4.0_dp*g(:)*bas%f(:,i)*psi(:)*dx)
        grad(i) = grad(i) + sum(bas%f(:,i)*psi(:)**3.0_dp*dx)-2.0_dp*mu*sum(bas%f(:,i)*psi(:)*dx)
        
        grad(i) = grad(i) + 1.0_dp/(8.0_dp) *&
        sum(bas%dx * 4.0_dp * (bas%fder(:,i)*psi(:)+bas%f(:,i)*derivative(bas,psi))*&
        (2.0_dp*derivative(bas,psi)*psi)/psi(:)**2.0_dp )
        grad(i) = grad(i) - 2.0_dp/(8.0_dp) *&
         sum(bas%dx * bas%f(:,i)*abs(derivative(bas,psi**2.0_dp))**2.0_dp/(psi**3.0_dp))

        grad(i) = grad(i) - sum(nonloc(:,i)*bas%dx)
    end do
    !print *, "we here?"
end function
function get_vw_total_energy(bas,alpha,sn,vec,mu) result(k)
    type(basis_t)           :: bas
    real(dp)                :: sn(:),vec(:),alpha,beta, k
    real(dp)                :: psi(size(bas%x)),dx, mu
    k = 0.0_dp
    dx = bas%x(2)-bas%x(1)
    k = get_thomas_fermi_total_energy(bas,alpha,sn,vec,mu)
    psi = expand_psi_ofdft(bas,vec(:)+alpha*sn(:))
    
    k = k + sum(( abs( derivative(bas,psi(:)**2.0_dp) ) )**2.0_dp/(psi(:)**2.0_dp)*dx)*1.0_dp/(8.0_dp*9.0_dp)

end function
function get_vw_gradient(bas,vec,mu) result(grad)
    type(basis_t)           :: bas
    real(dp)                :: vec(:), mu, grad(bas%nfuns), psi(size(bas%x))
    integer                 :: i

    grad(:) = get_thomas_fermi_gradient(bas,vec,mu)
    psi(:) = expand_psi_ofdft(bas,vec)
    do i=1, bas%nfuns
        grad(i) = grad(i) + 1.0_dp/(8.0_dp*9.0_dp) *&
        sum(bas%dx * 4.0_dp * (bas%fder(:,i)*psi(:)+bas%f(:,i)*derivative(bas,psi))*&
        (2.0_dp*derivative(bas,psi)*psi)/psi(:)**2.0_dp )
        grad(i) = grad(i) - 2.0_dp/(8.0_dp*9.0_dp) *&
         sum(bas%dx * bas%f(:,i)*abs(derivative(bas,psi**2.0_dp))**2.0_dp/(psi**3.0_dp))
    end do
end function
function get_thomas_fermi_total_energy(bas,alpha,sn,vec,mu) result(k)
    type(basis_t)           :: bas
    real(dp)                :: sn(:),vec(:),k,alpha,beta, kindens(size(bas%x)), deldens(size(bas%x)),mu,const
    real(dp)                :: psisn(size(bas%x)),pi
    integer                 :: i,j
    pi = 4*atan(1.0_dp)
    const = pi**2.0_dp/6.0_dp
    beta = 0.25_dp
    psisn = expand_psi_ofdft(bas,vec(:)+alpha*sn(:))
    deldens(:) = beta*psisn(:)**4.0_dp
    kindens(:) = psisn(:)**6.0_dp * const
    k = sum(deldens+kindens)-mu*sum(psisn(:)**2.0_dp)
    k = k*(bas%x(2)-bas%x(1))
end function
function get_thomas_fermi_gradient(bas,vec,mu) result(grad)
    type(basis_t)   :: bas
    real(dp)        :: mu, vec(:), psi(size(bas%x)),grad(bas%nfuns),alpha,beta,pi
    integer         :: i
    psi = 0.0_dp
    pi = 4*atan(1.0_dp)
    alpha = pi**2.0_dp/6.0_dp
    beta = 0.25_dp
    psi(:) = expand_psi_ofdft(bas,vec(:))

    do i=1, bas%nfuns
        grad(i) = sum(  6*alpha*bas%f(:,i)*psi(:)**5.0_dp+4.0_dp*beta*bas%f(:,i)*psi(:)**3.0_dp )-2*mu*sum(bas%f(:,i)*psi)
        grad(i) = grad(i) * (bas%x(2)-bas%x(1))
    end do
end function
end module







