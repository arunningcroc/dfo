module basis
use types,only: dp
implicit none

type, public:: basis_t 
    integer              :: nfuns, nelectrons
    real(dp),allocatable :: x(:), f(:,:), fder(:,:)
    real(dp)             :: L,dx
    contains
        final            :: destroy_basis
    
end type

interface Basis_t
    module procedure :: create_basis
end interface

contains

function create_basis(nelectrons,nfuns,ngrid,L) result(bas)
	type(Basis_t)		        :: bas
    integer,intent(in)          :: nfuns,ngrid,nelectrons
    real(dp),intent(in)         :: L
    integer                     :: i,j
    real(dp),allocatable        :: f(:,:),x(:),fder(:,:)
    real(dp)                    :: pi,ep
    pi = 4*atan(1.0_dp)
    allocate(f(ngrid,nfuns),x(ngrid),fder(ngrid,nfuns))
    ep = (L-L/ngrid)/(ngrid)
	do i=1, ngrid
		x(i) = (i) * ep
	end do
    do i=1, nfuns
        f(:,i) = sqrt(2/L) * sin(i*pi*x(:)/L)!+1.0_dp
        fder(:,i) = sqrt(2/L) * i*pi/L *cos(i*pi*x(:)/L)
    end do
    bas%nfuns = nfuns
    bas%f = f
    bas%x = x
    bas%L = L
    bas%dx= x(2)-x(1)
    bas%fder = fder
    bas%nelectrons = nelectrons
end function

subroutine destroy_basis(this)
    type(Basis_t) :: this

    deallocate(this%f, this%x)
end subroutine

function expand_rho(bas,vec,of,norbs) result(nx)
    type(Basis_t),intent(in)     :: bas
    logical, intent(in)          :: of
    real(dp), intent(in)         :: vec(:,:)
    integer, intent(in)          :: norbs
    integer                      :: i,j
    real(dp)                     :: psi(size(bas%f(:,1)),norbs), nx(size(bas%f(:,1))),dx
    psi(:,:) = 0.0_dp
    dx = bas%x(2)-bas%x(1)
    do i=1, norbs
        do j=1, size(bas%f(1,:))
            psi(:,i) = psi(:,i) + vec(j,i) * bas%f(:,j)
        end do
    end do
    nx(:) = 0.0_dp
    if(of .eqv. .true.) then
        nx(:) = bas%nelectrons* psi(:,1)**2.0_dp
    else
        do i=1, norbs
            nx(:) = nx(:) + 2.0_dp*psi(:,i)**2.0_dp
        end do
    end if
end function
function expand_psi(bas,vec,norbs) result(psi)
    type(Basis_t),intent(in) :: bas
    real(dp),intent(in)      :: vec(:,:)
    integer                  :: i,j,norbs
    real(dp)                 :: psi(size(bas%f(:,1)),norbs)

    psi(:,:) = 0.0_dp
    do i=1, norbs
        do j=1, size(bas%f(1,:))
            psi(:,i) = psi(:,i) + vec(j,i) * bas%f(:,j)
        end do
    end do

end function
function expand_psi_ofdft(bas,vec) result(psi)
    type(basis_t),intent(in) :: bas
    real(dp),intent(in)      :: vec(:)
    integer                  :: i
    real(dp)                 :: psi(size(bas%x))

    psi(:) = 0.0_dp
    do i=1, bas%nfuns
        psi(:) = psi(:) + vec(i)*bas%f(:,i)
    end do
    

end function
function derivative(bas,nx) result(dnx)
    type(basis_t),intent(in) :: bas
    real(dp), intent(in)     :: nx(:)
    real(dp)                 :: dnx(size(bas%x)), dermat(size(bas%x),size(bas%x)),dx
    integer                  :: i,j
    dx = bas%x(2)-bas%x(1)
    dnx(:) = 0.0_dp
    do i=1, 3
        dnx(i) = -11.0_dp*nx(i) + 18.0_dp * nx(i+1) - 9.0_dp * nx(i+2) + 2.0_dp * nx(i+3)
        dnx(i) = dnx(i)/(6.0_dp*dx)
    end do
    do i=4, size(bas%x)-3
        dnx(i) = nx(i-3)-45.0_dp*nx(i-1)-60.0_dp*nx(i)+135.0_dp*nx(i+1)-36.0_dp*nx(i+2)+5.0_dp*nx(i+3)
        dnx(i) = dnx(i)/(120.0_dp * dx)
    end do
    do i=size(bas%x)-2,size(bas%x)
        dnx(i) = -2.0_dp * nx(i-3) + 9.0_dp * nx(i-2) - 18.0_dp * nx(i-1) + 11.0_dp*nx(i)
        dnx(i) = dnx(i)/(6.0_dp*dx)
    end do
end function
end module basis
