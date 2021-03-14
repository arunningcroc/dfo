program mainn
use types, only:dp
use hamilton
use solver, only: solve_eigenproblem, eigres_t,dft_t,solve_dft,solve_ofdft
use basis
use energies
implicit none
integer :: npeak
type(eigres_t) :: sol(7), sol2
type(basis_t) :: testb
type(dft_t)     :: res, res1(5)
real(dp):: H(1,100,100),a,b,L,V0,nx(2000),hart(5,2000), pi, psi(2000,5),nx3(2000),db,nx2(2000)
real(dp):: ergs(5)
integer :: i,j,k,norbs,nelectrons
norbs = 5
nelectrons=10

pi = 4.0_dp*atan(1.0_dp)
npeak = 10
a=1.0_dp
b=1.0_dp/6.0_dp
!cnt = 0.02_dp
L = 8.0_dp
V0 = 100.0_dp
H(:,:,:) = 0.0_dp
testb = basis_t(nelectrons,100,2000,L)
nx(:) = 0.0_dp

do i=1, norbs
    nx(:) = nx(:)+2.0_dp * 2.0_dp/testb%L * sin(i*pi/testb%L * testb%x(:))**2.0_dp
    psi(:,i) = sqrt(2.0_dp/testb%L) * sin(i*pi/testb%L * testb%x(:))
    ergs(i) = i**2.0_dp * pi**2.0_dp/(2.0_dp * L**2.0_dp)
end do


!res1(2) = solve_dft(testb,nx,hart(1,:),.false.,norbs,.true.,10)

!print *, sum((res1(1)%nx(:) -res1(2)%nx(:))**2*(testb%x(2)-testb%x(1)) )
!print *, res1(2)%energies
!call get_pauli_beta_rho(H,testb,sol2%energies(2),db,rhobeta)
!call get_pauli_exact(H(1,:,:),testb,res1(2)%psi,res1(2)%energies,res1(2)%nx)
call get_pauli_exact(H(1,:,:),testb,psi,ergs,nx)
call get_exchange(H(1,:,:),testb,nx)
call get_hartree(H(1,:,:),testb,nx)
call get_zero_energy(H(1,:,:),testb)
sol2 = eigres_t(1,size(testb%f(1,:)))
sol2 = solve_eigenproblem(H(1,:,:),testb,1)
print *, "hey"
nx2 = expand_rho(testb,sol2%vec,.true.,1)
H = 0.0_dp
call get_pauli_exact_symmetry(H(1,:,:),testb,psi,ergs,nx)
call get_exchange(H(1,:,:),testb,nx)
call get_hartree(H(1,:,:),testb,nx)
call get_zero_energy(H(1,:,:),testb)
sol2 = eigres_t(1,size(testb%f(1,:)))
sol2 = solve_eigenproblem(H(1,:,:),testb,1)
nx3 = expand_rho(testb,sol2%vec,.true.,1)
!call get_pauli_continuous_appro(H(1,:,:),testb,res1(2)%nx)
!print*, "Kinetic", get_kinetic_energy(testb,res1(2)%psi)
!print*, "Exchange", get_exchange_energy(testb,res1(2)%nx)
!print*, "Interaction", get_delta_energy(testb,res1(2)%nx)
!print*, "Total", get_kinetic_energy(testb,res1(2)%psi)+get_exchange_energy(testb,res1(2)%nx)+get_delta_energy(testb,res1(2)%nx)
!call get_pauli_appro(H(1,:,:),testb,res1(1)%nx)
!print *, res1(2)%energies
open(12, file="rho.dt")
do i=1, size(testb%x)
    write(12, *) testb%x(i),nx(i),nx2(i),nx3(i)
end do
close(12)


end program