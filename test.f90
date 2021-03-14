program testi
use types, only:dp
use basis
use hamilton
use solver
use optimize, only: deq_optimize,optimizer
use energies
!use optimizer
implicit none

abstract interface
    function func(bas,alpha,sn,vec,mu) result(k)
        import Basis_t
        import dp
        type(basis_t) :: bas
        real(kind(0.d0))  :: sn(:),vec(:),k,alpha,beta, kindens(size(bas%x)), deldens(size(bas%x)),mu
        real(dp)          :: psisn(size(bas%x)),pi
        integer           :: i,j 
    end function
    function gs(bas,vec,mu) result(grad)
        import dp
        import Basis_t
        type(Basis_t)       :: bas
        real(dp)            :: vec(:), psi(size(bas%x)),grad(bas%nfuns),alpha,beta,pi,mu
        integer             :: i
    end function
end interface

integer :: numfail, norbs, nfuns,nbar,i
real(dp):: ergs, H(100,100), L, xi(10), b, pi, nx(1000), x(2000), nx2(2000)
procedure(func),pointer :: f_ptr => null()
procedure(gs),pointer   :: g_ptr => null()
type(eigres_t) :: sol2
type(basis_t) :: testb
type(dft_t)     :: res, res2
nfuns = 100
numfail = 0
norbs = 10
nbar = 10
pi = 4.0_dp*atan(1.0_dp)
b = 1.0_dp/6.0_dp
L = 10.0_dp
testb = basis_t(2*norbs,nfuns,1000,L)
if(abs(sum(testb%f(:,5)**2.0_dp)*(testb%x(2)-testb%x(1))-1.0_dp) < 1e-2_dp ) then
    print *, "basis functions OK"
else
    numfail = numfail + 1 
    print *, "failure at basis functions"
end if
if(abs(sum(derivative(testb,testb%x**2.0_dp))-sum(2.0_dp*testb%x(:)))*(testb%x(2)-testb%x(1))<1e-2_dp) then
    print *, "derivatives ok"
else
    numfail = numfail +1
    print *, "derivatives failed"
end if
H = 0.0_dp
call get_zero_energy(H, testb)

sol2 = solve_eigenproblem(H,testb,norbs)
if(abs(sol2%energies(2)-0.197_dp)<1e-2 .and. abs(sol2%energies(3)-0.444_dp)<1e-2) then
    print *, "eigensolver first check ok"
else
    numfail =numfail + 1
    print *, "Eigensolver failed basic test 1"
end if
H = 0.0_dp
call get_zero_energy(H,testb)
call ext_kronig(H,nbar,100.0_dp,1.0_dp,b,L)
sol2 = solve_eigenproblem(H, testb, norbs)
if(abs(sol2%energies(9)-sol2%energies(10)) > 15.0_dp) then
    print *, "eigensolver second check OK; found band gap"
else
    numfail = numfail + 1
    print *, "Eigensolver failed KP band gap test"
end if
do i=1, norbs
    nx(:) = nx(:)+2.0_dp * 2.0_dp/testb%L * sin(i*pi/testb%L * testb%x(:))**2.0_dp
end do
res = solve_dft(testb,nx,nx,.false.,norbs,.false.,3)
res2 = solve_dft(testb,nx,nx,.false.,norbs,.true.,3)
if(abs(sum(res%nx(:)*(testb%x(2)-testb%x(1)))-sum(res2%nx(:)*(testb%x(2)-testb%x(1))))<1e-1_dp) then
    print *, "DIIS check ok"
else
    numfail = numfail + 1
    print *, "KS-DFT iteration failed DIIS test"
end if
f_ptr => get_thomas_fermi_total_energy
g_ptr => get_thomas_fermi_gradient
norbs = 5
L = 8.0_dp
testb = basis_t(10,60,2000,L)
call optimizer(testb,f_ptr,g_ptr)
open(12,file="rho-minim.dt")
do i=1,2000
    read(12,*) x(i), nx2(i)
end do
close(12)
do i=200, 1800
    if(abs(nx2(i)) < 1e-3_dp ) then
        print *, "failure at energy minimization"
        numfail = numfail +1
        exit
    end if
end do
end program
