program minimm
use types, only: dp
use basis
use optimize, only: deq_optimize,optimizer
use energies
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

type(basis_t)           :: bas
procedure(func),pointer :: f_ptr => null()
procedure(gs),pointer   :: g_ptr => null()
real(dp)                :: L
integer                 :: nel

nel = 10
L = 8.0_dp
!f_ptr => get_vw_total_energy
!g_ptr => get_vw_gradient
!f_ptr =>get_nagy_total_energy
!g_ptr => get_nagy_gradient
f_ptr => get_vw_total_energy
g_ptr => get_vw_gradient
!f_ptr => get_tf_pot
!g_ptr => get_tf_pot_grad
bas = basis_t(nel,60,2000,L)
!f_ptr => get_thomas_fermi_total_energy
!g_ptr => get_thomas_fermi_gradient
!call deq_optimize(bas,f_ptr,g_ptr)
call optimizer(bas,f_ptr,g_ptr)
!L = get_nagy_total_energy(bas,0.0_dp,0.0_dp)
end program