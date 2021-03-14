module optimize
use types, only: dp
use basis, only: basis_t, expand_psi_ofdft
implicit none
public deq_optimize, optimizer
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

contains
function gss(nb,f,a,b,sn,vec,mu) result(gold)
    !Golden Section line search copied from wikipedia
    type(Basis_t)    :: nb
    procedure(func)  :: f
    real(dp)         :: a,b,sn(:),c,d,gr,tol,vec(:),gold,mu
    tol = 1e-35_dp
    gr = ( sqrt(5.0_dp) +1.0_dp)/2.0_dp
    c = b - (b-a) / gr
    d = a + (b-a) / gr
    do while(abs(c-d) > tol)
        if(f(nb,c,sn,vec,mu) < f(nb,d,sn,vec,mu) ) then
            b = d
        else 
            a = c
        end if

        c = b - (b - a) / gr
        d = a + (b - a) / gr
    end do

    gold = (b+a)/2.0_dp    
end function
function beta(g,gprev) result(b)
    !Beta parameter for conjugate gradient according to Polak-Ribiere formula
    real(dp) :: b, g(:),gprev(:)

    b = sum(-g(:)*(-g(:)+gprev(:) ))/(sum(gprev(:)*gprev(:)))
    
end function
subroutine optimizer(bas,defunc,gradfunc)
    !Optimizes electron density. Mu, the Lagrange multiplier ensuring the number of electrons
    !stays constant, is changed with the bisection method and the density is then minimized.
    !Once the correct number of electrons is found, the calculation is converged.
    type(Basis_t)     :: bas
    procedure(func)   :: defunc
    procedure(gs)     :: gradfunc
    real(dp)          :: mu, a,b,c, optres
    a = 0.0_dp
    b = 50.0_dp
    c = 25.0_dp
    optres = 10.0_dp
    !goto 10
    do while(abs(optres) > 1e-2_dp)
        c = (a+b)/2.0_dp
        optres = deq_optimize(bas,defunc,gradfunc,c)
        if(optres < 0.0_dp) then
            a=c
        else
            b=c
        end if
        print *, a,b
    end do
    print *,"mu now", c
    !10 continue
    !optres = deq_optimize(bas,defunc,gradfunc,a)
    
end subroutine
function deq_optimize(bas,defunc,gradfunc,mu) result(nelectron)
    !Conjugate-gradient optimization. Returns the number of electrons after minimization.
    type(Basis_t)     :: bas
    real(dp)          :: grad(bas%nfuns),gradprev(bas%nfuns),b,stp,etp
    real(dp)          :: sn(bas%nfuns),snprev(bas%nfuns), curr(bas%nfuns),prev(bas%nfuns),gamma,mu, psi(size(bas%x))
    real(dp)          :: rho(size(bas%x)), prevenergy,currenergy,nelectron,pi
    integer           :: i,j
    procedure(func)   :: defunc
    procedure(gs)     :: gradfunc
    !mu = 8.7_dp
    pi = 4.0_dp*atan(1.0_dp)
    stp = 0.0_dp
    etp = 10.0_dp
    prev(:) = 0.0_dp
    prev(:) = 0.05_dp
    !pi = defunc(bas,0.0_dp,prev,prev,mu)
    psi(:) = expand_psi_ofdft(bas,prev)
    !print *, (bas%x(2)-bas%x(1))*sum(psi(:)**2.0_dp)
    curr(:) = 0.0_dp
    sn(:) = 0.0_dp
    grad(:) = gradfunc(bas,prev,mu)
    snprev(:) = -grad(:)
    !print *, grad
    gamma = gss(bas, defunc, stp,etp,-grad(:),prev(:),mu)
    print *, gamma
    !gamma = 0.0001_dp
    stp = 0.0_dp
    etp = 10.0_dp
    curr(:) = prev(:) - gamma*grad(:)
  
    prevenergy = defunc(bas,stp,sn(:),curr(:),mu)
    do i=1, 100000
        psi(:) = expand_psi_ofdft(bas,curr)
        !print *, sum(psi(:)),"lol"
        !print *, "mi"
        gradprev(:) = grad(:)
        grad(:)     = gradfunc(bas,curr,mu)
        b           = beta(grad,gradprev)
        b           = max(0.0_dp,b)
        !b           = 0.0_dp
        sn(:)       = -grad(:) + b*snprev(:)
        snprev(:)   = sn(:)
        gamma       = gss(bas, defunc, stp,etp,sn(:),curr(:),mu)
        print *, gamma
        !These are changed inside GSS, so best to reset them
        stp         = 0.0_dp
        etp         = 10.0_dp
        !
        curr(:)     = prev(:) + gamma*sn(:)

        !bas%varcoeff(:) = curr(:)
        prev(:)     = curr(:)
        
        !if(defunc(basis,stp,sn(:),r) < 1e-5_dp) then
        print *, i, sum(grad*grad)
        if(sum(grad(:)*grad(:)) < 1e-8_dp .or. &
          (abs(defunc(bas,stp,sn(:),curr(:),mu)-prevenergy )<1e-8_dp .and.sum(grad(:)*grad(:)) < 2.0_dp ) ) then
            !print *, "stopped",i
            
            
            
            open(12, file="rho-minim.dt")
            psi = expand_psi_ofdft(bas,curr(:))
            rho = psi**2.0_dp
            nelectron= sum(rho)*(bas%x(2)-bas%x(1))-bas%nelectrons
            print *, "electron difference",nelectron
            do j=1, size(bas%x)
                write(12, *) bas%x(j), rho(j)
            end do
            !print *, curr
            close(12)
            exit
        end if
        prevenergy = defunc(bas,stp,sn(:),curr(:),mu)
    end do
    !print *, curr(:)
end function

end module