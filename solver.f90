module solver
use types, only: dp
use basis, only: Basis_t, expand_rho
use hamilton
implicit none
type, public :: eigres_t
    integer :: norbs
    real(dp), allocatable :: energies(:), vec(:,:)
end type


type, public :: dft_t
    real(dp), allocatable :: energies(:), psi(:,:), nx(:)
end type
interface eigres_t
    module procedure :: create_eigres
end interface

interface dft_t
    module procedure :: create_dft
end interface
contains 

function create_eigres(norbs,nfuns) result(newres)
    type(eigres_t)       :: newres 
    integer, intent(in)  :: norbs, nfuns
    allocate(newres%energies(norbs),newres%vec(nfuns,norbs))
    newres%norbs = norbs
    newres%energies = 0.0_dp
    newres%vec = 0.0_dp
end function
function create_dft(bas,norbs,gridpoints, energies, vec, nx) result(df)
    type(basis_t),intent(in):: bas
    type(dft_t)             :: df
    integer                 :: norbs, gridpoints
    real(dp)                :: vec(:,:), energies(:),nx(:)
    real(dp)                :: psi(gridpoints,norbs)
    integer                 :: i,j

    allocate(df%energies(norbs),df%psi(gridpoints,norbs),df%nx(gridpoints))
    df%energies = energies
    psi(:,:) = 0.0_dp
    do i=1, norbs
        do j=1, size(bas%f(1,:))
            psi(:,i) = psi(:,i) + vec(j,i) * bas%f(:,j)
        end do
    end do
    df%psi = psi
    df%nx = nx

end function
function solve_eigenproblem(H,b,norbs) result(sol)
    type(eigres_t)            :: sol
    type(basis_t),intent(in)  :: b
    integer, intent(in)       :: norbs
    real(dp)                  :: H(:,:), A(b%nfuns,b%nfuns), W(b%nfuns), VU, VL
    character                 :: JOBZ, RANG,UPLO
    integer                   :: N, LDA, IL, IU,M,LDZ,LWORK, LIWORK, INFO
    real(dp)                  :: ABSTOL
    real(dp),allocatable      :: Z(:,:),works(:)
    integer, allocatable      :: ISUPPZ(:), IWORK(:)


    sol = eigres_t(norbs, b%nfuns)

    A(:,:) =H(:,:)
    JOBZ = 'V'
    RANG = 'I'
    UPLO = 'U'
    N = size(H(:,1))
    LDA = N
    VL = 0.0_dp
    VU = 0.0_dp
    IL = 1
    IU = norbs
    ABSTOL = -1.0_dp
    LDZ = N
    allocate(Z(LDZ,N))
    allocate(ISUPPZ(N))
    LWORK = -1
    LIWORK = -1
    allocate(works(1),IWORK(1))

    call DSYEVR(JOBZ,RANG,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,ISUPPZ,WORKS,LWORK,IWORK,LIWORK,INFO)
    LWORK = min(26*N, int(works(1)))
    LIWORK = min(26*N,IWORK(1))
    deallocate(works,IWORK)
    allocate(works(LWORK))
    allocate(IWORK(LIWORK))
    call DSYEVR(JOBZ,RANG,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,ISUPPZ,WORKS,LWORK,IWORK,LIWORK,INFO)
    
    sol%energies(1:norbs) = W(1:norbs)
    sol%vec(:,1:norbs)    = Z(:,1:norbs)

end function
function diis(bas,nx_diff,nx_diis) result(frho)
    type(Basis_t),intent(in):: bas
    real(dp), intent(in)    :: nx_diis(:,:),nx_diff(:,:)
    real(dp)                :: frho(size(bas%x)), B(size(nx_diis(:,1))+1,size(nx_diis(:,1))+1),dx, RHS(size(nx_diis(:,1))+1)
    integer                 :: i,j, N, NRHS, LDA, IPIV(size(B(:,1)) ), LDB, INFO
    dx = bas%x(2)-bas%x(1)
    do i=1,size(nx_diis(:,1))
        do j=1, size(nx_diis(:,1))
            B(i,j) = sum(nx_diff(i,:)*nx_diff(j,:)*dx)
        end do

    end do
    B(:,size(nx_diis(:,1))+1) = -1.0_dp
    B(size(nx_diis(:,1))+1,:) = -1.0_dp
    B(size(nx_diis(:,1))+1,size(nx_diis(:,1))+1) = 0.0_dp
    RHS = 0.0_dp
    RHS(size(RHS)) = -1.0_dp
    NRHS=1
    LDA = size(B(:,1)) 
    N=LDA
    LDB = size(B(:,1))
    call dgesv(N,NRHS,B,LDA,IPIV,RHS,LDB,INFO)
    if(INFO /= 0) then
        print *, "oh no, DIIS failed at DGESV"
    end if
    frho = 0.0_dp
    do i=1, size(nx_diis(:,1))
        frho(:) = frho(:)+RHS(i)*nx_diis(i,:)
    end do
    
end function
function solve_dft(bas,initrho,beta,of,norbs,diisflag,ndiis) result(res)
    type(basis_t),intent(in)    :: bas
    logical,intent(in)          :: of, diisflag
    real(dp), intent(in)        :: initrho(:)
    type(eigres_t)              :: sol
    type(dft_t)                 :: res
    integer, intent(in)         :: norbs,ndiis
    real(dp)                    :: preve, alpha, H(size(bas%f(1,:)),size(bas%f(1,:) )), nx(size(bas%f(:,1) )),beta(:)
    real(dp)                    :: vec(size(bas%f(1,:)),norbs), ediff, nx_old(size(bas%f(:,1) )),st,en,dx
    real(dp)                    :: nx_diis(ndiis,size(bas%x)), frho(size(bas%x)),nx_diff(ndiis,size(bas%x)), mu
    integer                     :: i, maxiter,j,ierr
    maxiter = 3000
    dx = bas%x(2)-bas%x(1)
    preve = -10000.0_dp
    H = 0.0_dp
    alpha = 0.01_dp
    ediff = 0.0_dp
    sol = eigres_t(norbs,size(bas%f(1,:)))
    sol%vec(:,:) = 0.0_dp
    nx_old = 0.0_dp
    nx(:) = initrho(:)
    mu = 15.0_dp
    print *, diisflag
    if(diisflag .eqv. .false.) print *, "noh"
    do i=1,maxiter
        H = 0.0_dp
        call get_hartree(H,bas,nx)
        call get_exchange(H,bas,nx)
        call get_zero_energy(H,bas)
        !call ext_kronig(H,10,100.0_dp,1.0_dp,1.0_dp/6.0_dp,bas%L)
        if(of .eqv. .true.) then
            call get_pauli_appro(H,bas,nx,mu)
            !call get_funny_fullpot(H, bas, nx)
            !call get_thomas_fermi(H,bas,nx)
        end if
        sol = solve_eigenproblem(H,bas,norbs)
        ediff = abs(sol%energies(1) - preve)
        preve = sol%energies(1)
        mu = sol%energies(1)
        !print *, i
        print *, sol%energies(1), sum(dx*(nx(:)-nx_old(:))**2),i
        if(sum(dx*(nx(:)-nx_old(:))**2) < 1e-8_dp) then
            print *, "Calculation has converged"
            !nx = expand_rho(bas,sol%vec,of,norbs)
            res = dft_t(bas,norbs, size(bas%f(:,1)) ,sol%energies,sol%vec,nx )
            frho = diis(bas,nx_diff,nx_diis)
            open(20,file="rhodiis.dt",STATUS='OLD', IOSTAT=IERR, ERR=100)
            100 print *, "hm"

            do j=1,size(bas%x)
                write(20,*) bas%x(j), frho(j)
            end do
            close(20)
            return
        end if
        nx_old(:) = nx(:)
        if((i < ndiis+1) .and. (diisflag .eqv. .true.)) then
            print *, "LINEAR MIXING"
            nx(:) =(1-alpha)*nx_old + (alpha) * expand_rho(bas,sol%vec,of,norbs)
            nx_diff(i,:) = nx(:)-nx_old(:)
            nx_diis(i,:) = nx(:)
        elseif((i>= ndiis+1) .and. (diisflag .eqv. .true.) .and. (mod(i,ndiis) /=  0)) then
            print *, "LINEAR MIXING"
            nx(:) =(1-alpha)*nx_old + (alpha) * expand_rho(bas,sol%vec,of,norbs)
            do j=1, ndiis-1
                nx_diff(j,:) = nx_diff(j+1,:)
                nx_diis(j,:) = nx_diis(j+1,:)
            end do
            nx_diff(ndiis,:) = nx(:)-nx_old(:)
            nx_diis(ndiis,:) = nx(:)
        elseif((i>= ndiis+1) .and. (diisflag .eqv. .true.) .and. (mod(i,ndiis) ==  0)) then
            print *, "DIIS"
            !nx(:) =(1.0_dp-alpha)*nx_old(:)+ alpha*diis(bas, nx_diff,nx_diis)
            nx(:) = diis(bas, nx_diff,nx_diis)
            do j=1, ndiis-1
                nx_diff(j,:) = nx_diff(j+1,:)
                nx_diis(j,:) = nx_diis(j+1,:)
            end do
            nx_diff(ndiis,:) = nx(:)-nx_old(:)
            nx_diis(ndiis,:) = nx(:)
        elseif(diisflag .eqv. .false.) then
            print *, "LINEAR MIXING"
            nx(:) =(1-alpha)*nx_old + (alpha) * expand_rho(bas,sol%vec,of,norbs)
        end if
        print *, "------------------------"

    end do
    
end function
function solve_ofdft(bas,initrho, norbs,nbp, diisflag, ndiis,lindiis) result(res)
    type(basis_t),intent(in)    :: bas
    real(dp), intent(in)        :: initrho(:,:)
    integer, intent(in)         :: norbs,nbp,ndiis,lindiis
    logical, intent(in)         :: diisflag
    type(eigres_t)              :: sol(nbp)
    type(dft_t)                 :: res
    real(dp)                    :: preve, alpha, H(nbp,size(bas%f(1,:)),size(bas%f(1,:) )), nx(nbp,size(bas%f(:,1) ))
    real(dp)                    :: vec(size(bas%f(1,:)),norbs), ediff, nx_old(nbp,size(bas%f(:,1) )),st,en,dx,db
    real(dp)                    :: nx_diis(nbp,ndiis,size(bas%x)), frho(nbp,size(bas%x)),nx_diff(nbp,ndiis,size(bas%x))
    real(dp)                    :: pauli(nbp,size(bas%f(1,:)),size(bas%f(1,:) ))
    integer                     :: i, maxiter,j,k
    maxiter = 2000
    db = 0.01_dp
    preve = -10000.0_dp
    H = 0.0_dp
    dx = bas%x(2)-bas%x(1)
    alpha = 0.05_dp
    ediff = 0.0_dp
    do i=1, nbp
        sol(i) = eigres_t(norbs,size(bas%f(1,:)))
        sol(i)%vec(:,:) = 0.0_dp
    end do
    nx_old = 0.0_dp
    nx(:,:) = initrho(:,:)
    en = 2.86_dp
    
    do i=1, maxiter
        H(:,:,:) = 0.0_dp
        do j=1,nbp
            call get_zero_energy(H(j,:,:),bas)
            call get_hartree(H(j,:,:),bas,nx(1,:))
            call get_exchange(H(j,:,:),bas,nx(1,:))
        end do
        !call get_pauli_beta_exact(H(:,:,:),bas,(/ 4.69587_dp, 5.22630_dp /),db,nx(:,:))
        if(mod(i,1) == 0 .or. i ==1) then
            call get_pauli_beta_rho(pauli,bas,en,db,nx)
            H(:,:,:) = H(:,:,:) + pauli(:,:,:)
        else
            H(:,:,:) = H(:,:,:) + pauli(:,:,:)
        end if
        do j=1, nbp
            sol(j) = solve_eigenproblem(H(j,:,:),bas,1)
        end do
        
        print *, sol(1)%energies(1),sum(dx*(nx(1,:)-nx_old(1,:))**2),i,"hmm"
        !en = sol(1)%energies(1)
        en = sol(1)%energies(1)
        if(sum(dx*(nx(1,:)-nx_old(1,:))**2) < 1e-8_dp) then
            print *, "Calculation has converged"
            nx(1,:) = expand_rho(bas,sol(1)%vec,.true.,norbs)
            res = dft_t(bas,norbs, size(bas%f(:,1)) ,sol(1)%energies,sol(1)%vec,nx(1,:) )
            return
        end if
        !DIIS
        nx_old(:,:) = nx(:,:)
        do k=1, nbp
            if(i < ndiis+1 .and. diisflag .eqv. .true.) then
                if (k==1) print *, "LINEAR MIXING"
                
                nx(k,:) =(1-alpha)*nx_old(k,:) + (alpha) * expand_rho(bas,sol(k)%vec,.true.,norbs)
                nx_diff(k,i,:) = nx(k,:)-nx_old(k,:)
                nx_diis(k,i,:) = nx(k,:)
            elseif(i>= ndiis+1 .and. diisflag .eqv. .true. .and. mod(i,lindiis) /=  0) then
                if (k==1) print *, "LINEAR MIXING"
                nx(k,:) =(1-alpha)*nx_old(k,:) + (alpha) * expand_rho(bas,sol(k)%vec,.true.,norbs)
                do j=1, ndiis-1
                    nx_diff(k,j,:) = nx_diff(k,j+1,:)
                    nx_diis(k,j,:) = nx_diis(k,j+1,:)
                end do
                nx_diff(k,ndiis,:) = nx(k,:)-nx_old(k,:)
                nx_diis(k,ndiis,:) = nx(k,:)
            elseif(i>= ndiis+1 .and. diisflag .eqv. .true. .and. mod(i,lindiis) == 0) then
                if(k==1) print *, "DIIS"
                !nx(k,:) =(1.0_dp-alpha)*nx_old(k,:)+ alpha*diis(bas, nx_diff(k,:,:),nx_diis(k,:,:))
                nx(k,:) = diis(bas, nx_diff(k,:,:),nx_diis(k,:,:))
                do j=1, ndiis-1
                    nx_diff(k,j,:) = nx_diff(k,j+1,:)
                    nx_diis(k,j,:) = nx_diis(k,j+1,:)
                end do
                nx_diff(k,ndiis,:) = nx(k,:)-nx_old(k,:)
                nx_diis(k,ndiis,:) = nx(k,:)
            else
                nx(k,:) =(1-alpha)*nx_old(k,:) + (alpha) * expand_rho(bas,sol(k)%vec,.true.,norbs)
            end if
        end do
        !LINEAR MIXING
        print *, "-------------"

    end do
end function
end module solver