
program hello
    implicit none
    include "mpif.h"
    integer,parameter::n_dim=20
    integer,parameter::n_eqb=15000
    integer,parameter::n_save=15000
    integer,parameter::rel=20
    real,parameter::pi=4.0*atan(1.0)
    real,parameter::J2=0.70
    integer::ranks,proc,ierr
    integer::i,j,k,l,k1,k2,k3,k4,a,b
    integer::nbr(8,n_dim**2),sb_lt(4,(n_dim**2)/4),nbr_full(8,2)
    real::T,T_in,mc,ovr,chiral,chi4,chi2
    real::Q(n_dim**2),U_all(rel)
    real,allocatable::U(:)
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,ranks,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,proc,ierr)
    allocate(U(rel/proc))
    do i=1,n_dim
        do j=1,n_dim
            nbr_full(1,1)=modulo(i-2,n_dim)+1
            nbr_full(1,2)=j

            nbr_full(2,1)=modulo(i-2,n_dim)+1
            nbr_full(2,2)=modulo(j-2,n_dim)+1

            nbr_full(3,1)=i
            nbr_full(3,2)=modulo(j-2,n_dim)+1

            nbr_full(4,1)=modulo(i,n_dim)+1
            nbr_full(4,2)=modulo(j-2,n_dim)+1

            nbr_full(5,1)=modulo(i,n_dim)+1
            nbr_full(5,2)=j

            nbr_full(6,1)=modulo(i,n_dim)+1
            nbr_full(6,2)=modulo(j,n_dim)+1

            nbr_full(7,1)=i
            nbr_full(7,2)=modulo(j,n_dim)+1

            nbr_full(8,1)=modulo(i-2,n_dim)+1
            nbr_full(8,2)=modulo(j,n_dim)+1

            do k=1,8
                nbr(k,j+n_dim*(i-1))=nbr_full(k,2)+n_dim*(nbr_full(k,1)-1)
            end do
        end do
    end do
    k1=1
    k2=1
    k3=1
    k4=1
    do i=1,n_dim
        do j=1,n_dim
            a=modulo(i,2)
            b=modulo(j,2)
            if(a/=0 .and. b/=0) then
                sb_lt(1,k1)=j+n_dim*(i-1)
                k1=k1+1
            elseif(a/=0 .and. b==0) then
                sb_lt(2,k2)=j+n_dim*(i-1)
                k2=k2+1
            elseif(a==0 .and. b/=0) then
                sb_lt(3,k3)=j+n_dim*(i-1)
                k3=k3+1
            else
                sb_lt(4,k4)=j+n_dim*(i-1)
                k4=k4+1
            end if
        end do
    end do
    do k=1,10
    do l=1,rel/proc
        call random_number(Q)
        Q=Q*2*pi-pi
        T=0.559+0.001*k
        k4=0.0
        k2=0.0
        chi2=0.0
        chi4=0.0
        T_in=1/T
        do i=1,n_eqb/2
            do j=1,n_dim**2
                Q(j)=mc(Q,j,nbr,T_in,J2,n_dim,pi)
            end do
            do j=1,n_dim**2
                Q(j)=ovr(Q,j,nbr,T_in,J2,n_dim,pi)
            end do
        end do

        do i=1,n_save/2
            do j=1,n_dim**2
                Q(j)=mc(Q,j,nbr,T_in,J2,n_dim,pi)
            end do
            do j=1,n_dim**2
                Q(j)=ovr(Q,j,nbr,T_in,J2,n_dim,pi)
            end do
            if(modulo(i,5)==0) then
                chi4=chi4+chiral(Q,n_dim,nbr)**4.0
                chi2=chi2+chiral(Q,n_dim,nbr)**2.0
            end if
        end do
        U(l)=1-0.10*n_save*chi4/(3*chi2**2)
    end do
    call MPI_Gather(U,rel/proc,MPI_REAL,U_all,rel/proc,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    if(ranks==0) then
        print*,sum(U_all)/rel,"      ",sqrt(sum(U_all**2)/rel-(sum(U_all)/rel)**2)
    end if
    end do
    deallocate(U)
    call MPI_Finalize(ierr)
end program

function mc(Q,ind,nbr,T_in,J2,n_dim,pi)
    implicit none
    integer::n_dim,ind
    integer::nbr(8,n_dim**2)
    real::T_in,J2,pi,wf,wi,Ei,Ef,r,mc
    real::Q(n_dim**2),p(8)
    wi=Q(ind)
    call random_number(wf)
    call random_number(r)
    wf=wf*2*pi-pi
    p=Q(nbr(1:8,ind))
    Ei=-1*(cos(wi-p(1))+cos(wi-p(3))+cos(wi-p(5))+cos(wi-p(7)))+J2*(cos(wi-p(2))+cos(wi-p(4))+cos(wi-p(6))+cos(wi-p(8)))
    Ef=-1*(cos(wf-p(1))+cos(wf-p(3))+cos(wf-p(5))+cos(wf-p(7)))+J2*(cos(wf-p(2))+cos(wf-p(4))+cos(wf-p(6))+cos(wf-p(8)))
    if(exp(T_in*(Ei-Ef))>r) then
        mc=wf
    else
        mc=wi
    end if
end function

function ovr(Q,ind,nbr,T_in,J2,n_dim,pi)
    implicit none
    integer::n_dim,ind
    integer::nbr(8,n_dim**2)
    real::T_in,J2,pi,wi,ovr,hloc_x,hloc_y
    real::Q(n_dim**2),p(8)
    wi=Q(ind)
    p=Q(nbr(1:8,ind))
    hloc_x=-cos(p(1))-cos(p(3))-cos(p(5))-cos(p(5))+J2*(cos(p(2))+cos(p(4))+cos(p(6))+cos(p(8)))
    hloc_y=-sin(p(1))-sin(p(3))-sin(p(5))-sin(p(5))+J2*(sin(p(2))+sin(p(4))+sin(p(6))+sin(p(8)))
    ovr=2*atan2(hloc_y,hloc_x)-wi
end function

function chiral(Q,n_dim,nbr)
    implicit none
    integer::n_dim,i
    integer::nbr(8,n_dim**2)
    real::chiral,qi,qj,qk,ql
    real::Q(n_dim**2)
    chiral=0.0
    do i=1,n_dim**2
        qi=Q(i)
        qj=Q(nbr(7,i))
        qk=Q(nbr(8,i))
        ql=Q(nbr(1,i))
        chiral=chiral+0.25*(cos(qi-qj)-cos(qi-ql)-cos(qk-qj)+cos(qk-ql))
    end do
    chiral=chiral/n_dim**2
end function
