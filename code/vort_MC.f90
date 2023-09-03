

program calling_func
   implicit none

   interface nbr_vec
        function nbr_vec(n_dim)
            integer ::n_dim
            integer,dimension(8,n_dim**2) ::nbr_vec
        end function nbr_vec
   end interface

   interface sub_latt
    function sub_latt(n_dim)
        integer ::n_din,sub_latt(4,n_dim**2/4)
    end function sub_latt
   end interface


   integer,parameter::rel=10
   integer,parameter::n_dim=20
   integer,parameter::n_save=50000
   integer,parameter::n_eqb=50000
   real,parameter::pi=4*atan(1.0)
   real,parameter::J2=0.7
   integer::i,j,k,l
   real :: T_arr(100)=[(i,i=1,100,1)]*0.01+0.2
   real :: E,avK2,avK4,chi,chiral,t_start,t_end
   integer::nbr_all(8,n_dim**2),latt(4,n_dim**2/4),vortex,vor
   real:: Q(n_dim**2),U(rel),V(rel)
   
   nbr_all=nbr_vec(n_dim)
   latt=sub_latt(n_dim)
   do k=10,70,5
    !table=HB(J2,n_lf,n_q,T_arr(K))
    do j=1,rel

        avK2=0
        avK4=0
        vor=0
        chi=0
        E=0
        call random_number(Q)
        Q=2*pi*Q-pi
        do i=1,n_eqb
            !call mc_sweep(Q,T_arr(k),J2,E,n_dim,nbr_all,pi,0.0,[(l,l=1,n_dim**2,1)])
            !call Heat_bath(Q,T_arr(k),J2,table,n_dim,nbr_all,n_q,n_lf,pi)
            call parallel_mc_sweep(Q,T_arr(k),J2,E,n_dim,nbr_all,pi,0.0,latt)



        end do
        do i=1,n_save
            !call mc_sweep(Q,T_arr(k),J2,E,n_dim,nbr_all,pi,0.0,[(l,l=1,n_dim**2,1)])
            call parallel_mc_sweep(Q,T_arr(k),J2,E,n_dim,nbr_all,pi,0.0,latt)
            !call Heat_bath(Q,T_arr(k),J2,table,n_dim,nbr_all,n_q,n_lf,pi)
            !avM2=avM2+sum(cos(Q))**2 +sum(sin(Q))**2
            !avMx=avMx+sum(cos(Q))
            !avMy=avMy+sum(sin(Q))
            chi=chiral(Q,nbr_all,n_dim)
            vor=vor+vortex(Q,nbr_all,latt,n_dim,pi)
            avk2=avk2+chi**2
            avk4=avk4+chi**4
            !avE=avE+E
        end do
        !prop(j)=(avM2-(avMx**2+avMy**2)/n_save)/((n_dim**2)*n_save*T_arr(k))
        U(j)=1-n_save*avK4/(3*(avK2**2))
        V(j)=vor/n_save

   end do
   write(*,*) T_arr(k),sum(V)/rel,sum(U)/(rel),sqrt((sum(U**2)-(sum(U)**2)/rel)/rel)
   end do



end program calling_func


    function nbr_vec(n_dim)
        implicit none
        integer::n_dim,i,j,k,l,nbr_vec(8,n_dim**2),nbr(n_dim,n_dim,8,2)

        do i=1,n_dim
            do j=1,n_dim
                nbr(i,j,1,1)=modulo(i-2,n_dim)+1
                nbr(i,j,1,2)=modulo(j-2,n_dim)+1

                nbr(i,j,2,1)=modulo(i-2,n_dim)+1
                nbr(i,j,2,2)=j

                nbr(i,j,3,1)=modulo(i-2,n_dim)+1
                nbr(i,j,3,2)=modulo(j,n_dim)+1

                nbr(i,j,4,1)=i
                nbr(i,j,4,2)=modulo(j,n_dim)+1

                nbr(i,j,5,1)=modulo(i,n_dim)+1
                nbr(i,j,5,2)=modulo(j,n_dim)+1

                nbr(i,j,6,1)=modulo(i,n_dim)+1
                nbr(i,j,6,2)=j

                nbr(i,j,7,1)=modulo(i,n_dim)+1
                nbr(i,j,7,2)=modulo(j-2,n_dim)+1

                nbr(i,j,8,1)=i
                nbr(i,j,8,2)=modulo(j-2,n_dim)+1
            end do
        end do

        do i=1,n_dim
            do j=1,n_dim
                k=i+n_dim*(j-1)
                do l=1,8
                    nbr_vec(l,k)=nbr(i,j,l,1)+n_dim*(nbr(i,j,l,2)-1)
                end do
            end do
        end do
    end function nbr_vec

function sub_latt(n_dim)
    implicit none
    integer::n_dim,i,j,k1,k2,k3,k4,mi,mj,sub_latt(4,(n_dim**2)/4)
    k1=1
    k2=1
    k3=1
    k4=1
    do i=1,n_dim
        do j=1,n_dim
            mi=modulo(i,2)
            mj=modulo(j,2)
            if(mi==0 .and. mj==0) then
                sub_latt(1,k1)=i+(j-1)*n_dim
                k1=k1+1
            elseif(mi/=0 .and. mj==0) then
                sub_latt(2,k2)=i+(j-1)*n_dim
                k2=k2+1
            elseif(mi==0 .and. mj/=0) then
                sub_latt(3,k3)=i+(j-1)*n_dim
                k3=k3+1
            else
                sub_latt(4,k4)=i+(j-1)*n_dim
                k4=k4+1
            end if
        end do
    end do
end function sub_latt


function delE1(i,nbr,Qmat,J2,w,n_dim,twist)
    implicit none
    real,dimension(n_dim**2) :: Qmat
    integer ::i, n_dim
    real ::q, qr, ql, qd, qu, Ei, Ef, w, delE1, J2, qn, qs, qw, qe,twist
    integer,dimension(8,n_dim**2) :: nbr
    q=Qmat(i)
    qr=Qmat(nbr(4,i))+twist
    ql=Qmat(nbr(8,i))-twist
    qu=Qmat(nbr(2,i))
    qd=Qmat(nbr(6,i))
    qs=Qmat(nbr(5,i))+twist
    qn=Qmat(nbr(1,i))-twist
    qw=Qmat(nbr(3,i))+twist
    qe=Qmat(nbr(7,i))-twist
    Ei=-1*(cos(q-qr)+cos(q-ql)+cos(q-qu)+cos(q-qd)) +J2*(cos(q-qn)+cos(q-qs)+cos(q-qe)+cos(q-qw))
    Ef=-1*(cos(w-qr)+cos(w-ql)+cos(w-qu)+cos(w-qd)) +J2*(cos(w-qn)+cos(w-qs)+cos(w-qe)+cos(w-qw))
    delE1=Ef-Ei
end function delE1


subroutine parallel_mc_sweep(Q,T,J2,E,n_dim,all_nbr,pi,twist,sub_latt)
    implicit none
    integer:: n_dim,i,j
    integer::all_nbr(8,n_dim**2),sub_latt(4,n_dim**2/4)
    real :: T,pi,E,J2,twist,dE,w,p,delE1
    real,dimension(n_dim**2) :: Q
    do i=1,4
        do j=1,n_dim**2/4
            call random_number(w)
            call random_number(p)
            w=2*pi*w-pi
            dE=delE1(sub_latt(i,j),all_nbr,Q,J2,w,n_dim,twist)
            if (exp(-1*dE/T)>p) then
                Q(sub_latt(i,j))=w
                E=E+dE
            end if
        end do
    end do
end subroutine parallel_mc_sweep



function chiral(Q,adj,n_dim)
    implicit none
    integer::n_dim,i
    real,dimension(n_dim**2) :: Q
    integer,dimension(8,n_dim**2) :: adj
    real :: chiral,k,qi,qj,qk,ql
    chiral=0.0
    do i=1,n_dim**2
            qi=Q(i)
            qj=Q(adj(4,i))
            qk=Q(adj(5,i))
            ql=Q(adj(6,i))
            k=(cos(qi-qj)-cos(qi-ql)-cos(qk-qj)+cos(qk-ql))/4
            chiral=chiral+k
        end do
end function chiral

function pr(q1,q2,pi)
    implicit none 
    real::q1,q2,pi,pr
    if(q1-q2>pi) then
        pr=q1-q2-2*pi
    elseif(q1-q2<-1*pi) then
        pr=q1-q2+2*pi
    else
        pr=q1-q2
    end if
end function pr


function vortex(Q,adj,sblt,n_dim,pi)
    implicit none
    integer:: n_dim,i,loc,vortex,vorticity
    real::Q(n_dim**2),q1,q2,q3,q4,pr,curl,pi
    integer::adj(8,n_dim**2),sblt(4,n_dim**2/4)
    vortex=0
    vorticity=0
    do i=1,n_dim**2/4
        loc=sblt(4,i)
        q1=Q(loc)
        q2=pr(Q(adj(3,loc)),-1*pi,pi)
        q3=Q(adj(5,adj(3,loc)))
        q4=pr(Q(adj(5,loc)),-1*pi,pi)
        curl=pr(q2,q1,pi)+pr(q3,q2,pi)+pr(q4,q3,pi)+pr(q1,q4,pi)
        if(curl>0) then
            vortex=vortex+1
            vorticity=vorticity+1
        elseif(curl<0) then
            vortex=vortex+1
            vorticity=vorticity-1
        end if
    end do
end function vortex










