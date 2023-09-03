

program calling_func
   implicit none
   include "mpif.h"
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


   integer,parameter::rel=20
   integer,parameter::n_dim=20
   integer,parameter::n_save=50000
   integer,parameter::n_eqb=50000
   real,parameter::pi=4*atan(1.0)
   real,parameter::J2=0.7
   integer::i,j,k,ranks,n_size,n_proc,ierr,istat,isbl,ir,lc
   real :: T_arr(100)=[(i,i=1,100,1)]*0.0002+0.55
   real :: twist(11)=[(i,i=-5,5,1)]*pi/100
   real :: avK2,avK4,chi,chiral,flip,t_start,t_end
   integer::nbr_all(8,n_dim**2),latt(4,n_dim**2/4)
   real:: Q(n_dim**2),U(rel),send_buff(n_dim**2)
   real,allocatable::recv_buff(:)


   nbr_all=nbr_vec(n_dim)
   latt=sub_latt(n_dim)
   


   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_Comm_WORLD,ranks,ierr)
   call MPI_Comm_size(MPI_Comm_WORLD,n_size,ierr)
   allocate(recv_buff(n_size*n_dim**2))
   n_proc=n_dim**2/(4*n_size)
   call MPI_Barrier(MPI_Comm_WORLD,ierr)
   do k=51,60
    do j=1,rel
        avK2=0.0
        avK4=0.0
        chi=0.0
        if(ranks==0) then
            call random_number(Q)
            Q=2*pi*Q-pi
        end if
  
        call MPI_Bcast(Q,n_dim**2,MPI_REAL,0,MPI_Comm_WORLD,ierr)
	
        do i=1,n_eqb
            do isbl=1,4
                send_buff=0
                if(ranks==0) then
                    do ir=1,n_proc
			
                        Q(latt(isbl,ir))=Q(latt(isbl,ir))+flip(latt(isbl,ir),Q,T_arr(k),J2,n_dim,nbr_all,pi,0.0)
                    end do
		
		    
                else 
                    do ir=n_proc*ranks+1,n_proc*(ranks+1)
			
                        send_buff(latt(isbl,ir))=flip(latt(isbl,ir),Q,T_arr(k),J2,n_dim,nbr_all,pi,0.0)
			
                    end do
		
                end if
                call MPI_Gather(send_buff,n_dim**2,MPI_REAL,recv_buff,n_dim**2,MPI_REAL,0,MPI_Comm_WORLD,ierr)
		
                if(ranks==0) then
                    do ir=1,n_size-1
                        Q=Q+recv_buff(ir*n_dim**2+1:(ir+1)*n_dim**2)
			
                    end do
		
		end if
		    call MPI_Bcast(Q,n_dim**2,MPI_REAL,0,MPI_Comm_WORLD,ierr)
		
            end do
        end do
        do i=1,n_save
	!do lc=1,10
            do isbl=1,4
                send_buff=0
                if(ranks==0) then
                    do ir=1,n_proc
			
                        Q(latt(isbl,ir))=Q(latt(isbl,ir))+flip(latt(isbl,ir),Q,T_arr(k),J2,n_dim,nbr_all,pi,0.0)
                    end do

		    
                else 
                    do ir=n_proc*ranks+1,n_proc*(ranks+1)
			
                        send_buff(latt(isbl,ir))=flip(latt(isbl,ir),Q,T_arr(k),J2,n_dim,nbr_all,pi,0.0)
			
                    end do
                end if
                call MPI_Gather(send_buff,n_dim**2,MPI_REAL,recv_buff,n_dim**2,MPI_REAL,0,MPI_Comm_WORLD,ierr)
                if(ranks==0) then
                    do ir=1,n_size-1
                        Q=Q+recv_buff(ir*n_dim**2+1:(ir+1)*n_dim**2)
			
                    end do
		end if
		    call MPI_Bcast(Q,n_dim**2,MPI_REAL,0,MPI_Comm_WORLD,ierr)
            end do
            !end do
            if(ranks==0) then
                chi=chiral(Q,nbr_all,n_dim)
                avk2=avk2+chi**2
                avk4=avk4+chi**4
		
            end if
        end do
        if(ranks==0) then
            U(j)=1-n_save*avK4/(3*(avK2**2))
        end if

   end do
   if(ranks==0) then
    write(*,*) T_arr(k),sum(U)/(rel),sqrt((sum(U**2)-(sum(U)**2)/rel)/rel)
   end if
   end do
   call MPI_Barrier(MPI_Comm_WORLD,ierr)
   deallocate(recv_buff)
   call MPI_Finalize(ierr)
 

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

function flip(ind,Q,T,J2,n_dim,all_nbr,pi,twist)
    implicit none
    integer:: n_dim,i,j,ind
    integer::all_nbr(8,n_dim**2)
    real :: T,pi,J2,twist,dE,w,p,delE1,flip
    real,dimension(n_dim**2) :: Q
            call random_number(w)
            call random_number(p)
            w=2*pi*w-pi
            flip=0.0
            dE=delE1(ind,all_nbr,Q,J2,w,n_dim,twist)
            if (exp(-1*dE/T)>p) then
                flip=w-Q(ind)
            end if
end function flip



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










