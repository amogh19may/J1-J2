

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


   integer,parameter::n_dim=8
   integer,parameter::n_save=5000
   integer,parameter::n_eqb=5000
   real,parameter::pi=4*atan(1.0)
   real,parameter::J2=0.7
   integer::i,j,k,ranks,n_size,n_proc,ierr,istat,isbl,ir,lovr,lmc,lm2
   real :: T_arr(100)=[(i,i=1,100,1)]*0.001+0.5
   real :: M1x,M1y,M2x,M2y,M,chi,chiral,flip,t_start,t_end,loc_field,E,energy
   integer::nbr_all(8,n_dim**2),latt(4,n_dim**2/4)
   real:: Q(n_dim**2),send_buff(n_dim**2),prop(6),prop_buff(24)
   real,allocatable::recv_buff(:)



   nbr_all=nbr_vec(n_dim)
   latt=sub_latt(n_dim)
  


   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_Comm_WORLD,ranks,ierr)
   call MPI_Comm_size(MPI_Comm_WORLD,n_size,ierr)
   allocate(recv_buff(n_size*n_dim**2))
   n_proc=n_dim**2/(4*n_size)
	
   call MPI_Barrier(MPI_Comm_WORLD,ierr)
   t_start=MPI_Wtime()

   do k=1,5
	if(ranks==0) then
	end if
        chi=0.0
	M=0.0
	prop=0.0
        if(ranks==0) then
            call random_number(Q)
            Q=2*pi*Q-pi
        end if

        call MPI_Bcast(Q,n_dim**2,MPI_REAL,0,MPI_Comm_WORLD,ierr)

    do i=1,n_eqb/10
    !do lc=1,10
 	!start of MC step
    	do lmc=1,5 
            do isbl=1,4
                send_buff=0.0
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
 	   end do !end of lmc loop
            !end of MC step

            !start of OVR
          do lovr=1,5
            do isbl=1,4
	        send_buff=0
	        recv_buff=0
                do ir=n_proc*ranks+1,n_proc*(ranks+1)
                    send_buff(latt(isbl,ir))=loc_field(latt(isbl,ir),Q,J2,nbr_all,n_dim)
                end do
                call MPI_Gather(send_buff,n_dim**2,MPI_REAL,recv_buff,n_dim**2,MPI_REAL,0,MPI_Comm_WORLD,ierr)
                if(ranks==0) then
                    do ir=0,n_size-1
                        Q=Q+recv_buff(ir*n_dim**2+1:(ir+1)*n_dim**2)
                    end do
                end if
                call MPI_Bcast(Q,n_dim**2,MPI_REAL,0,MPI_Comm_WORLD,ierr)
            end do
         end do !end of lovr loop
            !end of OVR
          !end do !end of lc loop
            !equillibrium attained (hopefully!!)
        end do

        !start generating measurement configs
        do i=1,n_save/10
	!do lc=1,10

            !start of MC sweep
	do lmc=1,5	
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
	    end do !end of lmc loop 	
            !end of MC sweep

            !start of OVR
	do lovr=1,5
            do isbl=1,4
		recv_buff=0
		send_buff=0
                do ir=n_proc*ranks+1,n_proc*(ranks+1)
                    send_buff(latt(isbl,ir))=loc_field(latt(isbl,ir),Q,J2,nbr_all,n_dim)
                end do
                call MPI_Gather(send_buff,n_dim**2,MPI_REAL,recv_buff,n_dim**2,MPI_REAL,0,MPI_Comm_WORLD,ierr)
                if(ranks==0) then
                    do ir=0,n_size-1
                        Q=Q+recv_buff(ir*n_dim**2+1:(ir+1)*n_dim**2)
                    end do
                end if
                call MPI_Bcast(Q,n_dim**2,MPI_REAL,0,MPI_Comm_WORLD,ierr)
            end do !end of isbl loop
	end do !end of lovr loop
            !end of OVR

        !end do !end of lc loop
            !calculate chirality
            if(ranks==0) then
                chi=chiral(Q,nbr_all,n_dim)
		prop(1)=prop(1)+chi
		prop(2)=prop(2)+chi**2
		prop(3)=prop(3)+chi**3
		prop(4)=prop(4)+chi**4
		prop(5)=prop(5)+chi**6
		prop(6)=prop(6)+chi**8
	    else if(ranks==1) then
                chi=chiral(Q,nbr_all,n_dim)
		E=energy(Q,nbr_all,n_dim,J2)
		prop(1)=prop(1)+chi*E
		prop(2)=prop(2)+chi*(E**2)
		prop(3)=prop(3)+(chi**2)*E
		prop(4)=prop(4)+(chi**2)*(E**2)
		prop(5)=prop(5)+E	
	    else if(ranks==2) then
		M1x=0.0
		M1y=0.0
		M2x=0.0
		M2y=0.0
		do lm2=1,(n_dim**2)/4
		    M1x = M1x+cos(Q(latt(1,lm2)))-cos(Q(latt(4,lm2)))
		    M2x = M2x+cos(Q(latt(2,lm2)))-cos(Q(latt(3,lm2)))
		    M1y = M1y+sin(Q(latt(1,lm2)))-sin(Q(latt(4,lm2)))
		    M2y = M2y+sin(Q(latt(2,lm2)))-sin(Q(latt(3,lm2)))
		end do
		M=sqrt((abs(M1x)+abs(M2x))**2+(abs(M1y)+abs(M2y))**2)
		prop(1)=prop(1)+M
		prop(2)=prop(2)+M**2
		prop(3)=prop(3)+M**3
		prop(4)=prop(4)+M**4
		prop(5)=prop(5)+M**6
		prop(6)=prop(6)+M**8
	    else if(ranks==3) then
		M1x=0.0
		M1y=0.0
		M2x=0.0
		M2y=0.0
		E=energy(Q,nbr_all,n_dim,J2)
		do lm2=1,(n_dim**2)/4
		    M1x = M1x+cos(Q(latt(1,lm2)))-cos(Q(latt(4,lm2)))
		    M2x = M2x+cos(Q(latt(2,lm2)))-cos(Q(latt(3,lm2)))
		    M1y = M1y+sin(Q(latt(1,lm2)))-sin(Q(latt(4,lm2)))
		    M2y = M2y+sin(Q(latt(2,lm2)))-sin(Q(latt(3,lm2)))
		end do
		M=sqrt((abs(M1x)+abs(M2x))**2+(abs(M1y)+abs(M2y))**2)
		prop(1)=prop(1)+M*E
		prop(2)=prop(2)+M*(E**2)
		prop(3)=prop(3)+(M**2)*E
		prop(4)=prop(4)+(M**2)*(E**2)
		prop(5)=prop(5)+E**2
            end if  
        
	end do !end of n_save loop 
        call MPI_Gather(prop,6,MPI_REAL,prop_buff,6,MPI_REAL,0,MPI_Comm_WORLD,ierr)
	


   if(ranks==0) then
	print*,T_arr(k),"	",10*prop_buff/n_save		
   end if

   end do

   call MPI_Barrier(MPI_Comm_WORLD,ierr)
   t_end=MPI_Wtime()
   write(6,*) t_end-t_start
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
    real :: T,pi,J2,twist,dE,w,p,delE1,flip,Tin
    real,dimension(n_dim**2) :: Q
            call random_number(w)
            call random_number(p)
            w=2*pi*w-pi
	    Tin=1/T
            flip=0.0
            dE=delE1(ind,all_nbr,Q,J2,w,n_dim,twist)
            if (exp(-Tin*dE)>p) then
                flip=w-Q(ind)
            end if
end function flip

function loc_field(i,Qmat,J2,nbr,n_dim)
        implicit none
        real::q,qr,ql,qu,qd,qs,qn,qw,qe,sx,sy,loc_field,J2
        real::Qmat(n_dim**2)
        integer::i,n_dim,nbr(8,n_dim**2)
        q=Qmat(i)
        qr=Qmat(nbr(4,i))
        ql=Qmat(nbr(8,i))
        qu=Qmat(nbr(2,i))
        qd=Qmat(nbr(6,i))
        qs=Qmat(nbr(5,i))
        qn=Qmat(nbr(1,i))
        qw=Qmat(nbr(3,i))
        qe=Qmat(nbr(7,i))
        sx= -1*(cos(qr)+cos(ql)+cos(qu)+cos(qd))+J2*(cos(qs)+cos(qn)+cos(qw)+cos(qe))
        sy= -1*(sin(qr)+sin(ql)+sin(qu)+sin(qd))+J2*(sin(qs)+sin(qn)+sin(qw)+sin(qe))
        loc_field= 2*(atan2(sy,sx)-q)
end function loc_field

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

function energy(Q,adj,n_dim,J2)
	implicit none
	real::Q(n_dim**2),J2,energy,qe,q1,q2,q3,q4
	integer::adj(8,n_dim**2),n_dim,i
	energy=0
	do i=1,n_dim**2
	    qe=Q(i)
	    q1=Q(adj(1,i))
	    q2=Q(adj(2,i))
	    q3=Q(adj(3,i))
	    q4=Q(adj(4,i))
	    energy=energy-(cos(qe-q2)+cos(qe-q4))+J2*(cos(qe-q1)+cos(qe-q3))
	end do
end function energy	




