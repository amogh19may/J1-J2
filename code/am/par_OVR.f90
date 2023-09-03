

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


   integer,parameter::rel=2
   integer,parameter::n_dim=16
   integer,parameter::n_save=1000
   integer,parameter::n_eqb=1000
   real,parameter::pi=4.0*atan(1.0)
   real,parameter::J2=0.7
   integer::i,j,k,ranks,n_size,n_proc,ierr,istat,isbl,ir,lovr,lmc,lm2
   real :: T_arr(100)=[(i,i=1,100,1)]*0.0002+0.55
   real :: avp1,avp2,avp3,avE,avV,M1x,M1y,M2x,M2y,M,chi,E,chiral,energy,flip,loc_field
   integer::nbr_all(8,n_dim**2),latt(4,n_dim**2/4)
   real:: Q(n_dim**2),prop1(rel),prop2(rel),send_buff(n_dim**2),stag(n_dim**2),prop_send(4),prop_recv(16)
   real,allocatable::recv_buff(:)

   nbr_all=nbr_vec(n_dim)
   latt=sub_latt(n_dim)



   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_Comm_WORLD,ranks,ierr)
   call MPI_Comm_size(MPI_Comm_WORLD,n_size,ierr)
   allocate(recv_buff(n_size*n_dim**2))
   n_proc=n_dim**2/(4*n_size)


   do k=51,60
	if(ranks==0) then
	end if
    do j=1,rel
	avp1=0.0
        avp2=0.0
        avp3=0.0
        chi=0.0
	M=0.0
        if(ranks==0) then
            call random_number(Q)
            Q=2*pi*Q-pi
        end if

        call MPI_Bcast(Q,n_dim**2,MPI_REAL,0,MPI_Comm_WORLD,ierr)

    do i=1,n_eqb/10

 	!start of MC step
    	do lmc=1,10
            do isbl=1,4
                send_buff=0.0
		recv_buff=0.0
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
	        send_buff=0.0
	        recv_buff=0.0
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
                send_buff=0.0
		recv_buff=0.0
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
		recv_buff=0.0
		send_buff=0.0
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
                chi=chiral(Q,nbr_all,n_dim)/(n_dim**2)
		avp1=avp1+chi
                avp2=avp2+chi**2
                avp3=avp3+chi**4
	    else if(ranks==1) then
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
		M=sqrt((abs(M1x)+abs(M2x))**2+(abs(M1y)+abs(M2y))**2)/(n_dim**2)
		avp1=avp1+M
		avp2=avp2+M**2
		avp3=avp3+M**4
             else if(ranks==2) then
		chi=chiral(Q,nbr_all,n_dim)/(n_dim**2)
		E=energy(Q,nbr_all,n_dim,J2)/(n_dim**2)
		avp1=avp1+E
		avp2=avp2+chi**2
		avp3=avp3+E*chi**2
	     else if(ranks==3) then
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
		M=sqrt((abs(M1x)+abs(M2x))**2+(abs(M1y)+abs(M2y))**2)/(n_dim**2)
		E=energy(Q,nbr_all,n_dim,J2)/(n_dim**2)
		avp1=avp1+E
		avp2=avp2+M**2
		avp3=avp3+E*M**2
            end if

	end do !end of n_save loop


        if(ranks==0 .or. ranks==1) then
            prop1(j)=1-n_save*avp3/(30.0*(avp2**2))
	    prop2(j)=(10.0*avp2/n_save-(10.0*avp1/n_save)**2)*n_dim**2/(T_arr(k))
	else if(ranks==2 .or. ranks==3) then
	    prop1(j)=avp3/avp2-10.0*avp1/n_save
            prop2(j)=10.0*avp2/n_save
        end if


   end do !end of rel loop

    prop_send(1)=sum(prop1)/rel
    prop_send(2)=sqrt((sum(prop1**2)-(sum(prop1)**2)/rel)/rel)
    prop_send(3)=sum(prop2)/rel
    prop_send(4)=sqrt((sum(prop2**2)-(sum(prop2)**2)/rel)/rel)

    !communicate and print
    call MPI_Gather(prop_send,4,MPI_REAL,prop_recv,4,MPI_REAL,0,MPI_Comm_WORLD,ierr)
    if(ranks==0) then
	write(*,*) T_arr(k),prop_recv
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
	energy=0.0
	do i=1,n_dim**2
	    qe=Q(i)
	    q1=Q(adj(1,i))
	    q2=Q(adj(2,i))
	    q3=Q(adj(3,i))
	    q4=Q(adj(4,i))
	    energy=energy-(cos(qe-q2)+cos(qe-q4))+J2*(cos(qe-q1)+cos(qe-q3))
	end do
end function energy




