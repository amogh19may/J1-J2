program send
    implicit none
    include "mpif.h"
    integer::ranks,ierr,a(4),i,b,recv_buff(16),send_buff(4)
    integer::state(MPI_STATUS_SIZE)
    a=[1,2,3,4]
    send_buff=0
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,ranks)
    do b=1,3
    if(ranks==0) then
	a(1)=0
    else  
	send_buff(ranks+1)=-ranks-1
        
    end if
    call MPI_Gather(send_buff,4,MPI_INT,recv_buff,4,MPI_INT,0,MPI_Comm_WORLD,ierr)
    if(ranks==0) then
	do i=1,3
	    a=a+recv_buff(4*(i)+1:4*(i+1))
	end do
    end if
    call MPI_Bcast(a,4,MPI_INT,0,MPI_Comm_WORLD,ierr)
    write(6,*) ranks," has ", a
    end do
  call MPI_finalize(ierr)
end program

