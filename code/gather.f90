program hello

   implicit none
   include "mpif.h"
    
   integer n_ranks, rank, sender, ierr
   character(len=40) send_message
   character, dimension(:), allocatable :: receive_message

   ! First call MPI_Init
   call MPI_Init(ierr)

   ! Get my rank and the number of ranks
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, n_ranks, ierr)

   ! Allocate space for all received messages in receive_message
   allocate ( receive_message(n_ranks*40) )

   ! Use gather to send all messages to rank 0
   write(send_message,*) "Hello World, I'm rank", rank
   call MPI_Gather( send_message, 40, MPI_CHARACTER, receive_message, 40, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
   
   if (rank == 0) then
       do sender = 0, n_ranks-1
           write(6,*) receive_message(40*sender+1: 40*(sender+1))
       end do
   end if

   ! Free memory and finalise
   deallocate( receive_message )
   call MPI_Finalize(ierr)
end