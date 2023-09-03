program hello

   implicit none
   include "mpif.h"
    
   integer, parameter :: n_numbers=4
   integer i
   integer rank, n_ranks, neighbour, ierr, req
   integer status(MPI_STATUS_SIZE)
   integer send_message(n_numbers)
   integer recv_message(n_numbers)

   ! First call MPI_Init
   call MPI_Init(ierr)

   ! Get my rank and the number of ranks
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, n_ranks, ierr)

   ! Check that there are exactly two ranks
   if (n_ranks .NE. 2) then
        write(6,*) "This example requires exactly two ranks"
        error stop
   end if
   recv_message=[4,3,2,1]
   ! Call the other rank the neighbour
   if (rank == 0) then
       neighbour = 1
   else
       neighbour = 0
   end if

   ! Generate numbers to send
   do i = 1, n_numbers
       send_message(i) = i;
   end do

   ! Send the message to other rank
   call MPI_ISend( send_message, n_numbers, MPI_INTEGER, neighbour, 0, MPI_COMM_WORLD, req,ierr )

   ! Receive the message from the other rank

   write(6,*) recv_message
   call MPI_IRecv( recv_message, n_numbers, MPI_INTEGER, neighbour, 0, MPI_COMM_WORLD,req, ierr )
   write(6,*) recv_message
   write(6,*) "Message received by rank", rank

   ! Call MPI_Finalize at the end
   call MPI_Finalize(ierr)
end