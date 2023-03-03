!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  The subroutines below are adapted from the Computational Toolkit provided in: 
!   
!      Schuster, G. T. (2017). Seismic inversion. Society of Exploration Geophysicists.
!
!  We kindly thank Prof. Schuster for allowing us to use these useful and efficient Fortran 
!  subroutines. Please Cite the book above if you use these subroutines.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



module mmi_mpi

use mpi
implicit none
!include 'mpif.h'
integer  :: root, rank, nsize, tag, status(MPI_STATUS_SIZE), ierr, sendtag, recvtag, currentShot, &
            source_tag, data_tag, terminate_tag, result_tag, slave, count

contains

!------------------------------------------------------------------------------
subroutine start_mpi

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nsize,ierr)

end subroutine start_mpi

!------------------------------------------------------------------------------
subroutine stop_mpi

call MPI_Finalize(ierr)

end subroutine stop_mpi

!-----------------------------------------------------------------------------------------
! Static load balancing
!
! Purpose: Assign starting and ending shot numbers to an MPI process
!
! Modified from the original version written by Jianming Sheng
!
! is1 = starting shot number of the original list of shot numbers
! is2 = ending shot number of the original list of shot numbers
! is_begin = starting shot number for an MPI process
! is2_end = ending shot number for an MPI process
!
subroutine get_assigned(is1, is2, is_begin, is_end)

integer, intent(in)  :: is1, is2
integer, intent(out) :: is_begin, is_end
integer              :: nnn, ii1, ii2

    if(nsize == 1)then
        is_begin = is1
        is_end = is2
        return
    endif
    
    nnn = is2-is1+1
    ii1 = int(nnn/nsize)
    if(ii1 == 0)then
        is_begin = is1+rank
        is_end = is_begin
        return
    endif
    
    ii2=nnn-ii1*nsize
    if(rank < (nsize-ii2))then
        is_begin = is1+rank*ii1
        is_end   = is1+(rank+1)*ii1-1
    else 
        is_begin = is1+ii2+rank*ii1+rank-nsize
        is_end   = is1+ii2+(rank+1)*ii1+rank-nsize
    endif
end subroutine get_assigned

end module mmi_mpi

