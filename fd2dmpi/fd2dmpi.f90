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


! 
! Main program for forward and adjoint modeling.
!
program fd2dmpi

use solver
use parser

implicit none

character*256 :: parfile

    call readCmd('par',parfile,'fdparfile')
    call readParFile(parfile,'jobtype',jobtype)

    if (jobtype(1:16) == 'forward_modeling') then
        call forward_modeling(parfile)
    else if (jobtype(1:16) == 'adjoint_modeling') then
        call adjoint_modeling(parfile)
    endif
 
end program fd2dmpi

