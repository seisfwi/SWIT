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
    
    ! Check the command line argument
    if (iargc() < 1) then
        write(*,*) 'Error: No parameter file specified.'
        write(*,*) '    Usage: fd2dmpi config=solver.config'
        stop 
    endif

    ! Read the parameter file
    call readCmd('config', parfile, 'solver.config')
    call readParFile(parfile,'jobtype',jobtype)
    
    ! Forward modeling
    if (jobtype(1:7) == 'forward') then
        call forward_modeling(parfile)

    ! Adjoint modeling
    else if (jobtype(1:7) == 'adjoint') then
        call adjoint_modeling(parfile)

    ! Gradient computation
    else if (jobtype(1:8) == 'gradient') then
        call gradient_computing(parfile)
    endif
 
end program fd2dmpi

