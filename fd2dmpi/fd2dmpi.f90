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

