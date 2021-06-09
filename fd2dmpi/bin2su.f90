! Program for converting a binary file into a SU format.
!
program bin2su

use datatype
use io
use mmi_mpi
use parser
use string
use su

implicit none

type(acquisition) :: coord
character(200)    :: coordfile, csg_in, data_out, input, output
integer           :: is, ig, nt_in, nt_out, iargc
real, allocatable :: csg(:,:)

if (iargc() == 0) then
  write(*,*) 'Usage: bin2su.x coordfile csg_in data_out nt_in nt_out'
  goto 999
endif

! Read input parameters from command line
call readCmd('coordfile', coordfile)
call readCmd('csg_in', csg_in)
call readCmd('data_out', data_out)
call readCmd('nt_in', nt_in)
call readCmd('nt_out', nt_out)

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nsize,ierr)

call readcoordfile(coordfile, coord)

! Allocate memory
allocate(csg(nt_in,coord%ngmax))

! Convert data
do is=1,coord%ns
  write(*,*) 'is = ', is

  ! Read .bin file
  call filename(input, csg_in, is, '.bin');
  call read_binfile(input, csg, nt_in, coord%ng(is))

  ! Write .su file
  call filename(output, data_out, is, '.su');
  call write_sufile(output, csg(1:nt_out,1:coord%ng(is)), nt_out, coord%ng(is), 1.0, 1.0)
enddo

call MPI_Finalize(MPI_COMM_WORLD, ierr)

999 continue

end program bin2su

