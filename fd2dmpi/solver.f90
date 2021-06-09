! module for 2D applications including forward and adjoint modeling
!
module solver

use fdcore
use datatype

implicit none

type(param),       private              :: par
type(acquisition), private              :: coord
logical,           private              :: store_snap, message
integer,           private              :: ix, iz, it, is, is1, is2
integer,           private, allocatable :: fs(:)
real,              private, allocatable :: s(:), c(:,:), den(:,:), adj_src(:,:)

                                           
contains

!-----------------------------------------------------------------------------------------
subroutine forward_modeling(parfile)

use io
use math
use mmi_mpi
use pml
use source

character(len=*), intent(in) :: parfile
character(len=200)      :: str

call start_mpi

! Read input parameters
call readparamfile(parfile, par)
! PML setting: damp & damp_global
call init_pml(par%nx, par%nz, par%npml)
! Read acquisition geometry data
call readcoordfile(par%coordfile, coord)

! Memory allocations
allocate(s(par%nt))
allocate(fs(nx_pml))
allocate(c(nz_pml,nx_pml), den(nz_pml,nx_pml))

! Set up free surface
call free_surface(par, fs, npml)
! Read velocity model
call readvelfile(par,c,npml,nx_pml,nz_pml)
! Read density model
call read_densityfile(par,c,den,npml,nx_pml,nz_pml)
! Determine minimum velocity
par%vmin = min_value(c(iz1:iz2,ix1:ix2),par%nz,par%nx,fs-npml)
! Setup PML damping coefficient
call setup_pml(par%dx, par%vmin)

if (rank==0) then
  write(*,*) '************************************************************'
  write(*,*) ' '
  write(*,*) ' Forward modeling using a time-domain staggered-42-FD method'
  write(*,*) ' '
  write(*,*) '************************************************************'
endif
! Forward modeling using dynamic load balancing
message = .true.
store_snap = .false.
if (par%store_snap .eq. 1) store_snap = .true.

call get_assigned(1, coord%ns, is1, is2)
do is = is1, is2, 1
  ! Read source file
  call filename(str,'./parfile/forward_source/src',is,'.bin')
  call read_binfile(str, s, par%nt)
  if (message) then
    write(*,*) 'Process ', rank, ', shot', is, 'source max value = ', maxval(s)
    flush(6)
  endif
  ! staggered FD scheme 
  call forward_modeling_staggered42_is(is, par, coord, s, c, den, fs, nx_pml, nz_pml, npml, damp_global, store_snap)

enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!   
call MPI_Barrier(MPI_COMM_WORLD, ierr)

deallocate(c, den, s, fs, damp, damp_global)
999 continue
call stop_mpi

end subroutine forward_modeling


!-----------------------------------------------------------------------------------------
subroutine adjoint_modeling(parfile)

    use io
    use math
    use mmi_mpi
    use pml
    use source
    
    character(len=*), intent(in) :: parfile
    character(len=200)      :: str
    
    call start_mpi
    ! Read input parameters
    call readparamfile(parfile, par)
    ! PML setting: damp & damp_global
    call init_pml(par%nx, par%nz, par%npml)
    ! Read acquisition geometry data
    call readcoordfile(par%coordfile, coord)
    ! Memory allocations
    allocate(adj_src(par%nt, coord%ngmax))
    allocate(fs(nx_pml))
    allocate(c(nz_pml,nx_pml), den(nz_pml,nx_pml))

    ! Set up free surface
    call free_surface(par, fs, npml)
    ! Read velocity model
    call readvelfile(par,c,npml,nx_pml,nz_pml)
    ! Read density model
    call read_densityfile(par,c,den,npml,nx_pml,nz_pml)
    
    ! Determine minimum velocity
    par%vmin = min_value(c(iz1:iz2,ix1:ix2),par%nz,par%nx,fs-npml)
    ! Setup PML damping coefficient
    call setup_pml(par%dx, par%vmin)
 
    if (rank==0) then
      write(*,*) '************************************************************'
      write(*,*) ' '
      write(*,*) ' Adjoint modeling using a time-domain staggered-42-FD method'
      write(*,*) ' '
      write(*,*) '************************************************************'
    endif
    ! Adjoint modeling using dynamic load balancing
    message = .true.
    store_snap = .false.
    if (par%store_snap .eq. 1) store_snap = .true.

    call get_assigned(1, coord%ns, is1, is2)
    do is = is1, is2, 1
      ! Read adjoint source
      call filename(str,'./parfile/adjoint_source/src',is,'.bin')
      call read_binfile(str, adj_src, par%nt, coord%ngmax) 
      if (message) then
        write(*,*) 'Process ', rank, ', shot', is
        flush(6)
      endif
      call backward_modeling_staggered42_is(is, par, coord, adj_src, c, den, fs, &
                                            nx_pml, nz_pml, npml, damp_global, ix1, ix2, iz1, iz2, store_snap)
      
    enddo 

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!   
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    deallocate(c,den,fs,damp, damp_global, adj_src)
    999 continue
    
    call stop_mpi
    
    end subroutine adjoint_modeling
    




end module solver
