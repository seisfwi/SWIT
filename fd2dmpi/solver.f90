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
real,              private, allocatable :: s(:), c(:,:), den(:,:), s_adj(:,:), &
                                           p_end(:,:), u_end(:,:), w_end(:,:), boundary(:,:), &
                                           p0_bk(:,:), p_bk(:,:), u_bk(:,:), w_bk(:,:), &
                                           dg(:,:), dfw(:,:), dbk(:,:), g(:,:), fw(:,:), bk(:,:)

                                           
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
    write(*,*) ' Forward modeling ---- time-domain staggered-42-FD          '
    write(*,*) ' '
    write(*,*) '************************************************************'
  endif

  ! Forward modeling using dynamic load balancing
  message = .true.
  store_snap = .false.
  if (par%store_snap .eq. 1) store_snap = .true.

  call get_assigned(1, coord%ns, is1, is2)

  do is = is1, is2, 1
    
    call filename(str,'./parfile/forward_source/src',is,'.bin')
    call read_binfile(str, s, par%nt)

    if (message) then
      write(*,*) 'Process ', rank, ', shot', is, 'source max value = ', maxval(s)
      flush(6)
    endif

    call staggered42_modeling_is(is, par, coord, s, c, den, fs, nx_pml, nz_pml, npml, damp_global, store_snap)
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
  allocate(s(par%nt))
  allocate(s_adj(par%nt, coord%ngmax))
  allocate(fs(nx_pml))
  allocate(c(nz_pml,nx_pml), den(nz_pml,nx_pml))

  allocate(boundary(par%nx*9+par%nz*6,par%nt))
  allocate(p_end(par%nz,par%nx))
  allocate(u_end(par%nz,par%nx))
  allocate(w_end(par%nz,par%nx))

  allocate(p0_bk(nz_pml,nx_pml))
  allocate(p_bk(nz_pml,nx_pml))
  allocate(u_bk(nz_pml,nx_pml))
  allocate(w_bk(nz_pml,nx_pml))

  allocate(dg(par%nz,par%nx))
  allocate(dfw(par%nz,par%nx))
  allocate(dbk(par%nz,par%nx))
  allocate(g(par%nz,par%nx))  
  allocate(fw(par%nz,par%nx))
  allocate(bk(par%nz,par%nx))

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
    write(*,*) ' Adjoint modeling ---- time-domain staggered-42-FD          '
    write(*,*) ' '
    write(*,*) '************************************************************'
  endif
  ! Adjoint modeling using dynamic load balancing
  message = .true.
  store_snap = .false.
  if (par%store_snap .eq. 1) store_snap = .true.

  call get_assigned(1, coord%ns, is1, is2)

  do is = is1, is2, 1
    
    ! Initialize
    s     = 0.0
    s_adj = 0.0
    boundary = 0.0
    p_end = 0.0
    u_end = 0.0
    w_end = 0.0

    dg   = 0.0
    dfw  = 0.0
    dbk  = 0.0

    p_bk = 0.0
    u_bk = 0.0
    w_bk = 0.0

    ! Read forward source 
    call filename(str,'./parfile/forward_source/src',is,'.bin')
    call read_binfile(str, s, par%nt)
    
    ! Read adjoint source
    call filename(str,'./parfile/adjoint_source/src',is,'.bin')
    call read_binfile(str, s_adj, par%nt, coord%ngmax) 
    
    ! Read boundary, p_end, u_end, w_end
    call filename(str, par%data_out, is, '_snapshot/boundary.bin')
    call read_binfile(str, boundary, par%nx*9+par%nz*6, par%nt)
    call filename(str, par%data_out, is, '_snapshot/p_end.bin')
    call read_binfile(str, p_end, par%nz, par%nx)
    call filename(str, par%data_out, is, '_snapshot/u_end.bin')
    call read_binfile(str, u_end, par%nz, par%nx)
    call filename(str, par%data_out, is, '_snapshot/w_end.bin')
    call read_binfile(str, w_end, par%nz, par%nx)

    do it=par%nt,2,-1

      p0_bk = p_bk
      
      call staggered42_reco_it(is, it, par, coord, s, c(iz1:iz2,ix1:ix2), den(iz1:iz2,ix1:ix2), fs, boundary, p_end, u_end, w_end)
      call staggered42_back_it(is, it, par, coord, s_adj, c, den, fs, nx_pml, nz_pml, npml, damp_global, p_bk, u_bk, w_bk)
    
      dg  = dg  + p_end * (p0_bk(iz1:iz2,ix1:ix2)-p_bk(iz1:iz2,ix1:ix2))
      dfw = dfw + p_end * p_end
      dbk = dbk + p_bk(iz1:iz2,ix1:ix2) *  p_bk(iz1:iz2,ix1:ix2)

      ! ! output wavefield
      ! if (mod(it, 10) ==0) then
      !   call filename(str, par%data_out, is, '_snapshot/')
      !   call filename(str, str, it, '_forward.bin')
      !   call write_binfile(str, p_end,  par%nz, par%nx)

      !   call filename(str, par%data_out, is, '_snapshot/')
      !   call filename(str, str, it, '_backward.bin')
      !   call write_binfile(str, p_bk(iz1:iz2,ix1:ix2),  par%nz, par%nx)

      ! endif

    enddo

    ! ! output gradient
    ! call filename(str, par%data_out, is,'_kernel_vp.bin')
    ! call write_binfile(str, dg,  par%nz, par%nx)

  enddo

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!   
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(dg,   g, par%nz*par%nx, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(dfw, fw, par%nz*par%nx, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(dbk, bk, par%nz*par%nx, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

  if (rank == 0) then
    
    ! perform the proper scale to the gradient 
    g = -2.0 * g / c(iz1:iz2,ix1:ix2)

    call filename(output, par%data_out, 0, '_kernel_vp.bin')
    call write_binfile(output,  g, par%nz, par%nx) 

    call filename(output, par%data_out, 0, '_illum_forw.bin')
    call write_binfile(output, fw, par%nz, par%nx) 

    call filename(output, par%data_out, 0, '_illum_back.bin')
    call write_binfile(output, bk, par%nz, par%nx) 

  endif

  deallocate(c, den, s, s_adj, fs, damp, damp_global)
  deallocate(boundary, p_end, u_end, w_end, p0_bk, p_bk, u_bk, w_bk, dg, dfw, dbk, g, fw, bk)

  999 continue
  call stop_mpi

  end subroutine adjoint_modeling
  

end module solver
