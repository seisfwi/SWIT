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
module solver

use fdcore
use datatype

implicit none

type(param),       private              :: par
type(acquisition), private              :: coord
logical,           private              :: store_snap, store_boundary
integer,           private              :: ix, iz, it, is, is1, is2, ig, igx, igz
integer,           private, allocatable :: fs(:)
real,              private, allocatable :: s(:), c(:,:), den(:,:), s_adj(:,:), csg(:,:), &
                                           p_end(:,:), u_end(:,:), w_end(:,:), boundary(:,:), &
                                           p0_bk(:,:), p_bk(:,:), u_bk(:,:), w_bk(:,:), &
                                           g(:,:), fw(:,:), bk(:,:)

                                           
contains

  !-----------------------------------------------------------------------------------------
  ! Forward modeling
  subroutine forward_modeling(parfile)

    ! Load modules
    use io
    use math
    use mmi_mpi
    use pml
    use source

    ! Parameters
    character(len=*), intent(in) :: parfile
    character(len=200) :: str

    ! Start MPI
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

    ! Save snapshots or not
    store_snap = .false.
    store_boundary = .false.
    if (par%store_snap .eq. 1) store_snap = .true.

    ! Forward modeling using dynamic load balancing
    call get_assigned(1, coord%ns, is1, is2)
    do is = is1, is2, 1
      
      ! Read source wavelet ./config/wavelet/src1.bin, src2.bin, ...
      call filename(str, par%sourcefile, is, '.bin')
      call read_binfile(str, s, par%nt)
      
      ! Forward modeling
      call staggered42_modeling_is(is, par, coord, s, c, den, fs, nx_pml, nz_pml, & 
                                   npml, damp_global, store_boundary, store_snap)
    
    enddo


    ! Synchronize all processes   
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! Deallocate memory
    deallocate(c, den, s, fs, damp, damp_global)

    999 continue

    ! Stop MPI
    call stop_mpi

  end subroutine forward_modeling


  !--------------------------------------------------------------------------
  ! Adjoint modeling
  subroutine adjoint_modeling(parfile)

    ! Load modules
    use io
    use math
    use mmi_mpi
    use pml
    use source
    
    ! Parameters
    character(len=*), intent(in) :: parfile
    character(len=200)      :: str
    
    ! Start MPI
    call start_mpi

    ! Read input parameters
    call readparamfile(parfile, par)

    ! PML setting: damp & damp_global
    call init_pml(par%nx, par%nz, par%npml)

    ! Read acquisition geometry data
    call readcoordfile(par%coordfile, coord)

    ! Memory allocations
    allocate(fs(nx_pml))
    allocate(c(nz_pml,nx_pml), den(nz_pml,nx_pml))
    allocate(p_bk(nz_pml,nx_pml))
    allocate(u_bk(nz_pml,nx_pml))
    allocate(w_bk(nz_pml,nx_pml))
    allocate(csg(par%nt,coord%ngmax))
    allocate(s_adj(par%nt, coord%ngmax))

    ! Set up free surface
    call free_surface(par, fs, npml)

    ! Read velocity model
    call readvelfile(par,c,npml,nx_pml,nz_pml)

    ! Read density model
    call read_densityfile(par,c,den,npml,nx_pml,nz_pml)
    
    ! Determine minimum velocity
    par%vmin = min_value(c(iz1:iz2,ix1:ix2), par%nz,par%nx,fs-npml)

    ! Setup PML damping coefficient
    call setup_pml(par%dx, par%vmin)

    ! Output adjoint modeling information
    if (rank==0) then
      write(*,*) '************************************************************'
      write(*,*) ' '
      write(*,*) ' Adjoint modeling ---- time-domain staggered-42-FD          '
      write(*,*) ' '
      write(*,*) '************************************************************'
    endif

    ! Adjoint modeling using dynamic load balancing
    store_snap = .false.
    if (par%store_snap .eq. 1) store_snap = .true.

    ! Adjoint modeling using dynamic load balancing
    call get_assigned(1, coord%ns, is1, is2)
    do is = is1, is2, 1
      
      ! Initialize
      s_adj = 0.0
      p_bk = 0.0
      u_bk = 0.0
      w_bk = 0.0
      csg = 0.0
      
      ! Read adjoint sources: ./config/wavelet/src1_adj.bin, src2_adj.bin, ...
      call filename(str, par%sourcefile, is, '_adj.bin')
      call read_binfile(str, s_adj, par%nt, coord%ngmax) 

      ! Time loop
      do it=par%nt,2,-1
        
        ! Adjoint modeling
        call staggered42_back_it(is, it, par, coord, s_adj, c, den, fs, nx_pml, &
                       nz_pml, npml, damp_global, p_bk, u_bk, w_bk, store_snap)

        ! Output pressure seismogram
        do ig=1,coord%ng(is)
          igx = npml + int(coord%xg(is,ig)/par%dx)+1
          igz = npml + int(coord%zg(is,ig)/par%dx)+1
          if (igz == fs(igx)) igz = igz + 1
          csg(it,ig) = p_bk(igz,igx)
        enddo

        ! Store wavefield snapshots
        if (mod(it, 10) == 0 .and. store_snap) then
          call filename(str, par%data_out, is, '/pressure')
          call filename(str, str, it, '_snapshot_adj.bin')
          call write_binfile(str, p_bk, par%nz, par%nx)
        endif

      enddo

      ! Export seismogram in SU format or binary format
      if (par%fileformat(1:2) == 'su') then
        call filename(output, par%data_out, is, '_sg.su')
        call write_sufile(output, csg(:,1:coord%ng(is)), par%nt, coord%ng(is), &
                          par%dt, 1.0, coord%xs(is), coord%xg(is,:))
      else
        call filename(output, par%data_out, is, '_sg.bin')
        call write_binfile(output, csg(:,1:coord%ng(is)), par%nt, coord%ng(is))
      endif

    enddo

    ! Synchronize all processes   
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! Deallocate memory
    deallocate(c, den, s_adj, fs, damp, damp_global)
    deallocate(p_bk, u_bk, w_bk)

    999 continue

    ! Stop MPI
    call stop_mpi

    end subroutine adjoint_modeling
    


  !-----------------------------------------------------------------------------------------
  ! This subroutine computes the gradient of the objective function based on the adjoint method
  subroutine gradient_computing(parfile)

    ! Use modules
    use io
    use math
    use mmi_mpi
    use pml
    use source
    
    ! Local variables
    character(len=*), intent(in) :: parfile
    character(len=200)      :: str
    
    ! Set up MPI
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

    ! Output information 
    if (rank==0) then
      write(*,*) '************************************************************'
      write(*,*) ' '
      write(*,*) ' Gradient Ccomputing ---- time-domain staggered-42-FD       '
      write(*,*) ' '
      write(*,*) '************************************************************'
    endif

    ! Adjoint modeling using dynamic load balancing
    store_snap = .false.
    store_boundary = .true.

    ! Loop over sources
    call get_assigned(1, coord%ns, is1, is2)
    do is = is1, is2, 1

      ! Initialize
      s = 0.0
      s_adj = 0.0
      boundary = 0.0
      p_end = 0.0
      u_end = 0.0
      w_end = 0.0
      p_bk = 0.0
      u_bk = 0.0
      w_bk = 0.0
      g = 0.0
      fw = 0.0
      bk = 0.0

      ! Read forward source: ./config/wavelet/src1.bin, src2.bin, ...
      call filename(str, par%sourcefile, is, '.bin')
      call read_binfile(str, s, par%nt)
      
      ! Forward modeling
      call staggered42_modeling_is(is, par, coord, s, c, den, fs, nx_pml, nz_pml, &
                  npml, damp_global, store_boundary, store_snap)
    
      ! Read adjoint source
      call filename(str, par%sourcefile, is, 'adj_.bin')
      call read_binfile(str, s_adj, par%nt, coord%ngmax) 
      
      ! Read boundary, p_end, u_end, w_end
      call filename(str, par%data_out, is, '/boundary.bin')
      call read_binfile(str, boundary, par%nx*9+par%nz*6, par%nt)
      call filename(str, par%data_out, is, '/p_end.bin')
      call read_binfile(str, p_end, par%nz, par%nx)
      call filename(str, par%data_out, is, '/u_end.bin')
      call read_binfile(str, u_end, par%nz, par%nx)
      call filename(str, par%data_out, is, '/w_end.bin')
      call read_binfile(str, w_end, par%nz, par%nx)

      ! Time loop 
      do it=par%nt,2,-1

        p0_bk = p_bk
        
        ! Backward modeling
        call staggered42_reco_it(is, it, par, coord, s, c(iz1:iz2,ix1:ix2), &
                      den(iz1:iz2,ix1:ix2), fs, boundary, p_end, u_end, w_end)
        
        ! Adjoint modeling
        call staggered42_back_it(is, it, par, coord, s_adj, c, den, fs, nx_pml, & 
                      nz_pml, npml, damp_global, p_bk, u_bk, w_bk, store_snap)

        ! Gradient
        g  = g  + p_end * (p0_bk(iz1:iz2,ix1:ix2)-p_bk(iz1:iz2,ix1:ix2))
        fw = fw + p_end * p_end
        bk = bk + p_bk(iz1:iz2,ix1:ix2) *  p_bk(iz1:iz2,ix1:ix2)

      enddo

      ! Write gradient and illumination

      ! perform the proper scale to the gradient 
      g = -2.0 * g / c(iz1:iz2,ix1:ix2)

      call filename(output, par%data_out, is, '/vp_gradient.bin')
      call write_binfile(output,  g, par%nz, par%nx) 

      call filename(output, par%data_out, is, '/forward_illumination.bin')
      call write_binfile(output, fw, par%nz, par%nx) 

      call filename(output, par%data_out, is, '/adjoint_illumination.bin')
      call write_binfile(output, bk, par%nz, par%nx) 

    enddo

    ! Synchronize all processes   
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! Deallocate memory
    deallocate(c, den, s, s_adj, fs, damp, damp_global)
    deallocate(boundary, p_end, u_end, w_end, p0_bk, p_bk, u_bk, w_bk, g, fw, bk)

    999 continue

    ! Stop MPI
    call stop_mpi

  end subroutine gradient_computing

end module solver

