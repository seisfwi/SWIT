! module for 2D forward modeling using a time-domain FD method
!
module fdcore

use global
implicit none

integer, private              :: ix, iz, it, is, is1, is2, isx, isz, ig, igx, igz, i, j, jj
real,    private              :: alpha, beta, a, d1, d2, kappa
real,    private, allocatable :: u(:,:), w(:,:), p(:,:), p0(:,:), csg(:,:),  &
                                 p_for(:,:), kernel(:,:), illum_forward(:,:), illum_backward(:,:)

contains

!-----------------------------------------------------------------------------------------
subroutine forward_modeling_staggered42_is(is, par, coord, s, c, den, fs, nx_pml, nz_pml, npml, damp, store_snap)
  use io
  use datatype
  use math
  use source
  use string

  type(param), intent(in)       :: par
  type(acquisition), intent(in) :: coord
  logical, intent(in)           :: store_snap
  integer, intent(in)           :: is, fs(:), nx_pml, nz_pml, npml
  real, intent(in)              :: damp(:,:), c(:,:), den(:,:)
  real, intent(inout)           :: s(:)

  ! Memory allocations
  allocate(u(nz_pml,nx_pml-1))
  allocate(w(nz_pml-1,nx_pml))
  allocate(p(nz_pml,nx_pml))
  allocate(csg(par%nt,coord%ngmax))

  ! Initialize to zero
  u = 0.0
  w = 0.0
  p = 0.0
  csg = 0.0

  ! Setup source and receiver position
  isx = npml+int(coord%xs(is)/par%dx)+1
  isz = npml+int(coord%zs(is)/par%dx)+1
  if (isz == fs(isx)) isz = isz + 1

  ! Time marching loop
  do it=2,par%nt
    ! Update the pressure field
    do ix=3,nx_pml-2
      do iz=fs(ix)+1,nz_pml-2
        alpha = den(iz,ix)*c(iz,ix)*c(iz,ix)*par%dt/(par%dx)
        kappa = damp(iz,ix)*par%dt
        p(iz,ix) = (1.0-kappa)*p(iz,ix) - alpha* &
                  (c1_staggered*(u(iz,ix)-u(iz,ix-1) + w(iz,ix)-w(iz-1,ix)) + &
                  c2_staggered*(u(iz,ix+1)-u(iz,ix-2) + w(iz+1,ix)-w(iz-2,ix)))
      enddo
    enddo

    ! Add source
    beta = par%dt
    p(isz,isx) = p(isz,isx) + beta*s(it)

    ! Free surface
    if (par%free_surface == 1) then
      p(fs-1,:) = -p(fs+1,:)
    endif

    ! Update horizontal particle velocity: u
    do ix=2,nx_pml-2
      do iz=fs(ix),nz_pml-1
        alpha = par%dt/(den(iz,ix)*par%dx)
        kappa = 0.5*(damp(iz,ix)+damp(iz,ix+1))*par%dt
        u(iz,ix) = (1.0-kappa)*u(iz,ix) - alpha*(c1_staggered*(p(iz,ix+1)-p(iz,ix)) &
                                                +c2_staggered*(p(iz,ix+2)-p(iz,ix-1)))
      enddo
    enddo
    
    ! Update verticle particle velocity: w
    do ix=2,nx_pml-1
      do iz=fs(ix),nz_pml-2
        alpha = par%dt/(den(iz,ix)*par%dx)
        kappa = 0.5*(damp(iz,ix)+damp(iz+1,ix))*par%dt
        w(iz,ix) = (1.0-kappa)*w(iz,ix) - alpha*(c1_staggered*(p(iz+1,ix)-p(iz,ix)) &
                                                +c2_staggered*(p(iz+2,ix)-p(iz-1,ix)))
      enddo
    enddo

    ! Free surface boundary condition for velocity field
    if (par%free_surface.eq.1) then
      do ix=1,nx_pml
        w(fs(ix)-1,ix) = w(fs(ix),ix)
      enddo
    endif

    ! Output pressure seismogram
    do ig=1,coord%ng(is)
      igx = npml + int(coord%xg(is,ig)/par%dx)+1
      igz = npml + int(coord%zg(is,ig)/par%dx)+1
      if (igz == fs(igx)) igz = igz + 1
      csg(it,ig) = p(igz,igx)
    enddo
    
    ! Save snapshot
    if (mod(it,par%store_step)==0 .and. store_snap) then
      ! save path name
      call filename(output, par%data_out, is, '_snapshot/')
      call filename(output, output, it, '_snapshot.bin')
      call write_binfile(output, p, nz_pml, nx_pml)

    endif

  enddo  ! End of time-marching loop

  ! Export seismogram
  if (par%fileformat(1:2) == 'su') then
    call filename(output, par%data_out, is, '_sg.su')
    call write_sufile(output, csg(:,1:coord%ng(is)), par%nt, coord%ng(is), par%dt, 1.0, &
      coord%xs(is), coord%xg(is,:))
  else
    call filename(output, par%data_out, is, '_sg.bin')
    call write_binfile(output, csg(:,1:coord%ng(is)), par%nt, coord%ng(is))
  endif

  deallocate(u, w, p, csg)

  999 continue

end subroutine forward_modeling_staggered42_is



!-----------------------------------------------------------------------------------------
subroutine backward_modeling_staggered42_is(is, par, coord, adj_src, c, den, fs, &
                                            nx_pml, nz_pml, npml, damp, ix1, ix2, iz1, iz2, store_snap)
  use io
  use datatype
  use math
  use source
  use string

  type(param), intent(in)       :: par
  type(acquisition), intent(in) :: coord
  logical, intent(in)           :: store_snap
  integer, intent(in)           :: is, fs(:), nx_pml, nz_pml, npml, ix1, ix2, iz1, iz2
  real, intent(in)              :: damp(:,:), c(:,:), den(:,:)
  real, intent(inout)           :: adj_src(:,:)

  ! Memory allocations
  allocate(u(nz_pml,nx_pml-1))
  allocate(w(nz_pml-1,nx_pml))
  allocate(p(nz_pml,nx_pml))
  allocate(p0(nz_pml,nx_pml))
  allocate(p_for(nz_pml,nx_pml))
  allocate(kernel(par%nz,par%nx))
  allocate(illum_forward(par%nz,par%nx))
  allocate(illum_backward(par%nz,par%nx))

  ! Initialize to zero
  u = 0.0
  w = 0.0
  p = 0.0
  p0 = 0.0
  p_for  = 0.0
  kernel = 0.0
  illum_forward  = 0.0
  illum_backward = 0.0

  ! Setup source and receiver position
  isx = npml+int(coord%xs(is)/par%dx)+1
  isz = npml+int(coord%zs(is)/par%dx)+1
  if (isz == fs(isx)) isz = isz + 1

  ! Time marching loop
  
  p0 = p
  
  do it=par%nt,2,-1
    ! Update the pressure field
    do ix=3,nx_pml-2
      do iz=fs(ix)+1,nz_pml-2
        alpha = den(iz,ix)*c(iz,ix)*c(iz,ix)*par%dt/(par%dx)
        kappa = damp(iz,ix)*par%dt
        p(iz,ix) = (1.0-kappa)*p(iz,ix) - alpha* &
                  (c1_staggered*(u(iz,ix)-u(iz,ix-1) + w(iz,ix)-w(iz-1,ix)) + &
                  c2_staggered*(u(iz,ix+1)-u(iz,ix-2) + w(iz+1,ix)-w(iz-2,ix)))
      enddo
    enddo

   ! Add source
    do ig=1,coord%ng(is)
      igx = npml + int(coord%xg(is,ig)/par%dx)+1
      igz = npml + int(coord%zg(is,ig)/par%dx)+1
      beta = par%dt
  
      ! Geophone is shifted downward 1 node if it is on the free surface.
      if (par%data(1:8) == 'pressure') then
        if (igz .gt. fs(igx)) then
          p(igz,igx) = p(igz,igx) + beta*adj_src(it,ig)
        else
          p(fs(igx)+1,igx) = p(fs(igx)+1,igx) + beta*adj_src(it,ig)
        endif
      else
        if (igz .gt. fs(igx)) then
          w(igz,igx) = w(igz,igx) + beta*adj_src(it,ig)
        else
          w(fs(igx)+1,igx) = w(fs(igx)+1,igx) + beta*adj_src(it,ig)
        endif
      endif
    enddo

    ! Free surface
    if (par%free_surface == 1) then
      p(fs-1,:) = -p(fs+1,:)
    endif

    ! Update horizontal particle velocity: u
    do ix=2,nx_pml-2
      do iz=fs(ix),nz_pml-1
        alpha = par%dt/(den(iz,ix)*par%dx)
        kappa = 0.5*(damp(iz,ix)+damp(iz,ix+1))*par%dt
        u(iz,ix) = (1.0-kappa)*u(iz,ix) - alpha*(c1_staggered*(p(iz,ix+1)-p(iz,ix)) &
                                                +c2_staggered*(p(iz,ix+2)-p(iz,ix-1)))
      enddo
    enddo

    ! Update verticle particle velocity: w
    do ix=2,nx_pml-1
      do iz=fs(ix),nz_pml-2
        alpha = par%dt/(den(iz,ix)*par%dx)
        kappa = 0.5*(damp(iz,ix)+damp(iz+1,ix))*par%dt
        w(iz,ix) = (1.0-kappa)*w(iz,ix) - alpha*(c1_staggered*(p(iz+1,ix)-p(iz,ix)) &
                                                +c2_staggered*(p(iz+2,ix)-p(iz-1,ix)))
      enddo
    enddo

    ! Free surface boundary condition for velocity field
    if (par%free_surface.eq.1) then
      do ix=1,nx_pml
        w(fs(ix)-1,ix) = w(fs(ix),ix)
      enddo
    endif
    
    if (mod(it,par%store_step) == 0) then 
      call filename(output, './data/syn/src', is, '_snapshot/')
      call filename(output, output, it, '_snapshot.bin')
      call read_binfile(output, p_for, nz_pml,nx_pml) 
      kernel = kernel + p_for(iz1:iz2,ix1:ix2)* (p0(iz1:iz2,ix1:ix2) - p(iz1:iz2,ix1:ix2))
      illum_forward  = illum_forward  + p_for(iz1:iz2,ix1:ix2) * p_for(iz1:iz2,ix1:ix2)
      illum_backward = illum_backward + p(iz1:iz2,ix1:ix2) * p(iz1:iz2,ix1:ix2)
    endif
    
    ! ! Save adjoint snapshot
    ! if (mod(it,par%store_step)==0 .and. store_snap) then
    !   call filename(output, par%data_out, is, '/')
    !   call filename(output, output, it, '_snapshot_adj.bin')
    !   call write_binfile(output, p, nz_pml, nx_pml)
    ! endif

  enddo  ! End of time-marching loop

  ! Write kernel.
  call filename(output, par%data_out, is, '_kernel_vp.bin')
  call write_binfile(output, kernel, par%nz, par%nx) 
  call filename(output, par%data_out, is, '_illum_forw.bin')
  call write_binfile(output, illum_forward,  par%nz, par%nx) 
  call filename(output, par%data_out, is, '_illum_back.bin')
  call write_binfile(output, illum_backward, par%nz, par%nx) 

  deallocate(u, w, p, p0, p_for, kernel, illum_forward, illum_backward)
  999 continue

end subroutine backward_modeling_staggered42_is


end module fdcore









