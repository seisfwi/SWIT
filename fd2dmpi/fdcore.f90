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


! module for 2D forward modeling using a time-domain FD method
!
module fdcore

use global
implicit none

integer, private              :: ix, iz, it, is, is1, is2, isx, isz, ig, igx, igz, i, j, jj
real,    private              :: alpha, beta, a, d1, d2, kappa
real,    private, allocatable :: u(:,:), w(:,:), p(:,:), csg(:,:), &
                                 p_end(:,:), u_end(:,:), w_end(:,:), boundary(:,:)


contains


!-----------------------------------------------------------------------------------------
subroutine staggered42_modeling_is(is, par, coord, s, c, den, fs, nx_pml, nz_pml, npml, damp, store_boundary)

  use io
  use datatype
  use math
  use source
  use string

  type(param), intent(in)       :: par
  type(acquisition), intent(in) :: coord
  logical, intent(in)           :: store_boundary
  integer, intent(in)           :: is, fs(:), nx_pml, nz_pml, npml
  real, intent(in)              :: damp(:,:), c(:,:), den(:,:)
  real, intent(inout)           :: s(:)

  ! Memory allocations
  allocate(u(nz_pml,nx_pml-1))
  allocate(w(nz_pml-1,nx_pml))
  allocate(p(nz_pml,nx_pml))
  allocate(csg(par%nt,coord%ngmax))

  if (store_boundary) then 
    allocate(boundary(par%nx*9+par%nz*6,par%nt))
    allocate(p_end(par%nz,par%nx))
    allocate(u_end(par%nz,par%nx))
    allocate(w_end(par%nz,par%nx))
  endif

  u = 0.0
  w = 0.0
  p = 0.0
  csg = 0.0
  if (store_boundary) then 
    boundary = 0.0
    p_end = 0.0
    u_end = 0.0
    w_end = 0.0
  endif

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

    if (store_boundary) then
      jj = 0
      ! Bottom wavefields
      do j=npml+1,nx_pml-npml
        boundary(jj+1,it) = p(nz_pml-npml  ,j)
        boundary(jj+2,it) = p(nz_pml-npml-1,j)
        boundary(jj+3,it) = w(nz_pml-npml-1,j)
        jj = jj+3
      enddo
      ! Left and right wavefields
      do j=npml+1,nz_pml-npml
        boundary(jj+1,it) = p(j,npml+1)
        boundary(jj+2,it) = p(j,npml+2)
        boundary(jj+3,it) = u(j,npml+1)
        boundary(jj+4,it) = p(j,nx_pml-npml)
        boundary(jj+5,it) = p(j,nx_pml-npml-1)
        boundary(jj+6,it) = u(j,nx_pml-npml-1)
        jj = jj+6
      enddo
      ! Top wavefields
      do j=npml+1,nx_pml-npml
        boundary(jj+1,it) = p(fs(j)  ,j)
        boundary(jj+2,it) = p(fs(j)+1,j)
        boundary(jj+3,it) = p(fs(j)+2,j)
        boundary(jj+4,it) = u(fs(j)  ,j)
        boundary(jj+5,it) = u(fs(j)+1,j)
        boundary(jj+6,it) = w(fs(j)  ,j)
        jj = jj+6
      enddo
    endif
  enddo  ! End of time-marching loop

  if (store_boundary) then

    p_end = p(npml+1:nz_pml-npml,npml+1:nx_pml-npml)
    u_end = u(npml+1:nz_pml-npml,npml+1:nx_pml-npml)
    w_end = w(npml+1:nz_pml-npml,npml+1:nx_pml-npml)
    
    ! save boundary, p_end, u_end, w_end
    call filename(output, par%data_out, is, '_snapshot/boundary.bin')
    call write_binfile(output, boundary, par%nx*9+par%nz*6, par%nt)
    call filename(output, par%data_out, is, '_snapshot/p_end.bin')
    call write_binfile(output, p_end, par%nz, par%nx)
    call filename(output, par%data_out, is, '_snapshot/u_end.bin')
    call write_binfile(output, u_end, par%nz, par%nx)
    call filename(output, par%data_out, is, '_snapshot/w_end.bin')
    call write_binfile(output, w_end, par%nz, par%nx)

  endif

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
  if (store_boundary) then
    deallocate(boundary, p_end, u_end, w_end)
  endif 
  
  999 continue

end subroutine staggered42_modeling_is


!-----------------------------------------------------------------------------------------
subroutine staggered42_reco_it(is, it, par, coord, s, c, den, fs, boundary, p, u, w)

  use io
  use datatype
  use pml
  use string
  
  type(param), intent(in)       :: par
  type(acquisition), intent(in) :: coord
  integer, intent(in)           :: fs(:), is, it
  real, intent(in)              :: s(:), c(:,:), den(:,:)
  real, intent(inout)           :: p(:,:), u(:,:), w(:,:)
  real, intent(in)              :: boundary(:,:)
  
  ! Setup source and receiver position
  isx = int(coord%xs(is)/par%dx)+1
  isz = int(coord%zs(is)/par%dx)+1
  
  ! Update horizontal particle velocity: u
  do ix=2,par%nx-2
    do iz=fs(ix+npml)-npml+1,par%nz-1
      alpha = par%dt/(den(iz,ix)*par%dx)
      u(iz,ix) = u(iz,ix) + alpha*(c1_staggered*(p(iz,ix+1)-p(iz,ix)) &
                                  +c2_staggered*(p(iz,ix+2)-p(iz,ix-1)))
    enddo
  enddo
  ! Update verticle particle velocity: w
  do ix=2,par%nx-1
    do iz=fs(ix+npml)-npml+1,par%nz-2
      alpha = par%dt/(den(iz,ix)*par%dx)
      w(iz,ix) = w(iz,ix) + alpha*(c1_staggered*(p(iz+1,ix)-p(iz,ix)) &
                                  +c2_staggered*(p(iz+2,ix)-p(iz-1,ix)))
    enddo
  enddo
  
  ! Insert boundary values
  jj = 0
  do j=1,par%nx
    w(par%nz-1,j) = boundary(jj+3,it-1)
    jj = jj+3
  enddo
  do j=1,par%nz
    u(j,1       ) = boundary(jj+3,it-1)
    u(j,par%nx-1) = boundary(jj+6,it-1)
    jj = jj+6
  enddo
  
  ! Subtract source
  if (isx > 2 .and. isx < par%nx-1 .and. isz > 3 .and. isz < par%nz-1) then
    beta = par%dt
    if (isz .gt. fs(isx+npml)-npml) then
      p(isz,isx) = p(isz,isx) - beta*s(it+1)
    else
      p(fs(isx+npml)-npml+1,isx) = p(fs(isx+npml)-npml+1,isx) - beta*s(it+1)
    endif
  endif
  
  ! Update the pressure field
  do ix=3,par%nx-2
    do iz=fs(ix+npml)-npml+2,par%nz-2
      alpha = den(iz,ix)*c(iz,ix)*c(iz,ix)*par%dt/(par%dx)
      p(iz,ix) = p(iz,ix) + alpha* &
                (c1_staggered*(u(iz,ix)-u(iz,ix-1) + w(iz,ix)-w(iz-1,ix)) + &
                 c2_staggered*(u(iz,ix+1)-u(iz,ix-2) + w(iz+1,ix)-w(iz-2,ix)))
    enddo
  enddo
  
  jj=0
  do j=1,par%nx
    p(par%nz  ,j) = boundary(jj+1,it-1)
    p(par%nz-1,j) = boundary(jj+2,it-1)
    jj=jj+3
  enddo
  do j=1,par%nz
    p(j,1)        = boundary(jj+1,it-1)
    p(j,2)        = boundary(jj+2,it-1)
    p(j,par%nx  ) = boundary(jj+4,it-1)
    p(j,par%nx-1) = boundary(jj+5,it-1)
    jj=jj+6
  enddo
  do j=1,par%nx
    p(fs(j+npml)-npml  ,j) = boundary(jj+1,it-1)
    p(fs(j+npml)-npml+1,j) = boundary(jj+2,it-1)
    p(fs(j+npml)-npml+2,j) = boundary(jj+3,it-1)
    jj=jj+6
  enddo
  
end subroutine staggered42_reco_it
  

!-----------------------------------------------------------------------------------------
subroutine staggered42_back_it(is, it, par, coord, s, c, den, fs, nx_pml, nz_pml, npml, damp, p, u, w)
  
  use io
  use datatype
  use math
  
  type(param), intent(in)       :: par
  type(acquisition), intent(in) :: coord
  integer, intent(in)           :: fs(:), nx_pml, nz_pml, npml, is, it
  real, intent(in)              :: damp(:,:), c(:,:), s(:,:), den(:,:)
  real, intent(inout)           :: p(:,:), u(:,:), w(:,:)
  
  ! Setup source and receiver position
  isx = npml+int(coord%xs(is)/par%dx)+1
  isz = npml+int(coord%zs(is)/par%dx)+1
  
  ! Update horizontal particle velocity: u
  do ix=2,nx_pml-2
    do iz=2,nz_pml-1
      alpha = par%dt/(den(iz,ix)*par%dx)
      kappa = 0.5*(damp(iz,ix)+damp(iz,ix+1))*par%dt
      u(iz,ix) = (1.0-kappa)*u(iz,ix) - alpha*(c1_staggered*(p(iz,ix+1)-p(iz,ix)) &
                                                +c2_staggered*(p(iz,ix+2)-p(iz,ix-1)))
    enddo
  enddo
  ! Update verticle particle velocity: w
  do ix=2,nx_pml-1
    do iz=2,nz_pml-2
      alpha = par%dt/(den(iz,ix)*par%dx)
      kappa = 0.5*(damp(iz,ix)+damp(iz+1,ix))*par%dt
      w(iz,ix) = (1.0-kappa)*w(iz,ix) - alpha*(c1_staggered*(p(iz+1,ix)-p(iz,ix)) &
                                                +c2_staggered*(p(iz+2,ix)-p(iz-1,ix)))
    enddo
  enddo
  
  ! Free surface boundary condition for velocity field
  if (par%free_surface.eq.1) then
    do ix=2,nx_pml-1
      w(fs(ix)-1,ix) = w(fs(ix),ix)
    enddo
  endif
  
  ! Update the pressure field
  do ix=3,nx_pml-2
    do iz=3,nz_pml-2
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
    ! if (par%data(1:8) == 'pressure') then
      if (igz .gt. fs(igx)) then
        p(igz,igx) = p(igz,igx) + beta*s(it,ig)
      else
        p(fs(igx)+1,igx) = p(fs(igx)+1,igx) + beta*s(it,ig)
      endif
    ! else
    !   if (igz .gt. fs(igx)) then
    !     w(igz,igx) = w(igz,igx) + beta*s(it,ig)
    !   else
    !     w(fs(igx)+1,igx) = w(fs(igx)+1,igx) + beta*s(it,ig)
    !   endif
    ! endif
  enddo
  
  ! Free surface
  do ix=1,nx_pml
  !  p(fs(ix),ix) = 0.0
    p(fs(ix)-1,ix) = -p(fs(ix)+1,ix)
  enddo
  
end subroutine staggered42_back_it
 

end module fdcore









