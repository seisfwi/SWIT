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



! module for defining PML arrays
!
module pml

implicit none
integer :: npml, nx_pml, nz_pml, iz1, iz2, ix1, ix2
real, allocatable :: damp(:), damp_global(:,:)

contains

subroutine init_pml(nx,nz,n_pml)

integer, intent(in) :: nx, nz, n_pml

npml = n_pml
nx_pml = nx+2*npml
nz_pml = nz+2*npml
ix1 = npml+1
ix2 = npml+nx
iz1 = npml+1
iz2 = npml+nz
allocate(damp(npml))
allocate(damp_global(nz_pml,nx_pml))

end subroutine init_pml

!-------------------------------------------------------------------
subroutine setup_damp_compact(nx_compact, nz_compact, dx, vmin, damp_compact)

integer, intent(in) :: nx_compact, nz_compact
real, intent(in)    :: dx, vmin
real, intent(out)   :: damp_compact(:,:)
real                :: a, xa, kappa
integer             :: ix, iz, nx_compact_pml, nz_compact_pml

nx_compact_pml = nx_compact + 2*npml
nz_compact_pml = nz_compact + 2*npml
damp_compact = 0.0
a = (npml-1)*dx
kappa = 15.0*vmin/(2.0*a)
do ix=1,npml
  xa = real(ix-1)*dx/a
  damp(ix) = kappa*xa*xa
enddo
do ix=1,npml
  do iz=1,nz_compact_pml
    damp_compact(iz,npml-ix+1) = damp(ix)
    damp_compact(iz,nx_compact_pml+ix-1-npml) = damp(ix)
  enddo
enddo
do iz=1,npml
  do ix=1+(npml-(iz-1)),nx_compact_pml-(npml-(iz-1))
    damp_compact(npml-iz+1,ix) = damp(iz)
    damp_compact(nz_compact_pml+iz-1-npml,ix) = damp(iz)
  enddo
enddo

end subroutine setup_damp_compact

!-------------------------------------------------------------------
subroutine setup_pml(dx, vmin)

real, intent(in) :: dx, vmin
real             :: a, xa, kappa
integer          :: ix, iz


damp_global = 0.0
a = (npml-1)*dx

! Set the damping factor kappa
if (vmin .eq. 0.0) then
  kappa = 3.0*1000.*log(1000.0)/(2.0*a)  ! Adjust the damping effect.
else
  kappa = 3.0*vmin*log(1000.0)/(2.0*a)   ! Adjust the damping effect.
endif

do ix=1,npml
  xa = real(ix-1)*dx/a
  damp(ix) = kappa*xa*xa
enddo
do ix=1,npml
  do iz=1,nz_pml
    damp_global(iz,npml-ix+1) = damp(ix)
    damp_global(iz,nx_pml+ix-npml) = damp(ix)
  enddo
enddo
do iz=1,npml
  do ix=1+(npml-(iz-1)),nx_pml-(npml-(iz-1))
    damp_global(npml-iz+1,ix) = damp(iz)
    damp_global(nz_pml+iz-npml,ix) = damp(iz)
  enddo
enddo

end subroutine setup_pml

!------------------------------------------------------------------------------------
! m is the model to be padded
subroutine padmodel(par,m,npml,nx_pml,nz_pml)

use datatype
!use global

type(param),intent(in)  :: par
integer,    intent(in)  :: npml, nx_pml, nz_pml
real,       intent(out) :: m(:,:)
integer :: ix, iz

! Extrapolate velocity in PML regions
do ix=1,npml
  m(npml+1:npml+par%nz,ix) = m(npml+1:npml+par%nz,npml+1)
  m(npml+1:npml+par%nz,nx_pml-npml+ix) = m(npml+1:npml+par%nz,nx_pml-npml)
enddo
do iz=1,npml
  m(iz,:) = m(npml+1,:)
  m(nz_pml-npml+iz,:) = m(nz_pml-npml,:)
enddo

end subroutine padmodel

end module pml

