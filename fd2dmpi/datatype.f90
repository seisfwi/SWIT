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


! module for defining data structures
!
module datatype

type param
  integer            :: nx, nz, nt, free_surface, npml, store_snap, store_step
  real               :: dx, dt, f, vmin, vmax, xmin, xmax
  character(len=200) :: coordfile, data_out, velfile, densityfile, fileformat, sourcefile, data, sourcetype
                        
end type param

type acquisition
  ! Data format: t(is,ig), xg(iz,ix), zg(iz,ix)

  integer :: ns, ngmax
  integer, allocatable :: ng(:)
  real, allocatable :: xs(:), zs(:), xg(:,:), zg(:,:), t(:,:)
end type acquisition

type wavelet_param
  character(len=200) :: fileformat, inputprefix, outputprefix
  integer            :: ns, ng, nt, iref
  real               :: dx, dt, ds, dg
end type wavelet_param

end module datatype

