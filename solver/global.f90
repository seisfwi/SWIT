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


! module defining global variables
!
module global

implicit none

! Constants
real(4), parameter :: pi        = 3.14159265358979
real(4), parameter :: twopi     = 2.0*pi

! Type of job
character(len=200) :: jobtype, input, output

! Record length unit
integer, parameter :: i4 = 4

! FD coefficients
real, parameter :: c1_2nd_order = -2.5
real, parameter :: c2_2nd_order = 4.0/3.0
real, parameter :: c3_2nd_order = -1.0/12.0

real, parameter :: c1_staggered = 9.0/8.0
real, parameter :: c2_staggered = -1.0/24.0

real, parameter :: prefact = 1.0/180.0
real, parameter :: c1_6th = -490
real, parameter :: c2_6th = 270.0
real, parameter :: c3_6th = -27
real, parameter :: c4_6th = 2
end module global

