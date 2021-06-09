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
!时间二阶、空间四阶同位网格差分系数
real, parameter :: c1_2nd_order = -2.5
real, parameter :: c2_2nd_order = 4.0/3.0
real, parameter :: c3_2nd_order = -1.0/12.0

!时间二阶、空间四阶交错网格差分系数
real, parameter :: c1_staggered = 9.0/8.0
real, parameter :: c2_staggered = -1.0/24.0

real, parameter :: prefact = 1.0/180.0
real, parameter :: c1_6th = -490
real, parameter :: c2_6th = 270.0
real, parameter :: c3_6th = -27
real, parameter :: c4_6th = 2
end module global

