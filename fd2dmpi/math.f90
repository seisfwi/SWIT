! module containing mathematical utility functions/subroutines
!
module math

implicit none

interface interp1
  module procedure interp1_v1
  module procedure interp1_v2
end interface

interface dot
  module procedure dot_1d
  module procedure dot_2d
end interface

interface norm
  module procedure norm_vector
  module procedure norm_matrix
end interface

interface norm2
  module procedure norm2_vector
  module procedure norm2_matrix
  module procedure norm2_matrix_skip
end interface

interface normalize
  module procedure normalize_vector
  module procedure normalize_matrix
end interface

interface max_value
  module procedure max_value_2d
  module procedure max_value_2d_topo
end interface

interface min_value
  module procedure min_value_3arg
  module procedure min_value_4arg
  module procedure min_value_2d
  module procedure min_value_2d_topo
end interface

interface max_abs
  module procedure max_abs_1d
  module procedure max_abs_2d
end interface

interface min_abs
  module procedure min_abs_1d
  module procedure min_abs_2d
end interface

interface mean
  module procedure mean_1d
  module procedure mean_2d
end interface

interface mean_abs
  module procedure mean_abs_1d
  module procedure mean_abs_2d
end interface

interface smooth
  module procedure smooth_1d
  module procedure smooth_2d
  module procedure smooth_2d_wb
  module procedure smooth_2d_rect_filter
end interface

contains

!------------------------------------------------------------------------------
function norm_vector(x, n)

integer, intent(in) :: n
real, intent(in)    :: x(:)
integer             :: i
real                :: norm_vector

norm_vector = 0.0

do i=1,n
  norm_vector = norm_vector + x(i)*x(i)
enddo
norm_vector = sqrt(norm_vector)

end function norm_vector

!------------------------------------------------------------------------------
function norm2_vector(x, n)

integer, intent(in) :: n
real, intent(in)    :: x(:)
integer             :: i
real                :: temp, norm2_vector

norm2_vector = 0.0

do i=1,n
  temp = x(i)
  norm2_vector = norm2_vector + temp*temp
enddo

end function norm2_vector

!------------------------------------------------------------------------------
function norm_matrix(A, n1, n2)

integer, intent(in) :: n1, n2
real, intent(in)    :: A(:,:)
integer             :: i, j
real                :: norm_matrix

norm_matrix = 0.0

do i=1,n1
  do j=1,n2
    norm_matrix = norm_matrix + A(i,j)*A(i,j)
  enddo
enddo
norm_matrix = sqrt(norm_matrix)

end function norm_matrix

!------------------------------------------------------------------------------
function norm2_matrix(A, n1, n2)

integer, intent(in) :: n1, n2
real, intent(in)    :: A(:,:)
integer             :: i, j
real                :: norm2_matrix, temp

norm2_matrix = 0.0

do i=1,n1
  do j=1,n2
    temp = A(i,j)
    norm2_matrix = norm2_matrix + temp*temp
  enddo
enddo

end function norm2_matrix

!------------------------------------------------------------------------------
function norm2_matrix_skip(A, n1, n2, skip)

integer, intent(in) :: n1, n2, skip
real, intent(in)    :: A(:,:)
integer             :: i, j
real                :: norm2_matrix_skip, temp

norm2_matrix_skip = 0.0

do i=1,n1
  do j=1,n2,skip
    temp = A(i,j)
    norm2_matrix_skip = norm2_matrix_skip + temp*temp
  enddo
enddo

end function norm2_matrix_skip

!------------------------------------------------------------------------------
function min_abs_1d(x, n)

integer, intent(in) :: n
real, intent(in)    :: x(:)
integer             :: i
real                :: min_abs_1d

min_abs_1d = 10000000000.0

do i=1,n
  if (min_abs_1d > abs(x(i))) then
    min_abs_1d = abs(x(i))
  endif
enddo

end function min_abs_1d

!------------------------------------------------------------------------------
function min_abs_2d(x, n1, n2)

integer, intent(in) :: n1, n2
real, intent(in)    :: x(:,:)
integer             :: i, j
real                :: min_abs_2d

min_abs_2d = 100000000000000.0

do i=1,n1
  do j=1,n2
    if (min_abs_2d > abs(x(i,j))) then
      min_abs_2d = abs(x(i,j))
    endif
  enddo
enddo

end function min_abs_2d

!------------------------------------------------------------------------------
function max_abs_1d(x, n)

integer, intent(in) :: n
real, intent(in)    :: x(:)
integer             :: i
real                :: max_abs_1d

max_abs_1d = 0.0

do i=1,n
  if (max_abs_1d < abs(x(i))) then
    max_abs_1d = abs(x(i))
  endif
enddo

end function max_abs_1d

!------------------------------------------------------------------------------
function max_abs_2d(x, n1, n2)

integer, intent(in) :: n1, n2
real, intent(in)    :: x(:,:)
integer             :: i, j
real                :: max_abs_2d

max_abs_2d = 0.0

do i=1,n1
  do j=1,n2
    if (max_abs_2d < abs(x(i,j))) then
      max_abs_2d = abs(x(i,j))
    endif
  enddo
enddo

end function max_abs_2d

!------------------------------------------------------------------------------
function max_value_2d(x, n1, n2)

integer, intent(in) :: n1, n2
real, intent(in)    :: x(:,:)
integer             :: i, j
real                :: max_value_2d

max_value_2d = 0.0

do i=1,n1
  do j=1,n2
    if (max_value_2d < x(i,j)) then
      max_value_2d = x(i,j)
    endif
  enddo
enddo

end function max_value_2d

!------------------------------------------------------------------------------
function max_value_2d_topo(x, n1, n2, fs)

integer, intent(in) :: n1, n2, fs(:)
real, intent(in)    :: x(:,:)
integer             :: i, j
real                :: max_value_2d_topo

max_value_2d_topo = 0.0

do i=1,n2
  do j=1,n1
    if (j > fs(i)) then
      if (max_value_2d_topo < x(j,i)) then
        max_value_2d_topo = x(j,i)
      endif
    endif
  enddo
enddo

end function max_value_2d_topo

!------------------------------------------------------------------------------
function min_value_3arg(x1, x2, x3)

real, intent(in)    :: x1, x2, x3
real                :: min_value_3arg

min_value_3arg = x1
if (min_value_3arg > x2) min_value_3arg = x2
if (min_value_3arg > x3) min_value_3arg = x3

end function min_value_3arg

!------------------------------------------------------------------------------
function min_value_4arg(x1, x2, x3, x4)

real, intent(in)    :: x1, x2, x3, x4
real                :: min_value_4arg

min_value_4arg = x1
if (min_value_4arg > x2) min_value_4arg = x2
if (min_value_4arg > x3) min_value_4arg = x3
if (min_value_4arg > x4) min_value_4arg = x4

end function min_value_4arg

!------------------------------------------------------------------------------
function min_value_2d(x, n1, n2)

integer, intent(in) :: n1, n2
real, intent(in)    :: x(:,:)
integer             :: i, j
real                :: min_value_2d

min_value_2d = 10000000000.0

do i=1,n1
  do j=1,n2
    if (min_value_2d > x(i,j)) then
      min_value_2d = x(i,j)
    endif
  enddo
enddo

end function min_value_2d

!------------------------------------------------------------------------------
function min_value_2d_topo(x, n1, n2, fs)

integer, intent(in) :: n1, n2, fs(:)
real, intent(in)    :: x(:,:)
integer             :: i, j
real                :: min_value_2d_topo

min_value_2d_topo = 10000000000.0

do i=1,n1
  do j=1,n2
    if (i > fs(j)) then
      if (min_value_2d_topo > x(i,j)) then
        min_value_2d_topo = x(i,j)
      endif
    endif
  enddo
enddo

end function min_value_2d_topo

!------------------------------------------------------------------------------
function vmax_abs(x, n)

integer, intent(in) :: n
complex, intent(in) :: x(:)
integer             :: i
real                :: vmax_abs

vmax_abs = 0.0

do i=1,n
  if (vmax_abs < abs(x(i))) then
    vmax_abs = abs(x(i))
  endif
enddo

end function vmax_abs

!------------------------------------------------------------------------------
function maxvalue(x, n)

integer, intent(in) :: n
real, intent(inout) :: x(:)
integer             :: i
real                :: maxvalue

maxvalue = x(1)
do i=2,n
  if (maxvalue < x(i)) then
    maxvalue = x(i)
  endif
enddo

end function maxvalue

!------------------------------------------------------------------------------
subroutine normalize_vector(x, n)

integer, intent(in) :: n
real, intent(inout) :: x(:)
integer             :: i

x(1:n) = x(1:n)/norm(x,n)

end subroutine normalize_vector

!------------------------------------------------------------------------------
subroutine normalize_matrix(A, n1, n2)

integer, intent(in) :: n1, n2
real, intent(inout) :: A(:,:)
integer             :: i

A(1:n1,1:n2) = A(1:n1,1:n2)/norm(A,n1,n2)

end subroutine normalize_matrix

!------------------------------------------------------------------------------
subroutine smooth_1d(x, n, w)

integer, intent(in) :: n, w
real, intent(inout) :: x(:)
integer             :: i, j, k, half
real                :: factor, sumx

factor = 1.0/real(w)
half = (w-1)/2
do i=1,n
  k = 0
  sumx = 0.0
  do j=max(1,i-half),min(n,i+half)
    k = k+1
    sumx = sumx + x(j)
  enddo
  x(i) = sumx/real(k)
enddo

end subroutine smooth_1d

!------------------------------------------------------------------------------
! Square smoothing filter
subroutine smooth_2d(x, n1, n2, hw)

integer, intent(in) :: n1, n2, hw
real, intent(inout) :: x(:,:)
integer             :: i, j, k, i1, j1
real                :: sumx
real, allocatable   :: temp(:,:)

allocate(temp(n1,n2))
temp = x
do i=1,n1
  do j=1,n2
    k = 0
    sumx = 0.0
    do i1=max(1,i-hw),min(n1,i+hw)
      do j1=max(1,j-hw),min(n2,j+hw)
        k = k+1
        sumx = sumx + temp(i1,j1)
      enddo
    enddo
!    if (k > 0) then
      x(i,j) = sumx/real(k)
!    endif
  enddo
enddo
deallocate(temp)

end subroutine smooth_2d

!------------------------------------------------------------------------------
! Square smoothing filter
subroutine smooth_2d_pad(m, n1, n2, npad)

integer, intent(in) :: n1, n2, npad
real, intent(inout) :: m(:,:)
integer             :: i, j, i1, j1, npad2, n
real                :: sumx
real, allocatable   :: temp(:,:)

npad2 = 2*npad
n = (n1+npad2)*(n2+npad2)
allocate(temp(n1+npad2,n2+npad2))
call padding_2d(m, n1, n2, npad, temp)

do i=1,n1
  do j=1,n2
    sumx = 0.0
    do i1=i,i+npad2
      do j1=j,j+npad2
        sumx = sumx + temp(i1,j1)
      enddo
    enddo
    m(i,j) = sumx/real(n)
  enddo
enddo
deallocate(temp)

end subroutine smooth_2d_pad

!------------------------------------------------------------------------------
! Square smoothing filter
subroutine smooth_2d_wb(x, nx, nz, hw, izwb)

integer, intent(in) :: nx, nz, hw, izwb(:)
real, intent(inout) :: x(:,:)
integer             :: ix, iz, ixx, izz, k
real                :: sumx
real, allocatable   :: temp(:,:)

allocate(temp(nz,nx))
temp = x
do ix=1,nx
  do iz=1,nz
    if (iz >= izwb(ix)) then
      k = 0
      sumx = 0.0
      do ixx=max(1,ix-hw),min(nx,ix+hw)
        do izz=max(izwb(ixx),iz-hw),min(nz,iz+hw)
          k = k+1
          sumx = sumx + temp(izz,ixx)
        enddo
      enddo
      if (k > 0) then
        x(iz,ix) = sumx/real(k)
      endif
    endif
  enddo
enddo
deallocate(temp)

end subroutine smooth_2d_wb

!------------------------------------------------------------------------------
! Rectangular smoothing filter
subroutine smooth_2d_rect_filter(x, n1, n2, w1, w2)

integer, intent(in) :: n1, n2, w1, w2
real, intent(inout) :: x(:,:)
integer             :: i, j, k, i1, j1, half1, half2
real                :: sumx
real, allocatable   :: temp(:,:)

allocate(temp(n1,n2))
temp = x
half1 = int((w1-1)/2)
half2 = int((w2-1)/2)
do i=1,n1
  do j=1,n2
    k = 0
    sumx = 0.0
    do i1=max(1,i-half1),min(n1,i+half1)
      do j1=max(1,j-half2),min(n2,j+half2)
        k = k+1
        sumx = sumx + temp(i1,j1)
      enddo
    enddo
    if (k > 0) then
      x(i,j) = sumx/real(k)
    endif
  enddo
enddo
deallocate(temp)

end subroutine smooth_2d_rect_filter

!------------------------------------------------------------------------------
subroutine diff(x, n, dx)

integer, intent(in) :: n
real, intent(in)    :: dx
real, intent(inout) :: x(:)
integer             :: i

do i=1,n-1
  x(i) = (x(i+1)-x(i))/dx
enddo
x(n) = x(n-1)

end subroutine diff

!------------------------------------------------------------------------------
function dot_1d(x1, x2, n)

integer, intent(in) :: n
real, intent(in)    :: x1(:), x2(:)
real                :: dot_1d
integer             :: i

dot_1d = 0.0
do i=1,n
  dot_1d = dot_1d + x1(i)*x2(i)
enddo

end function dot_1d

!------------------------------------------------------------------------------
function dot_2d(x1, x2, n1, n2)

integer, intent(in) :: n1, n2
real, intent(in)    :: x1(:,:), x2(:,:)
real                :: dot_2d
integer             :: i, j

dot_2d = 0.0
do i=1,n2
  do j=1,n1
    dot_2d = dot_2d + x1(j,i)*x2(j,i)
  enddo
enddo

end function dot_2d

!------------------------------------------------------------------------------
function mean_1d(A, n)

integer, intent(in) :: n
real, intent(in)    :: A(:)
integer             :: i
real                :: mean_1d

mean_1d = 0.0
do i=1,n
  mean_1d = mean_1d + A(i)
enddo
mean_1d = mean_1d/real(n)

end function mean_1d

!------------------------------------------------------------------------------
function mean_2d(A, n1, n2)

integer, intent(in) :: n1, n2
real, intent(in)    :: A(:,:)
integer             :: i, j
real                :: mean_2d

mean_2d = 0.0
do i=1,n2
  do j=1,n1
    mean_2d = mean_2d + abs(A(j,i))
  enddo
enddo
mean_2d = mean_2d/real(n1*n2)

end function mean_2d

!------------------------------------------------------------------------------
function mean_abs_1d(A, n)

integer, intent(in) :: n
real, intent(in)    :: A(:)
integer             :: i
real                :: mean_abs_1d

mean_abs_1d = 0.0
do i=1,n
  mean_abs_1d = mean_abs_1d + abs(A(i))
enddo
mean_abs_1d = mean_abs_1d/real(n)

end function mean_abs_1d

!------------------------------------------------------------------------------
function mean_abs_2d(A, n1, n2)

integer, intent(in) :: n1, n2
real, intent(in)    :: A(:,:)
integer             :: i, j
real                :: mean_abs_2d

mean_abs_2d = 0.0
do i=1,n2
  do j=1,n1
    mean_abs_2d = mean_abs_2d + abs(A(j,i))
  enddo
enddo
mean_abs_2d = mean_abs_2d/real(n1*n2)

end function mean_abs_2d

!------------------------------------------------------------------------------
function quadratic_step(x2,x3,y1,y2,y3)

real, intent(in) :: x2, x3, y1, y2, y3
real             :: quadratic_step

quadratic_step = ((y1-y3)*x2*x2 + (y2-y1)*x3*x3)/((y1-y3)*x2 + (y2-y1)*x3)*0.5

end function quadratic_step

!------------------------------------------------------------------------------
function cubic_step(x2,x3,x4,y1,y2,y3,y4)

real, intent(in) :: x2, x3, x4, y1, y2, y3, y4
real             :: cubic_step, A(3,3), b(3), m, c(3)

A(1,1) = x2
A(1,2) = x2*x2
A(1,3) = x2*A(1,2)
A(2,1) = x3
A(2,2) = x3*x3
A(2,3) = x3*A(2,2)
A(3,1) = x4
A(3,2) = x4*x4
A(3,3) = x4*A(3,2)
b(1) = y2-y1
b(2) = y3-y1
b(3) = y4-y1

! Forward elimination
m = A(2,1)/A(1,1)
A(2,2) = A(2,2) - m*A(1,2)
A(2,3) = A(2,3) - m*A(1,3)
b(2) = b(2) - m*b(1)
m = A(3,1)/A(1,1)
A(3,2) = A(3,2) - m*A(1,2)
A(3,3) = A(3,3) - m*A(1,3)
b(3) = b(3) - m*b(1)
m = A(3,2)/A(2,2)
A(3,3) = A(3,3) - m*A(2,3)
b(3) = b(3) - m*b(2)

! Backward substitution
c(3) = b(3)/A(3,3)
c(2) = (b(2)-A(2,3)*c(3))/A(2,2)
c(1) = (b(1)-A(1,2)*c(2)-A(1,3)*c(3))/A(1,1)

cubic_step = c(2)*c(2) - 3.0*c(1)*c(3)
if (cubic_step < 0 .or. abs(c(3)) < 0.0000000001) then
  cubic_step = -1.0
else
  cubic_step = (-c(2)+sqrt(c(2)*c(2) - 3.0*c(1)*c(3)))/(3.0*c(3))
endif

end function cubic_step

!------------------------------------------------------------------------------
subroutine window(windowtype, n, w)

use global
character(len=*), intent(in) :: windowtype
integer, intent(in)          :: n
real, intent(out)            :: w(:)
integer                      :: i
real                         :: denom

if (windowtype(1:4) == 'hann') then
  denom = real(n-1)
  do i=1,n
    w(i) = 0.5*(1.0-cos(2.0*pi*real(i-1)/denom))
  enddo
endif

end subroutine window

!------------------------------------------------------------------------------
subroutine mute_iz(data,nx,nz,iz,taper)

integer, intent(in)    :: iz, nx, nz, taper
real   , intent(inout) :: data(:,:)
integer                :: ix, izmin, length, i
real, allocatable      :: w(:)

allocate(w(taper*2))

if (taper > nz) return

call window('hann',2*taper,w)
izmin = max(1,iz-taper+1)
length = iz - izmin+1

do ix=1,nx
  if (izmin > 1) then
    data(1:izmin-1,ix) = 0.0
  endif
  data(izmin:iz,ix) = data(izmin:iz,ix)*w(taper-length+1:taper)
enddo

deallocate(w)

end subroutine mute_iz

!------------------------------------------------------------------------------
subroutine mute(data,nx,nz,iz,taper)

integer, intent(in)    :: iz(:), nx, nz, taper
real   , intent(inout) :: data(:,:)
integer                :: ix, izmin, length, i
real, allocatable      :: w(:), wz(:)

allocate(w(taper*2))
allocate(wz(nz))

if (taper > nz) return

call window('hann',2*taper,w)

do ix=1,nx
  wz = 1.0
  izmin = max(1,iz(ix)-taper+1)
  length = iz(ix) - izmin+1
  if (izmin > 1) then
    data(1:izmin-1,ix) = 0.0
  endif
  data(izmin:iz(ix),ix) = data(izmin:iz(ix),ix)*w(taper-length+1:taper)
enddo

deallocate(w,wz)

end subroutine mute

!------------------------------------------------------------------------------
subroutine mute2(data,nx,nz,iz,taper)

integer, intent(in)    :: iz(:), nx, nz, taper
real   , intent(inout) :: data(:,:)
integer                :: ix, izmin, length, i, izix
real, allocatable      :: w(:)

allocate(w(taper*2))

if (taper > nz) return

call window('hann',2*taper,w)

do ix=1,nx
  izix = iz(ix)
  data(1:izix,ix) = 0.0
  data(izix+1:izix+taper,ix) = data(izix+1:izix+taper,ix)*w(1:taper)
enddo

deallocate(w)

end subroutine mute2

!------------------------------------------------------------------------------
subroutine mute_end(data,nx,nz,iz,taper)

integer, intent(in)    :: iz(:), nx, nz, taper
real   , intent(inout) :: data(:,:)
integer                :: ix, izmax, length, i
real, allocatable      :: w(:), wz(:)

allocate(w(taper*2))
allocate(wz(nz))

if (taper > nz) return

call window('hann',2*taper,w)

do ix=1,nx
  if (iz(ix) >1) then
    wz = 1.0
    izmax = min(nz,max(1,iz(ix)+taper-1))
    length = izmax - iz(ix) + 1
    if (izmax < nz) then
      data(izmax+1:nz,ix) = 0.0
    endif
    data(iz(ix):izmax,ix) = data(iz(ix):izmax,ix)*w(taper+1:taper+1+length)
  else
    data(:,ix) = 0.0
  endif
enddo

deallocate(w,wz)

end subroutine mute_end

!------------------------------------------------------------------------------
subroutine mute_1d(trace,n,imute,taper)

integer, intent(in)    :: imute, n, taper
real   , intent(inout) :: trace(:)
integer                :: ix, izmin, length, i
real, allocatable      :: w(:), wz(:)

allocate(w(taper*2))
allocate(wz(n))

if (taper > n) return

call window('hann',2*taper,w)

wz = 1.0
izmin = max(1,imute-taper+1)
length = imute - izmin+1
if (izmin > 1) then
  trace(1:izmin-1) = 0.0
endif
trace(izmin:imute) = trace(izmin:imute)*w(taper-length+1:taper)

deallocate(w,wz)

end subroutine mute_1d

!------------------------------------------------------------------------------
subroutine mute_1d_end(trace,n,imute,taper)

integer, intent(in)    :: imute, n, taper
real   , intent(inout) :: trace(:)
integer                :: ix, izmax, length, i
real, allocatable      :: w(:), wz(:)

allocate(w(taper*2))
allocate(wz(n))

if (taper > n) return

call window('hann',2*taper,w)

wz = 1.0
izmax = min(n,imute+taper-1)
if (izmax < n) then
  trace(izmax+1:n) = 0.0
endif
trace(imute-taper+1:imute) = trace(imute-taper+1:imute)*w(taper+1:2*taper)

deallocate(w,wz)

end subroutine mute_1d_end

!------------------------------------------------------------------------------
subroutine mute_1d_full(trace,n,taper)

integer, intent(in)    :: n, taper
real   , intent(inout) :: trace(:)
real, allocatable      :: w(:)

allocate(w(taper*2))

if (taper > n) return

call window('hann',2*taper,w)
trace(1:2*taper) = trace(1:2*taper)*w
trace(2*taper+1:n) = 0.0

deallocate(w)

end subroutine mute_1d_full

!------------------------------------------------------------------------------
! 1D Interpolation for uniform grids
subroutine interp1_v1(v_in,n1,d1,v_out,n2,d2)

integer, intent(in) :: n1, n2
real, intent(in)    :: v_in(:), d1, d2
real, intent(out)   :: v_out(:)
integer             :: i1, i2, ii1
real                :: x, h, alpha

do i2=1,n2
  x = real(i2-1)*d2
  i1 = int(x/d1)+1
  h = x - real(i1-1)*d1
  alpha = h/d1
  ii1 = i1+1
  if (ii1 <= n1) then
    v_out(i2) = (1.0-alpha)*v_in(i1) + alpha*v_in(i1+1)
  else
    v_out(i2) = v_in(i1)
  endif
enddo

end subroutine interp1_v1

!---------------------------------------------------------------
! 1D Interpolation for irregular grids
subroutine interp1_v2(x1,y1,n1,x2,y2,n2)

integer, intent(in) :: n1, n2
real, intent(in)    :: x1(:), y1(:), x2(:)
real, intent(out)   :: y2(:)
integer             :: i1, i2, ii1, i1p1
real                :: x, h, alpha

do i2=1,n2
  x = x2(i2)
  if (x <= x1(1)) then
    y2(i2) = y1(1)
    goto 200
  elseif (x >= x1(n1)) then
    y2(i2) = y1(n1)
    goto 200
  endif
  do i1=1,n1-1
    ii1 = i1
    if (x >= x1(i1) .and. x < x1(i1+1)) then
      goto 100
    endif
  enddo
  100 continue
  h = x - x1(i1)
  i1p1 = i1+1
  alpha = h/(x1(i1p1)-x1(i1))
  y2(i2) = (1.0-alpha)*y1(i1) + alpha*y1(i1p1)
  200 continue
enddo

end subroutine interp1_v2

!---------------------------------------------------------------
! 2D Interpolation for uniform grids
subroutine interp2(v_in,nx_in,nz_in,dx_in,v_out,nx_out,nz_out,dx_out)

integer, intent(in)  :: nx_in, nz_in, nx_out, nz_out
real, intent(in)     :: v_in(:,:), dx_in, dx_out
real, intent(out)    :: v_out(:,:)
integer              :: ix, iz, ix1, iz1, ixp1, izp1
real                 :: x, z, hx, hz, alpha_x, alpha_z

do ix=1,nx_out
  x = real(ix-1)*dx_out
  ix1 = int(x/dx_in)+1
  ixp1 = ix1+1
  hx = x - real(ix1-1)*dx_in
  alpha_x = hx/dx_in
  do iz=1,nz_out
    z = real(iz-1)*dx_out
    iz1 = int(z/dx_in)+1
    izp1 = iz1+1
    hz = z - real(iz1-1)*dx_in
    alpha_z = hz/dx_in

    if (ixp1 <= nx_in .and. izp1 <= nz_in) then
      v_out(iz,ix) = (1.0-alpha_z) * ((1.0-alpha_x)*v_in(iz1  ,ix1) + alpha_x*v_in(iz1  ,ix1+1)) + &
                          alpha_z * ((1.0-alpha_x)*v_in(iz1+1,ix1) + alpha_x*v_in(iz1+1,ix1+1))
    elseif (ixp1 <= nx_in .and. izp1 > nz_in) then
      v_out(iz,ix) = (1.0-alpha_x) * v_in(iz1  ,ix1) + alpha_x*v_in(iz1  ,ix1+1)
    elseif (ixp1 > nx_in .and. izp1 <= nz_in) then
      v_out(iz,ix) = (1.0-alpha_z) * v_in(iz1  ,ix1) + alpha_z * v_in(iz1+1,ix1)
    else
      v_out(iz,ix) = v_in(iz1  ,ix1)
    endif
  enddo
enddo

end subroutine interp2

!------------------------------------------------------------------------------
subroutine damp_gradient(g,nx,nz,npml,fs,taper)

integer, intent(in)    :: nx, nz, npml, fs(:), taper
real   , intent(inout) :: g(:,:)
integer                :: ix, izmin, length, ifs
real, allocatable      :: w(:)

allocate(w(taper*2))

if (taper > nz) return

call window('hann',2*taper,w)

do ix=1,nx
  ifs = fs(ix+npml)-npml
  if (ifs > 0) g(1:ifs,ix) = 0.0
  g(ifs+1:ifs+taper,ix) = g(ifs+1:ifs+taper,ix)*w(1:taper)
enddo

deallocate(w)

end subroutine damp_gradient

!------------------------------------------------------------------------------------
! Padding the input array m_in
subroutine padding_2d(m_in,n1,n2,npad,m_out)

integer,    intent(in)  :: n1, n2, npad
real,       intent(in)  :: m_in(:,:)
real,       intent(out) :: m_out(:,:)
integer                 :: i

m_out(npad+1:npad+n1,npad+1:npad+n2) = m_in
do i=1,npad
  m_out(i,:) = m_in(1,:)
  m_out(npad+n1+i,:) = m_in(n1,:)
enddo
do i=1,npad
  m_out(:,i) = m_in(:,1)
  m_out(:,npad+n2+i) = m_in(:,n2)
enddo

end subroutine padding_2d

end module math

