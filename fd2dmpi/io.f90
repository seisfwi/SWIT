! module containing I/O utilities
!
module io

use string
implicit none

interface read_sufile
  module procedure read_sufile1
  module procedure read_sufile_2d
end interface

interface read_binfile
  module procedure read_binfile_1d
  module procedure read_binfile_2d
  module procedure read_binfile_3d
end interface

interface write_sufile
  module procedure write_sufile_1d
  module procedure write_sufile_2d
  module procedure write_sufile_2d_coord
end interface

interface write_binfile
  module procedure write_binfile_1d
  module procedure write_binfile_2d
  module procedure write_binfile_3d
end interface

contains

!---------------------------------------------------------------------------------------------
! Rules for verifying the input parameters
!
! STRING:
!         - required parameter will have a default as 'n/a'
!         - optional parameter will have some default which is not 'n/a'
!
! INTEGER:
!         - required parameter will have a default as -1
!         - optional parameter will have some default which is not -1
!
! FLOATING-POINT:
!         - required parameter will have a default as -1.0
!         - optional parameter will have some default which is not -1.0
!
subroutine readparamfile(parfile, par)

use datatype
use parser
use mmi_mpi

character(len=*), intent(in) :: parfile
type(param), intent(out)     :: par
character(len=100)           :: tmp, velfile, coordfile
integer                      :: nx, nz, nt
real                         :: dx, dt, f

if (rank == 0) then
  ! Required String parameters
  call readParFile(parfile, 'COORD_FILE',       par%coordfile,          'n/a') !par%coordfile观测系统文件名
  call readParFile(parfile, 'DATA_OUT',         par%data_out,           'src') !par%data_out理论合成的地震记录文件名
  call readParFile(parfile, 'VEL_IN',           par%velfile,            'n/a') !par%velfile速度模型文件名
  call readParFile(parfile, 'DENSITYFILE',      par%densityfile,        'n/a') !par%densityfile密度模型文件名
  call readParFile(parfile, 'FILEFORMAT',       par%fileformat,         'su')  !par%fileformat文件格式，缺省值为su
  call readParFile(parfile, 'DATA',             par%data,               'pressure')
  call readParFile(parfile, 'SOURCETYPE',       par%sourcetype,         'normal')
  ! Required Integer parameters
  call readParFile(parfile, 'NX',               par%nx,                 -1)
  call readParFile(parfile, 'NZ',               par%nz,                 -1)
  call readParFile(parfile, 'NT_WORK',          par%nt,                 -1)
  call readParFile(parfile, 'FREESURFACE',      par%free_surface,        1)
  call readParFile(parfile, 'NPML',             par%npml,               40)
  call readParFile(parfile, 'STORE_SNAP',       par%store_snap,          0)
  call readParFile(parfile, 'STORE_STEP',       par%store_step,          5)
  ! Floating-point parameters
  call readParFile(parfile, 'DX',               par%dx,               -1.0)
  call readParFile(parfile, 'DT_WORK',          par%dt,               -1.0)
  call readParFile(parfile, 'FREQUENCY',        par%f,                 5.0)
  call readParFile(parfile, 'VMIN',             par%vmin,            330.0)
  call readParFile(parfile, 'VMAX',             par%vmax,           7000.0)
  call readParFile(parfile, 'XMIN',             par%xmin,              0.0)
  call readParFile(parfile, 'XMAX',             par%xmax,  real(par%nx-1)*par%dx)
endif

call MPI_BCAST(par%data_out,100,     MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%velfile,100,      MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%densityfile,100,  MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%fileformat,100,   MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%data,100,         MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%sourcetype,100,   MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%nx,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%nz,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%nt,1,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%free_surface,1,   MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%npml,1,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%store_snap,1,     MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%store_step,1,     MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

call MPI_BCAST(par%dx,1,             MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%dt,1,             MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%f,1,              MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%vmin,1,           MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%vmax,1,           MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%xmin,1,           MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(par%xmax,1,           MPI_REAL,0,MPI_COMM_WORLD,ierr)


end subroutine readparamfile

!-------------------------------------------------------------------------------
subroutine readcoordfile(coordfile, coord)

use datatype
use mmi_mpi

character(len=*), intent(in)   :: coordfile
type(acquisition), intent(out) :: coord
integer                        :: i, j, is0, is, ig, ns, ng, ngmax
real                           :: xs, zs, xg, zg, t

if (rank == 0) then
  open(10,file=coordfile,form='formatted')

  ! Determine the number of sources
  ns = 0
  read(10,*,end=100) is0, ig, xs, zs, xg, zg, t
  ns = 1
  do i=2,1000000000
    read(10,*,end=100) is, ig, xs, zs, xg, zg, t
    if (is.ne.is0) then
      is0 = is
      ns = ns + 1
    endif
  enddo
  100 continue
  rewind(10)
  coord%ns = ns
  allocate(coord%ng(ns))
  allocate(coord%xs(ns))
  allocate(coord%zs(ns))

  ! Determine the number of geophones for each shot
  read(10,*,end=200) is0, ig, xs, zs, xg, zg, t
  ns = 1
  ng = 1
  do i=2,1000000000
    read(10,*,end=200) is, ig, xs, zs, xg, zg, t
    if (is.ne.is0) then
      is0 = is
      coord%ng(ns) = ng
!    write(*,*) 'ng(',ns,') = ', coord%ng(ns)
      ns = ns + 1
      ng = 1
    else
      ng = ng + 1
    endif
  enddo
  200 continue
  coord%ng(ns) = ng
  !write(*,*) 'ng(',ns,') = ', coord%ng(ns)
  rewind(10)

  ! Determine the maximum number of geophones per shot
  coord%ngmax = 0
  do is=1,coord%ns
    if (coord%ngmax < coord%ng(is)) coord%ngmax = coord%ng(is)
  enddo
  !write(*,*) 'ng max = ', coord%ngmax

  allocate(coord%xg(is,coord%ngmax))
  allocate(coord%zg(is,coord%ngmax))
  allocate(coord%t(is,coord%ngmax))!每一炮的检波器个数

  ! Read source and receiver positions
  do i=1,coord%ns
    do j=1,coord%ng(i)
      read(10,*,end=300) is, ig, xs, zs, xg, zg, t
      coord%xs(i) = xs
      coord%zs(i) = zs
      coord%xg(i,j) = xg
      coord%zg(i,j) = zg
      coord%t(i,j) = t
    enddo
  enddo
  300 continue
  close(10)

endif

call MPI_BCAST(coord%ns,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (rank > 0) then
  allocate(coord%ng(coord%ns))
  allocate(coord%xs(coord%ns))
  allocate(coord%zs(coord%ns))
endif
call MPI_BCAST(coord%ng,coord%ns,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(coord%xs,coord%ns,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(coord%zs,coord%ns,MPI_REAL,0,MPI_COMM_WORLD,ierr)
call MPI_BCAST(coord%ngmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (rank > 0) then
  allocate(coord%xg(coord%ns,coord%ngmax))
  allocate(coord%zg(coord%ns,coord%ngmax))
  allocate(coord%t(coord%ns,coord%ngmax))
endif
do is=1,coord%ns
  call MPI_BCAST(coord%xg(is,:),coord%ng(is),MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(coord%zg(is,:),coord%ng(is),MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(coord%t(is,:),coord%ng(is),MPI_REAL,0,MPI_COMM_WORLD,ierr)
enddo

999 continue

end subroutine readcoordfile

!------------------------------------------------------------------------------------
subroutine readvelfile(par,c,npml,nx_pml,nz_pml)

use datatype
use global
use mmi_mpi

type(param),intent(in) :: par
real, intent(out) :: c(:,:)
integer, intent(in) :: npml, nx_pml, nz_pml
integer :: ix, iz

if (rank == 0) then
  ! Read velocity model
  open(10,file=par%velfile,access='direct',recl=i4*par%nz)
  do ix=1,par%nx
    read(10,rec=ix) c(npml+1:npml+par%nz,ix+npml)
  enddo
  close(10)

  ! Extrapolate velocity in PML regions
  do ix=1,npml
    c(npml+1:npml+par%nz,ix) = c(npml+1:npml+par%nz,npml+1)
    c(npml+1:npml+par%nz,nx_pml-npml+ix) = c(npml+1:npml+par%nz,nx_pml-npml)
  enddo
  do iz=1,npml
    c(iz,:) = c(npml+1,:)
    c(nz_pml-npml+iz,:) = c(nz_pml-npml,:)
  enddo
endif
do ix=1,nx_pml
  call MPI_BCAST(c(:,ix),nz_pml,MPI_REAL,0,MPI_COMM_WORLD,ierr)
enddo

end subroutine readvelfile


!------------------------------------------------------------------------------------
subroutine read_densityfile(par,c,den,npml,nx_pml,nz_pml)

use datatype
use global
use mmi_mpi

type(param),intent(in) :: par
real, intent(in)       :: c(:,:)
real, intent(out)      :: den(:,:)
integer, intent(in)    :: npml, nx_pml, nz_pml
integer                :: ix, iz


if (rank == 0) then
  ! Read density model
  call read_binfile_2d(par%densityfile,den(npml+1:npml+par%nz,npml+1:npml+par%nx),par%nz,par%nx)
  ! Extrapolate density in PML regions
  do ix=1,npml
    den(npml+1:npml+par%nz,ix) = den(npml+1:npml+par%nz,npml+1)
    den(npml+1:npml+par%nz,nx_pml-npml+ix) = den(npml+1:npml+par%nz,nx_pml-npml)
  enddo
  do iz=1,npml
    den(iz,:) = den(npml+1,:)
    den(nz_pml-npml+iz,:) = den(nz_pml-npml,:)
  enddo
endif

do ix=1,nx_pml
  call MPI_BCAST(den(:,ix),nz_pml,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
enddo


end subroutine read_densityfile


!------------------------------------------------------------------------------
subroutine read_sufile1(sufile, trace, nt, ntrace, sampling_interval)

use global
use su

character(len=*), intent(in) :: sufile
integer, intent(in) :: nt, ntrace
real, intent(out) :: sampling_interval
real, intent(out) :: trace(:,:)
integer :: nlen, i

nlen = nhead+nt
open(10,file=sufile,access='direct',recl=i4*nlen)
do i=1,ntrace
  read(10,rec=i,err=111) &
        tracl,tracr,fldr,tracf,ep,cdp,cdpt, &
        trid,nvs,nhs,duse,offset,gelev,selev,sdepth, &
        gdel,sdel,swdep,gwdep,scalel,scalco, &
        sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat, &
        tstat,laga,lagb,delrt,muts,mute,ns,dt, &
        gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
        stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
        lcs,hcs,year,day,hour,minute,sec,timbas, &
        trwf,grnors,grnofr,grnlof,gaps,otrav, &
        d1,f1,d2,f2,ungpow,unscale,mark,unass,trace(i,1:nt)
enddo
111 continue
close(10)

sampling_interval = real(dt)/1.0e6

end subroutine read_sufile1

!------------------------------------------------------------------------------
subroutine read_sufile_2d(sufile, trace, n1, n2, dd1, dd2)

use global
use su

character(len=*), intent(in) :: sufile
integer, intent(in) :: n1, n2
real, intent(out) :: dd1, dd2
real, intent(out) :: trace(:,:)
integer :: nlen, i

nlen = nhead+n1
open(10,file=sufile,access='direct',recl=i4*nlen)
do i=1,n2
  read(10,rec=i,err=111) &
        tracl,tracr,fldr,tracf,ep,cdp,cdpt, &
        trid,nvs,nhs,duse,offset,gelev,selev,sdepth, &
        gdel,sdel,swdep,gwdep,scalel,scalco, &
        sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat, &
        tstat,laga,lagb,delrt,muts,mute,ns,dt, &
        gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
        stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
        lcs,hcs,year,day,hour,minute,sec,timbas, &
        trwf,grnors,grnofr,grnlof,gaps,otrav, &
        d1,f1,d2,f2,ungpow,unscale,mark,unass,trace(1:n1,i)
enddo
111 continue
close(10)

dd1 = real(dt)/1.0e6
dd2 = d2
if (abs(dd2) < 1.0e-6) dd2 = 1

end subroutine read_sufile_2d

!-----------------------------------------------------------------------------------
subroutine write_sufile_1d(sufile, trace, n1, dd1)

use global
use su

character(len=*), intent(in) :: sufile
integer, intent(in) :: n1
real, intent(in) :: dd1, trace(:)
integer :: nlen

call clean_suheader

nlen = nhead+n1
ns = n1
ep = 1
fldr = 1
dt = int(dd1*1.0e6)
d1 = dd1
f1 = 0
tracl = 0
tracr = 0

open(10,file=sufile,access='direct',recl=nlen*i4,status='replace')
tracl = tracl + 1
tracr = tracr + 1
cdp = 1
cdpt = 1
write(10,rec=1) &
        tracl,tracr,fldr,tracf,ep,cdp,cdpt, &
        trid,nvs,nhs,duse,offset,gelev,selev,sdepth, &
        gdel,sdel,swdep,gwdep,scalel,scalco, &
        sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat, &
        tstat,laga,lagb,delrt,muts,mute,ns,dt, &
        gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
        stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
        lcs,hcs,year,day,hour,minute,sec,timbas, &
        trwf,grnors,grnofr,grnlof,gaps,otrav, &
        d1,f1,d2,f2,ungpow,unscale,mark,unass, trace(1:n1)
close(10)

end subroutine write_sufile_1d

!-----------------------------------------------------------------------------------
subroutine write_sufile_2d(sufile, trace, n1, n2, dd1, dd2)

use global
use su

character(len=*), intent(in) :: sufile
integer, intent(in) :: n1, n2
real, intent(in) :: dd1, dd2, trace(:,:)
integer :: nlen, i

call clean_suheader

nlen = nhead+n1
ns = n1
ep = 1
fldr = 1
dt = int(dd1*1.0e6)
d1 = dd1
d2 = dd2
f1 = 0
f2 = 1.0
tracl = 0
tracr = 0

open(10,file=sufile,access='direct',recl=nlen*i4,status='replace')
do i = 1,n2
  tracl = tracl + 1
  tracr = tracr + 1
  cdp = i
  cdpt = i
  write(10,rec=i) &
        tracl,tracr,fldr,tracf,ep,cdp,cdpt, &
        trid,nvs,nhs,duse,offset,gelev,selev,sdepth, &
        gdel,sdel,swdep,gwdep,scalel,scalco, &
        sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat, &
        tstat,laga,lagb,delrt,muts,mute,ns,dt, &
        gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
        stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
        lcs,hcs,year,day,hour,minute,sec,timbas, &
        trwf,grnors,grnofr,grnlof,gaps,otrav, &
        d1,f1,d2,f2,ungpow,unscale,mark,unass, trace(1:n1,i)
enddo
close(10)

end subroutine write_sufile_2d

!-----------------------------------------------------------------------------------
subroutine write_sufile_2d_offset(sufile, starting_record, is, trace, n1, n2, dd1, dd2, coord)

use global
use datatype
use su

type(acquisition), intent(in) :: coord
character(len=*), intent(in) :: sufile
integer, intent(in) :: starting_record, is, n1, n2
real, intent(in) :: dd1, dd2, trace(:,:)
integer :: nlen, i

call clean_suheader

nlen = nhead+n1
ns = n1
ep = 1
fldr = is
dt = int(dd1*1.0e6)
d1 = dd1
d2 = dd2
f1 = 0
f2 = 1.0
tracl = 0
tracr = 0
scalel = -10
scalco = -10

open(10,file=sufile,access='direct',recl=nlen*i4)
do i = 1,n2
  tracl = tracl + 1
  tracr = tracr + 1
  cdp = i
  cdpt = i
  sx    = nint(coord%xs(is)*10.0)
  sy    = 0
  selev = nint(coord%zs(is)*10.0)
  gx    = nint(coord%xg(is,i)*10.0)
  gy    = 0
  gelev = nint(coord%zg(is,i)*10.0)
  
  write(10,rec=starting_record+i) &
        tracl,tracr,fldr,tracf,ep,cdp,cdpt, &
        trid,nvs,nhs,duse,offset,gelev,selev,sdepth, &
        gdel,sdel,swdep,gwdep,scalel,scalco, &
        sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat, &
        tstat,laga,lagb,delrt,muts,mute,ns,dt, &
        gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
        stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
        lcs,hcs,year,day,hour,minute,sec,timbas, &
        trwf,grnors,grnofr,grnlof,gaps,otrav, &
        d1,f1,d2,f2,ungpow,unscale,mark,unass, trace(1:n1,i)
enddo
close(10)

end subroutine write_sufile_2d_offset

!-----------------------------------------------------------------------------------
subroutine write_sufile_2d_coord(sufile, trace, n1, n2, dd1, dd2, xs, xg)

use global
use su

character(len=*), intent(in) :: sufile
integer, intent(in) :: n1, n2
real, intent(in) :: dd1, dd2, trace(:,:), xs, xg(n2)
integer :: nlen, i

call clean_suheader

nlen = nhead+n1
ns = n1
ep = 1
fldr = 1
dt = int(dd1*1.0e6)
d1 = dd1
d2 = dd2
f1 = 0
f2 = 1.0
tracl = 0
tracr = 0
sx = nint(xs)

open(10,file=sufile,access='direct',recl=nlen*i4,status='replace')
do i = 1,n2
  tracl = tracl + 1
  tracr = tracr + 1
!  cdp = i
!  cdpt = i
  gx = nint(xg(i))
  cdp = nint(0.5*(xs+xg(i)))
  cdpt = cdp
  offset = gx - sx
  write(10,rec=i) &
        tracl,tracr,fldr,tracf,ep,cdp,cdpt, &
        trid,nvs,nhs,duse,offset,gelev,selev,sdepth, &
        gdel,sdel,swdep,gwdep,scalel,scalco, &
        sx,sy,gx,gy,counit,wevel,swevel,sut,gut,sstat,gstat, &
        tstat,laga,lagb,delrt,muts,mute,ns,dt, &
        gain,igc,igi,corr,sfs,sfe,slen,styp,stas, &
        stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf, &
        lcs,hcs,year,day,hour,minute,sec,timbas, &
        trwf,grnors,grnofr,grnlof,gaps,otrav, &
        d1,f1,d2,f2,ungpow,unscale,mark,unass, trace(1:n1,i)
enddo
close(10)

end subroutine write_sufile_2d_coord
!-----------------------------------------------------------------------------------

subroutine write_binfile_1d(filename, data, n)

use global

character(len=*), intent(in) :: filename
integer, intent(in) :: n
real, intent(in) :: data(:)
integer :: i,it

open(10,file=filename,access='direct',recl=n*i4,status='replace')
write(10,rec=1) data(1:n)
close(10)

end subroutine write_binfile_1d
!-----------------------------------------------------------------------------------

subroutine write_binfile_2d(filename, data, n1, n2)

use global

character(len=*), intent(in) :: filename
integer, intent(in) :: n1, n2
real, intent(in) :: data(:,:)
integer :: i,it

open(10,file=filename,access='direct',recl=n1*i4,status='replace')

do i = 1, n2
   write(10,rec=i)(data(it,i), it = 1, n1)
enddo
close(10)

end subroutine write_binfile_2d

!-----------------------------------------------------------------------------------
subroutine write_binfile_3d(filename, data, n1, n2, n3)

character(len=*), intent(in) :: filename
integer,          intent(in) :: n1, n2, n3
real,             intent(in) :: data(:,:,:)
integer,          parameter  :: i4=4
integer                      :: i1, i2, i3

open(10,file=filename,access='direct',recl=n1*i4,status='replace')

do i3=1,n3
  do i2=1,n2
    write(10,rec=i2+(i3-1)*n2) data(:,i2,i3)
  enddo
enddo
close(10)

end subroutine write_binfile_3d

!-----------------------------------------------------------------------------------
subroutine read_binfile_1d(filename, data, n)

use global

implicit none

character(len=*), intent(in) :: filename
integer, intent(in)          :: n
real                         :: data(:)

open(10,file=filename,access='direct',recl=n*i4)
read(10, rec=1) data(1:n)
close(10)

end subroutine read_binfile_1d

!-----------------------------------------------------------------------------------
subroutine read_binfile_2d(filename, data, n1, n2)

use global

implicit none

character(len=*), intent(in) :: filename
integer, intent(in) :: n1, n2
real                :: data(:,:)
integer :: i1,i

open(10,file=filename,access='direct',recl=n1*i4)

do i = 1, n2
   read(10, rec=i)(data(i1,i),i1=1,n1)
enddo
close(10)

end subroutine read_binfile_2d

!-----------------------------------------------------------------------------------
subroutine read_binfile_3d(filename, data, n1, n2, n3)

character(len=*), intent(in)  :: filename
integer,          intent(in)  :: n1, n2, n3
real,             intent(out) :: data(:,:,:)
integer,          parameter   :: i4=4
integer                       :: i1, i2, i3

open(10,file=filename,access='direct',recl=n1*i4,status='old')

do i3=1,n3
  do i2=1,n2
    read(10,rec=i2+(i3-1)*n2) data(:,i2,i3)
  enddo
enddo
close(10)

end subroutine read_binfile_3d

!-----------------------------------------------------------------------------------
subroutine clean_suheader

use su
implicit none

integer i

tracl = 0
tracr = 0
fldr = 0
tracf = 0
ep = 0
cdp = 0
cdpt = 0
trid = 0
nvs = 0
nhs = 0
duse = 0
offset = 0
gelev = 0
selev = 0
sdepth = 0
gdel = 0
sdel = 0
swdep = 0
gwdep = 0
scalel = 0
scalco = 0
sx = 0
sy = 0
gx = 0
gy = 0
counit = 0
wevel = 0
swevel = 0
sut = 0
gut = 0
sstat = 0
gstat = 0
tstat = 0
laga = 0
lagb = 0
delrt = 0
muts = 0
mute = 0
ns = 0
dt = 0
gain = 0
igc = 0
igi = 0
corr = 0
sfs = 0
sfe = 0
slen = 0
styp = 0
stas = 0
stae = 0
tatyp = 0
afilf = 0
afils = 0
nofilf = 0
nofils = 0
lcf = 0
hcf = 0
lcs = 0
hcs = 0
year = 0
day = 0
hour = 0
minute = 0
sec = 0
timbas = 0
trwf = 0
grnors = 0
grnofr = 0
grnlof = 0
gaps = 0
otrav = 0
d1 = 0.0
f1 = 0.0
d2 = 0.0
f2 = 0.0
ungpow = 0.0
unscale = 0.0
mark = 0
unass(:) = 0

end subroutine clean_suheader

!-----------------------------------------------------------------------------------
subroutine free_surface(par, fs, npml)

use datatype
use global
use mmi_mpi

type(param), intent(in) :: par
integer, intent(in)     :: npml
integer, intent(out)    :: fs(:)
integer                 :: nx_pml, ix
real                    :: topo(par%nx)

nx_pml = par%nx + 2*npml
if (rank == 0) then

  topo = 0.0
  ! Set up free surface
  fs(npml+1:nx_pml-npml) = npml+1+int(topo/par%dx)
  fs(1:npml) = fs(npml+1)
  fs(nx_pml-npml+1:nx_pml) = fs(nx_pml-npml)
endif
call MPI_BCAST(fs,nx_pml,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

end subroutine free_surface

end module io

