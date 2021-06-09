! module for SU file structure
!
module su

integer, parameter :: nhead=60

! 1-4	long tracl;	/* trace sequence number within line */
! 5-8	long tracr;	/* trace sequence number within reel */
! 9-12	long fldr;	/* field record number */
!13-16	long tracf;	/* trace number within field record */
!17-20	long ep;	/* energy source point number */
!21-24	long cdp;	/* CDP ensemble number */
!25-28	long cdpt;	/* trace number within CDP ensemble */

integer*4 :: tracl,tracr,fldr,tracf,ep,cdp,cdpt

!29-30	short trid;	/* trace identification code:
!			1 = seismic data
!			2 = dead
!			3 = dummy
!			4 = time break
!			5 = uphole
!			6 = sweep
!			7 = timing
!			8 = water break
!			9---, N = optional use (N = 32,767)
!
!			Following are CWP id flags:
!
!			 9 = autocorrelation
!
!			10 = Fourier transformed - no packing
!			     xr[0],xi[0], ..., xr[N-1],xi[N-1]
!
!			11 = Fourier transformed - unpacked Nyquist
!			     xr[0],xi[0],...,xr[N/2],xi[N/2]
!
!			12 = Fourier transformed - packed Nyquist
!	 		     even N:
!			     xr[0],xr[N/2],xr[1],xi[1], ...,
!				xr[N/2 -1],xi[N/2 -1]
!				(note the exceptional second entry)
!			     odd N:
!			     xr[0],xr[(N-1)/2],xr[1],xi[1], ...,
!				xr[(N-1)/2 -1],xi[(N-1)/2 -1],xi[(N-1)/2]
!				(note the exceptional second & last entries)
!
!			13 = Complex signal in the time domain
!			     xr[0],xi[0], ..., xr[N-1],xi[N-1]
!
!			14 = Fourier transformed - amplitude/phase
!			     a[0],p[0], ..., a[N-1],p[N-1]
!
!			15 = Complex time signal - amplitude/phase
!			     a[0],p[0], ..., a[N-1],p[N-1]
!
!			16 = Real part of complex trace from 0 to Nyquist
!
!			17 = Imag part of complex trace from 0 to Nyquist
!
!			18 = Amplitude of complex trace from 0 to Nyquist
!
!			19 = Phase of complex trace from 0 to Nyquist
!
!			21 = Wavenumber time domain (k-t)
!
!			22 = Wavenumber frequency (k-omega)
!
!			30 = Depth-Range (z-x) traces
!
!			101 = Seismic data packed to bytes (by supack1)
!			
!			102 = Seismic data packed to 2 bytes (by supack2)
!
!
!31-32	short nvs;	/* number of vertically summed traces (see vscode
!			   in bhed structure) */
!
!33-34	short nhs;	/* number of horizontally summed traces (see vscode
!			   in bhed structure) */
!
!35-36	short duse;	/* data use:
!				1 = production
!				2 = test */

integer*2 :: trid,nvs,nhs,duse 

!37-40	long offset;	/* distance from source point to receiver
!			   group (negative if opposite to direction
!			   in which the line was shot) */
!
!41-44	long gelev;	/* receiver group elevation from sea level
!			   (above sea level is positive) */
!
!45-48	long selev;	/* source elevation from sea level
!			   (above sea level is positive) */
!
!49-52	long sdepth;	/* source depth (positive) */
!
!53-56	long gdel;	/* datum elevation at receiver group */
!
!57-60	long sdel;	/* datum elevation at source */
!
!61-64	long swdep;	/* water depth at source */
!
!65-68	long gwdep;	/* water depth at receiver group */

integer*4 :: offset,gelev,selev,sdepth,gdel,sdel,swdep,gwdep

!69-70	short scalel;	/* scale factor for previous 7 entries
!			   with value plus or minus 10 to the
!			   power 0, 1, 2, 3, or 4 (if positive,
!			   multiply, if negative divide) */
!
!71-72	short scalco;	/* scale factor for next 4 entries
!			   with value plus or minus 10 to the
!			   power 0, 1, 2, 3, or 4 (if positive,
!			   multiply, if negative divide) */

integer*2 :: scalel,scalco

!73-76	long  sx;	/* X source coordinate */
!
!77-80	long  sy;	/* Y source coordinate */
!
!81-84	long  gx;	/* X group coordinate */
!
!85-88	long  gy;	/* Y source coordinate */

integer*4 :: sx,sy,gx,gy

!89-90	short counit;	/* coordinate units code:
!				for previoius four entries
!				1 = length (meters or feet)
!				2 = seconds of arc (in this case, the
!				X values are longitude and the Y values
!				are latitude, a positive value designates
!				the number of seconds east of Greenwich
!				or north of the equator */
!
!91-92	short wevel;	/* weathering velocity */
!
!93-94	short swevel;	/* subweathering velocity */
!
!95-96	short sut;	/* uphole time at source */
!
!97-98	short gut;	/* uphole time at receiver group */
!
!99-100	short sstat;	/* source static correction */
!
!101-102  short gstat;	/* group static correction */
!
!103-104  short tstat;	/* total static applied */
!
!105-106  short laga;	/* lag time A, time in ms between end of 240-
!			   byte trace identification header and time
!			   break, positive if time break occurs after
!			   end of header, time break is defined as
!			   the initiation pulse which maybe recorded
!			   on an auxiliary trace or as otherwise
!			   specified by the recording system */
!
!107-108  short lagb;	/* lag time B, time in ms between the time break
!		   and the initiation time of the energy source,
!		   may be positive or negative */
!
!109-110  short delrt;	/* delay recording time, time in ms between
!			   initiation time of energy source and time
!			   when recording of data samples begins
!			   (for deep water work if recording does not
!			   start at zero time) */
!
!111-112  short muts;	/* mute time--start */
!
!113-114  short mute;	/* mute time--end */
!
!115-116  unsigned short ns;	/* number of samples in this trace */
!
!117-118  unsigned short dt;	/* sample interval; in micro-seconds */

integer*2 :: counit,wevel,swevel,sut,gut,sstat,gstat,   &
             tstat,laga,lagb,delrt,muts,mute,ns,dt 

!119-120  short gain;	/* gain type of field instruments code:
!				1 = fixed
!				2 = binary
!				3 = floating point
!				4 ---- N = optional use */
!
!121-122  short igc;	/* instrument gain constant */
!
!123-124  short igi;	/* instrument early or initial gain */
!
!125-126  short corr;	/* correlated:
!				1 = no
!				2 = yes */
!
!127-128  short sfs;	/* sweep frequency at start */
!
!129-130  short sfe;	/* sweep frequency at end */
!
!131-132  short slen;	/* sweep length in ms */
!
!133-134  short styp;	/* sweep type code:
!				1 = linear
!				2 = cos-squared
!				3 = other */
!
!135-136  short stas;	/* sweep trace length at start in ms */
!
!137-138  short stae;	/* sweep trace length at end in ms */
!
!139-140  short tatyp;	/* taper type: 1=linear, 2=cos^2, 3=other */
!
!141-142  short afilf;	/* alias filter frequency if used */
!
!143-144  short afils;	/* alias filter slope */
!
!145-146  short nofilf;	/* notch filter frequency if used */
!
!147-148  short nofils;	/* notch filter slope */
!
!149-150  short lcf;	/* low cut frequency if used */
!
!151-152  short hcf;	/* high cut frequncy if used */
!
!153-154  short lcs;	/* low cut slope */
!
!155-156  short hcs;	/* high cut slope */
!
!157-158  short year;	/* year data recorded */
!
!159-160  short day;	/* day of year */
!
!161-162  short hour;	/* hour of day (24 hour clock) */
!
!163-164  short minute;	/* minute of hour */
!
!165-166  short sec;	/* second of minute */
!
!167-168  short timbas;	/* time basis code:
!				1 = local
!				2 = GMT
!				3 = other */
!
!169-170  short trwf;	/* trace weighting factor, defined as 1/2^N
!			   volts for the least sigificant bit */
!
!171-172  short grnors;	/* geophone group number of roll switch
!			   position one */
!
!173-174  short grnofr;	/* geophone group number of trace one within
!			   original field record */
!
!175-176  short grnlof;	/* geophone group number of last trace within
!			   original field record */
!
!177-178   short gaps;	/* gap size (total number of groups dropped) */
!
!179-180   short otrav;	/* overtravel taper code:
!				1 = down (or behind)
!				2 = up (or ahead) */

integer*2 :: gain,igc,igi,corr,sfs,sfe,slen,styp,stas,       &
             stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf,   &
             lcs,hcs,year,day,hour,minute,sec,timbas,        &
             trwf,grnors,grnofr,grnlof,gaps,otrav

!	/* local assignments */
!181-184  float d1;	/* sample spacing for non-seismic data */
!
!185-188  float f1;	/* first sample location for non-seismic data */
!
!189-192  float d2;	/* sample spacing between traces */
!
!193-196  float f2;	/* first trace location */
!
!197-200  float ungpow;	/* negative of power used for dynamic
!	   range compression */
!
!201-204  float unscale;	/* reciprocal of scaling factor to normalize
!	   range */

real*4 :: d1,f1,d2,f2,ungpow,unscale

!205-206	short mark;	/* mark selected traces */
!
!207-240  short unass[17];	/* unassigned--NOTE: last entry causes 
!			   a break in the word alignment, if we REALLY
!			   want to maintain 240 bytes, the following
!			   entry should be an odd number of short/UINT2
!			   OR do the insertion above the "mark" keyword
!			   entry */
integer*2 :: mark,unass(17)

end module su

