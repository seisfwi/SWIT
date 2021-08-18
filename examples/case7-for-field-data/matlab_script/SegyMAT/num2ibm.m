function b=num2ibm(x)
% num2ibm : convert IEEE 754 doubles to IBM 32 bit floating point format
%    b=num2ibm(x)
% x is a matrix of doubles
% b is a corresponding matrix of uint32
%
% The representations for NaN and inf are arbitrary
%
% See also ibm2num

% 
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%
% (C) Brian Farrelly, 22 October 2001
%  mailto:Brian.Farrelly@nho.hydro.com          Norsk Hydro Research Centre
%  phone +47 55 99 68 74                 (((                  Postboks 7190
%  fax   +47 55 99 69 70                2oooS                 N-5020 Bergen
%  home  +47 55 13 78 49                HYDRO                        Norway
%

b=repmat(uint32(0),size(x));
err=zeros(size(x));

%format long

x(x> 7.236998675585915e+75)= inf;        % change big numbers to infinity
x(x<-7.236998675585915e+75)=-inf;        % 7.236998675585915e+75 is
                                         %    ibm2num(uint32(hex2dec('7fffffff')) or
                                         %    ibm2num(num2ibm(inf))

[F E]=log2(abs(x));

e=E/4;                         % exponent of base 16
ec=ceil(e);                    % adjust upwards to integer
p=ec+64;                       % offset exponent

f=F.*2.^(-4*(ec-e));           % correct mantissa for fractional part of exponent
f=round(f*2^24);               % convert to integer. Roundoff here can be as large as
                               % 0.5/2^20 when mantissa is close to 1/16 so that
                               % 3 bits of signifance are lost.

p(f==2^24)=p(f==2^24)+1;       % Roundoff can cause f to be 2^24 for numbers just under a
f(f==2^24)=2^20;               % power of 16, so correct for this

%format hex
psi=uint32(p*2^24);            % put exponent in first byte of psi.
phi=uint32(f);                 % put mantissa into last 3 bytes of phi 

% make bit representation

b=bitor(psi,phi);                        % exponent and mantissa
b(x<0)=bitset(b(x<0),32);                % sign bit 
%format long

% special cases

b(x==0)          =uint32(0)                  ;         %  bias is incorrect for zero 
b(isnan(x))      =uint32(hex2dec('7fffffff'));         %  7.237005145973116e+75 in IBM format
b(isinf(x) & x>0)=uint32(hex2dec('7ffffff0'));         %  7.236998675585915e+75    ,,
b(isinf(x) & x<0)=uint32(hex2dec('fffffff0'));         % -7.236998675585915e+75    ,,
                                                       % Note that NaN > inf in IBM format

% check bit representation for normal cases 

checkx=ibm2num(b);                      % note that use of base 16 in IBM format
z=x==0;                                 % can lead to a loss of 3 bits of precision
err(z)=0;                               % compared with an IEEE single.
q=(checkx(~z)-x(~z))./x(~z);
err(~z) = abs(q) > 5e-7;                % this is almost reached with numbers
                                        % of the form 16^n + 0.5*16^(n-5) where
                                        % the mantissa is 100001 hex. Roundoff
                                        % error is then 0.5/16^5=0.5/2^20=4.7684e-7
                                                    
                                          
if any(err)
   warning('Conversion error in num2ibm for the following:')
   disp(x(logical(err)))
end   




