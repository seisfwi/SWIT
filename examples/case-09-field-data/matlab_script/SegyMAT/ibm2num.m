function x=ibm2num(b)
% ibm2num : convert IBM 32 bit floating point format to doubles
%    x=num2ibm(b)
% b is a matrix of uint32
% x is a corresponding matrix of doubles
%
%
% See also num2ibm

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



x=repmat(NaN,size(b));

sign=bitget(b,32);                            % get sign from first bit
sign=double(sign);

% format hex
exp=bitand(b,uint32(hex2dec('7f000000')));    % get exponent from first byte, last 7 bits
exp=bitshift(exp,-24);
%format long

exp=double(exp)- 64;                          % remove bias from exponent 

%format hex
frac=bitand(b,uint32(hex2dec('00ffffff')));   % get mantissa from last 3 bytes
%format long
frac=double(frac);
frac=frac/2^24;


x=(1-2*sign).*16.^exp .* frac;

err = frac==0 & (exp~=-64 | sign~=0);         % bias removal is incorrect for zero
if any(err)
   % TMH 19/06/2003
   disp(['WARNING : ',mfilename,' Invalid zero input --> Sure data are IBM FLOAT formatted ?'])	
   return;							     
   %warning('Invalid zero input in ibm2num for the following:')
   % format hex; disp(b(err)); format							      
end

err = frac~=0 & (frac<1/16 | frac>=1);
if any(err)
   % TMH 19/06/2003
   disp(['WARNING : ',mfilename,' Invalid mantissa input --> Sure data are IBM FLOAT formatted ?'])	
   return;
   % warning('Invalid mantissa input in ibm2num for the following:')
   % format hex; disp(b(err)); format							      
end   




