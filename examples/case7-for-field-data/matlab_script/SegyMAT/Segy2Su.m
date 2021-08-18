% Segy2Su : Converts SEGY file to SU format
% 
% Call : Segy2Su(filename,ReadSegyOption)
%    Replaces the filename suffix to '.su';
%    'ReadSegyOptions' are the same as to 'ReadSegy'
%
%  See also : ReadSegy
%

%
%
% (C) 2001,2002 Thomas Mejer Hansen, thomas@cultpenguin.com
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
function Segy2Su(filename,varargin);
[Data,SegyTraceHeaders,SegyHeader]=ReadSegy(filename,varargin);
SegyHeader.DataSampleFormat=5;       % IEEE
SegyHeader.SegyFormatRevisionNumber=100; % SEGY REVISION 1
WriteSuStructure([filename,'.su'],SegyHeader,SegyTraceHeaders,Data);
