% SegymatVersion - Returns the version and release date
%
% [ver,d]=SegymatVersion;
%

%
% (C) 2001-2011, Thomas Mejer Hansen, thomas.mejer.hansen@gmail.com
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

function [ver,d]=SegymatVersion;
  version='1.5.1';
  releasedate='October 28, 2011';
  
  SegymatVerbose(['This is SegyMAT version ',version,' - released ',releasedate],-1)    
  ver=version;
  d=releasedate;
  