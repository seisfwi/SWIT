% SegymatVerbose : Writes out verbose information to the screen
%
%
% Call : 
%   SegymatVerbose(text,verboselevel)
%   prints out 'text' to screen if verboselevel is higher than threshold
%   set in m-file.
%

%
% (C) 2001-2004, Thomas Mejer Hansen, tmh@gfy.ku.dk/thomas@cultpenguin.com
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


function SegymatVerbose(txt,level)

  if nargin==0, return, end
  if nargin==1, 
    level=1; 
  end
  
  VerboseLevel=0;
  
  % Only print information if at or above VerboseLevel
  if level<=VerboseLevel
    disp(sprintf('%s  : %s','SegyMAT',txt))
  end
  