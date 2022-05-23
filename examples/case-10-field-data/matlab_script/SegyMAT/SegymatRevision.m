% SegymatRevision - Returns the revision history
%
% Call : [Revision]=SegymatRevision
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
%

function [Rev]=SegymatRevision;
  
  revisionfile=[fileparts(which('ReadSegy')),'/REVISION'];
  
  fid=fopen(revisionfile,'r');
  cl=0;
  while 1
    cl=cl+1;
    currentline=fgetl(fid);
    if ~ischar(currentline), break, end
    st{cl} = currentline;
  end
  fclose(fid);
  
  if nargout==0,
    SegymatVerbose(['SegyMAT Revision History'])
    for i=1:length(st);
      SegymatVerbose(st{i})
    end
  else
    Rev=st;
  end
  