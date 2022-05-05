% GetSegyTraceHeaderInfo : Returns a array of a SEGY Y TraceHeader value
%
% Call :
% [value]=GetSegyHeaderInfo(SegyTraceHeaders,header)
%
% header is a header value like 'cdp','dt','TraceNumber'
%

%
% (C) 2001-2004, Thomas Mejer Hansen, tmh@gfy.ku.dk/thomas@cultpenguin.com
%
%
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
function [value]=GetSegyTraceHeaderInfo(SegyTraceHeaders,header)
  
  try
    value=[SegyTraceHeaders.(header)];
  catch
    SegymatVerbose(['Could not get  ',header,' header information from supplied SegyTraceHeaders']);
  end
