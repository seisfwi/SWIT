% GetSegyTrace : Reads a seg y trace, data and header
%
% [SegyTraceHeader,SegyData]=GetSegyTrace(segyid,TraceStart,DataFormat,ns);
%

%
% (C) 2001-2004 Thomas Mejer Hansen, tmh@gfy.ku.dk/thomas@cultpenguin.com
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

function [SegyTraceHeader,SegyData]=GetSegyTrace(segyid,TraceStart,DataFormat,ns,SegyTraceHeader);

if exist('DataFormat')==0, DataFormat='float32'; end
if exist('TraceStart')==0, TraceStart=ftell(segyid); end

if exist('SegyTraceHeader')
    if isempty('SegyTraceHeader');
        clear SegyTraceHeader;
    end
end

[SegyTraceHeader]=GetSegyTraceHeader(segyid,TraceStart,DataFormat,ns);
[SegyTraceData]=GetSegyDataTraceHeader(segyid,TraceStart,DataFormat,ns);
