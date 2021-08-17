% GetSegyTraceData : Get Segy trace data if filehandle
% 
% Call : 
%
%   tracedata=GetSegyTraceData(segyid,ns,SegyHeader,SkipData
%

%
% (C) 2001-2004 Thomas Mejer Hansen, tmh@gfy.ku.dk/thomas@cultpenguin.com
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
function tracedata=GetSegyTraceData(segyid,ns,SegyHeader,SkipData)
%
% Get Segy Trace Data. 
%
%

SkipData=[];  
  
if nargin==0,
    SegymatVerbose(exist('segyid'))
    SegymatVerbose([mfilename,' : SEGYID not specified - exiting'])
    tracedata=[];
    return
end
if nargin==1
  SegymatVerbose([mfilename,' : NS not specified - exiting'])
  tracedata=[];
  return
end
if nargin==2
    SegyHeader.DataFormat='float32'; 
    SegyHeader.BytesPerSample=32;
    SegyHeader.DataSampleFormat=5; % IEEE
    SegymatVerbose(['Dataformat not specified -  dataformat->',SegyHeader.DataFormat])
end;

if isempty(SkipData)==1, SkipData=0; end

Revision=SegyHeader.SegyFormatRevisionNumber;
if Revision>0, Revision=1; end
Format=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).format;  

BPS=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).bps;  
if (SkipData==1)
    SkipBytes=ns*BPS/8;
    fseek(segyid,SkipBytes,'cof');
    tracedata=[];
else  
    % SegymatVerbose([mfilename,' : ',Format,'. ns=',num2str(ns)])
    try 
      tracedata=fread(segyid,ns,Format);
    catch
      SegymatVerbose([mfilename,' : Error using fread - Possibly ''ns'' is negative -' ...
	    ' check byteorder-'])
      tracedata=[];
    end
    
    
    if (strcmp(Format,'uint32')==1), % IBM FLOATING POINT
        % CONVERT FROM FLOATING POINT
        verbose=1;
        if verbose>1, SegymatVerbose([mfilename,'Converting from IBM, DataFormat :',SegyHeader.DataFormat]); end
        try
          tracedata=ibm2num(uint32(tracedata));
        catch
          % SegymatVerbose([mfilename,' : SOMETHING BAD HAPPENED WHEN CONVERTING FROM IBM FLOATS TO IEEE. ARE YOU SURE DATA ARE IBM FLOAT FORMATTED ?' ])
          % tracedata=0.*tracedata;
          % return

        end
    end;
    
end
