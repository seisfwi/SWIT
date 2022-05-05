% SegyMAT_GAIN : Gain plugin for SegyMAT
%
% [Data,SegyTraceHeaders,SegyHeader]=SegyMAT_GAIN(Data,SegyTraceHeaders,SegyHeader,varargin);
%
% ex. AGC using AGC window of 100 ms :
% [Data]=SegyMAT_GAIN(Data,SegyTraceHeaders,SegyHeader,'agc',.1);
% ex. apply t^(pow), pow=2
% [Data]=SegyMAT_GAIN(Data,SegyTraceHeaders,SegyHeader,'pow',2);
% 
%
% (C) Thomas Mejer Hansen (thomas@cultpenguin.com), 2002
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
function [Data,SegyTraceHeaders,SegyHeader]=SegyMAT_GAIN(Data,SegyTraceHeaders,SegyHeader,varargin);

nv=0;
for i=1:length(varargin)
  if strcmp(varargin{i},'agc') 
    type=1;
    window=varargin{i+1};
  end
  if strcmp(varargin{i},'pow') 
    type=2;
    pow=varargin{i+1};
  end
end


if type==1
  disp([mfilename,' AGC'])
  for it=1:size(Data,2)
    disp([mfilename,' AGC trace : ',num2str(it)])
    
    % Trace Data
    TraceData=Data(:,it);
    
    % Window Length
    nsw=round(window./(SegyTraceHeaders(it).dt./1e+6));
    nshalf=floor(nsw/2);
    startsample=ceil(nsw/2);
    endsample=length(TraceData)-floor(nsw);
    for is=startsample:1:endsample;
      
      range=[is-nshalf+1:1:is+nshalf];
      
      gain=mean(abs(TraceData(range)));
      if gain~=0
	
        TraceData(is)=TraceData(is)./mean(abs(TraceData(range)));	

	% APPLY TO TOP
	if is==startsample
	  for i=[1:startsample-1];
	    TraceData(i)=TraceData(i)./mean(abs(TraceData(range)));
	  end
	end
	
	if is==endsample
	  for i=[endsample+1:1:length(TraceData)];
	    TraceData(i)=TraceData(i)./gain;
	  end
	end
	
	
      end
      
    end
    
    Data(:,it)=TraceData;
    
  end % END LOOP OVER TRACES
end % END TYPE


if type==2
  disp([mfilename,' POW'])
  for it=1:size(Data,2)
    % disp([mfilename,' POW trace : ',num2str(it)])
    t=[1:1:SegyTraceHeaders(it).ns]*SegyTraceHeaders(it).dt./1e+6+SegyTraceHeaders(it).DelayRecordingTime./1e+3;
    tp=t.^(pow)';
    Data(:,it)=Data(:,it).*tp;
  end
end