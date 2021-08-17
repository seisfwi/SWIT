% WriteSuStructure : writes data to disk using SU-CWP format
%
% EX
% WriteSuStructure('datacube.segy',SegyHeader,SegyTraceHeaders,Data);

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
%
function WriteSuStructure(filename,SegyHeader,SegyTraceHeaders,Data)


  segyid = fopen(filename,'w');   % USE LOCAL ENDIAN FLAG,
                                    % I.E. ENDIAN FORMAT NOT SPECIFIED 
% segyid = fopen(filename,'w','b'); % BIG ENDIAN
                                    


% FORCE THE USE OF IEEE				    
SegyHeader.SegyFormatRevisionNumber=100;
SegyHeader.DataSampleFormat=5; % IEEE 
    
                                    
ntraces=size(Data,2);
% hw=waitbar(0,['Writing to SU-file : ',filename]);
for i=1:ntraces;
    if (i/200)==round(i/200),
      SegymatVerbose(['writing trace ',num2str(i),' of ',num2str(ntraces),', filepos=',num2str(ftell(segyid))],1)
      SegymatVerbose(SegyTraceHeaders(i).ns,1)
    end
    PutSegyTrace(segyid,Data(:,i),SegyTraceHeaders(i),SegyHeader);
%   waitbar(i/ntraces,hw);
end
% close(hw)
fclose(segyid);                                  
