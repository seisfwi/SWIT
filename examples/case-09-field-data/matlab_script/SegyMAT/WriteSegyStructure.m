% WriteSegyStructure : writes data to disk using SEGY REV 0 and 1 standards.
%
% EX
% WriteSegyStructure('datacube.segy',SegyHeader,SegyTraceHeaders,Data);
%
% To force the use of SEG Y revision 0
% WriteSegyStructure('datacube.segy',SegyHeader,SegyTraceHeaders,Data,'revision',0);
% To force the use of SEG Y revision 1
% WriteSegyStructure('datacube.segy',SegyHeader,SegyTraceHeaders,Data,'revision',1);
% To force the data sampling format to be IBM Floating Point 
% WriteSegyStructure('datacube.segy',SegyHeader,SegyTraceHeaders,Data,'dsf',1);
% 
% To force the use of SEG Y revision 0 and data sampling format IEEE :
% WriteSegyStructure('datacube.segy',SegyHeader,SegyTraceHeaders,Data,'revision',1,'dsf',5);
%
% See the dokumentation for for proper values of 'dsf'
%
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
%
%
function SegyHeader=WriteSegyStructure(filename,SegyHeader,SegyTraceHeaders,Data,varargin)

for i=1:2:length(varargin)
  var=varargin{i};
  val=varargin{i+1};
  eval([var,'=[',num2str(val),'];']);
end  
  
  
if exist('revision')==1,
  if revision==0,
    SegyHeader.SegyFormatRevisionNumber=0;
  else
    SegyHeader.SegyFormatRevisionNumber=100;
  end  
  SegymatVerbose([mfilename,' : Using SEG Y revision ',num2str(revision)])
end

if exist('dsf'),
  SegyHeader.DataSampleFormat=dsf;
  SegymatVerbose([mfilename,' : Using Data Sample Format ',num2str(dsf)])
end

  
  

segyid = fopen(filename,'w','b');   % ALL DISK FILES ARE IN BIG
                                    % ENDIAN FORMAT, ACCORDING SEG
                                    % Y rev 1

% JUST SOME INFORMATION TO WRITE TO SCREEN :
Revision=SegyHeader.SegyFormatRevisionNumber;
if Revision>0, Revision=1; end
%Format=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).name;
%txt=['SegyRevision ',sprintf('%0.4g',Revision),', ',Format,'(',num2str(SegyHeader.DataSampleFormat),')'];
txt='';
				    
SegyHeader=PutSegyHeader(segyid,SegyHeader);

ntraces=size(Data,2);
hw=waitbar(0,['Writing to SEGY-file : ',filename,' - ',txt]);


for i=1:ntraces;
    if (i/100)==round(i/100),
      SegymatVerbose(['writing trace ',num2str(i),' of ',num2str(ntraces),', filepos=',num2str(ftell(segyid))])
      waitbar(i/ntraces,hw)
    end
    PutSegyTrace(segyid,Data(:,i),SegyTraceHeaders(i),SegyHeader);

end
close(hw)
fclose(segyid);                                  
