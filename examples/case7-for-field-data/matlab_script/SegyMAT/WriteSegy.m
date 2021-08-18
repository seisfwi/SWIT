% WriteSegy : writes data to disk using SEGY REV 1 standard.
%
% EX
% WriteSegy('datacube.segy',data,'dt',.004,'Inline3D',Inline,'Crossline3D',Crossline,'cdpX',X,'cdpY',Y);
%
% to use a specific SEG revision use :
% WriteSegy('test.segy',seisdata,'revision',0); % SEG-Y Revision 0
% WriteSegy('test.segy',seisdata,'revision',1); % SEG-Y Revision 1
% 
% to use a specific Data Sampling Format use :
% WriteSegy('test.segy',seisdata,'dsf',1); % IBM FLAOTING POINT
%
% Forice Revision 1 and IEEE Floating point :
% WriteSegy('test.segy',seisdata,'dsf',5,'revision',1); 
%
% See also : WriteSegyStructure, WriteSu, WriteSuStructure
%

%
% (C) 2001-2007, Thomas Mejer Hansen, tmh@gfy.ku.dk/thomas.mejer.hansen@gmail.com
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
function WriteSegy(filename,data,varargin);
    
    [ns,ntraces]=size(data);
    
    
    for i=1:2:length(varargin)
        var=varargin{i};
        val=varargin{i+1};
        eval([var,'=[',num2str(val),'];']);
    end  
    
    if exist('dt')~=1
        dt=0.004;
        SegymatVerbose([mfilename,' - No dt set. Using dt=',num2str(dt)],0)
    end
    
    %%%%%%%%%%
    % SET UP SegyHeader structure.
    if exist('TextualFileHeader'), SegyHeader.TextualFileHeader=TextualFileHeader; end
    if exist('Job')==1, SegyHeader.Job=Job; end;
    if exist('Line')==1, SegyHeader.Line=Line; end
    if exist('Reel')==1, SegyHeader.Reel=Reel; end
    if exist('DataTracePerEnsemble')==1, SegyHeader.DataTracePerEnsemble=DataTracePerEnsemble; end
    if exist('AuxiliaryTracePerEnsemble')==1, SegyHeader.AuxiliaryTracePerEnsemble=AuxiliaryTracePerEnsemble; end
    if exist('dt')==1, SegyHeader.dt=dt.*1e+6; end
    if exist('dtOrig')==1, SegyHeader.dtOrig=dtOrig; end
    if exist('ns')==1, SegyHeader.ns=ns; end
    if exist('nsOrig')==1, SegyHeader.nsOrig=nsOrig; end
    if exist('EnsembleFold')==1, SegyHeader.EnsembleFold=EnsembleFold;  end
    if exist('TraceSorting')==1, SegyHeader.TraceSorting=TraceSorting; end
    if exist('VerticalSumCode')==1, SegyHeader.VerticalSumCode=VerticalSumCode; end
    
    %%%%%%%%%%
    %% CHOOSE WHICH SEGY REVISION AND DATATYPE TO USE
    %% DEFAULT IS REVISION 1, 4-byte IEEE
    
    % IF A SPECFIC REVISION HAS BEEN CHOSEN, USE THAT
    if exist('revision')==1,
        if revision==1, 
            SegyHeader.SegyFormatRevisionNumber=100;
        else
            SegyHeader.SegyFormatRevisionNumber=0;
        end
        SegymatVerbose([mfilename,' :  Using user specified SEG Y revision : ',num2str(revision)],1)
    end
    
    
    % IF A SPECFIC DATA SAMPLING FORMAT HAS BEEN SELECTED USE THAT
    if exist('dsf')==1,
        SegyHeader.DataSampleFormat=dsf;
        SegymatVerbose([mfilename,' :  Using user specified Data Sample Format : ',num2str(revision)],1)
    end
    
    
    
    % UNCOMMENT THE FOLLOWING TWO LINES TO USE REVISION 1 (2002)
    %SegyHeader.SegyFormatRevisionNumber=100; % 2002 SEG Y STYLE
    %SegyHeader.DataSampleFormat=2; % '1'->4-byte IBM floating point 
    %                                 '2'->4-byte two's complement integer 
    %                                 '3'->2-byte two's complement integer
    %                                 '5'->4-byte IEEE floating point (default)
    %                                 '8'->1-byte two's complement integer
    
    % UNCOMMENT THE FOLLOWING TWO LINES TO USE REVISION 0 (1975)
    %SegyHeader.SegyFormatRevisionNumber=0; % 1975 SEoG Y STYLE
    %SegyHeader.DataSampleFormat=1; % '1'->4-byte IBM Floating Point
    
    % OPEN SEGY FILE HANDLE
    segyid = fopen(filename,'w','b');   % ALL DISK FILES ARE IN BIG
                                        % ENDIAN FORMAT, ACCORDING SEG
                                        % Y rev 1

    % Write SEGY HEADER
    SegyHeader=PutSegyHeader(segyid,SegyHeader);
    
    for i=1:ntraces;
        if (i/100)==round(i/100),
            SegymatVerbose(['writing trace ',num2str(i),' of ',num2str(ntraces)],0);
        end
        % Basic TraceHeader information
        
        % INITALIZE SEGY TRACE HEADER
        if exist('SegyTraceHeader')==0;
            SegyTraceHeader=InitSegyTraceHeader(ns,dt*1e+6);
        end
        
        if exist('TraceNumber')==0
            SegyTraceHeader.TraceNumber=i; 
            SegyTraceHeader.TraceSequenceFile=i; 
        end 
        
        % Update TraceHeader information if available
        if exist('cdpX')==1,SegyTraceHeader.cdpX = cdpX(i);end
        if exist('offset')==1,SegyTraceHeader.offset = offset(i);end
        if exist('cdpY')==1,SegyTraceHeader.cdpY = cdpY(i);end
        if exist('Inline3D')==1,SegyTraceHeader.Inline3D = Inline3D(i);end
        if exist('Crossline3D')==1,SegyTraceHeader.Crossline3D=Crossline3D(i);end

        if exist('YearDataRecorded')==1,SegyTraceHeader.YearDataRecorded=YearDataRecorded(i);end
        if exist('DayOfYear')==1,SegyTraceHeader.DayOfYear=DayOfYear(i);end
        if exist('HourOfDay')==1,SegyTraceHeader.HourOfDay=HourOfDay(i);end
        if exist('MinuteOfOur')==1,SegyTraceHeader.MinuteOfOur=MinuteOfOur(i);end
        if exist('SecondOfMinute')==1,SegyTraceHeader.SecondOfMinute=SecondOfMinute(i);end
        if exist('TimeBaseCode')==1,SegyTraceHeader.TimeBaseCode=TimeBaseCode(i);end

        
        
        
        % Write the Trace
        PutSegyTrace(segyid,data(:,i),SegyTraceHeader,SegyHeader);
    end
    
    % CLOSE SEGY FILE HANDLE                                 
    
    fclose(segyid);                                  
    