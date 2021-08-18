% ReadSu : Reads a SU formatted file (Seismic Unix)
%
% Call :
% [Data,SuTraceHeaders,SuHeader]=ReadSu(filename);
%
% To read in big endian format (default):
% [Data,SuTraceHeaders,SuHeader]=ReadSu(filename,'endian','b');
% To read in little endian format :
% [Data,SuTraceHeaders,SuHeader]=ReadSu(filename,'endian','l');
%
%
% To read in trace data as 'int32' :
% [Data,SuTraceHeaders,SuHeader]=ReadSu(filename,'DataFormat','int32');
% To read time slice 0.5<t<5 :
% [Data,SuTraceHeaders,SuHeader]=ReadSu(filename,'trange',.5,3);
% Skip every 5th trace :
% [Data,SuTraceHeaders,SuHeader]=ReadSu(filename,'jump',5);
% Read data in a CDP header range : 5000<cdp<5800 
% (change cdp to any other valid TraceHeader value)
% [Data,SuTraceHeaders,SuHeader]=ReadSu(filename,'minmax','cdp'5000,5800);
%
% Combine any combination of the above
% [Data,SuTraceHeaders,SuHeader]=ReadSu(filename,'jump',1,'minmax','cdp',5300,5400);
%
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
%
% 

%
% 0.1 : INitial Release
% 0.2 : Added SkipData var, to skip reading of data.
% 0.3 : May 01, 2002 
%       Added ability to read in ever 'jump' traces.
%       Added ability to read in time range.
%       Added abiliy to read in header range (ex. mincdp to maxcdp).
%

function [Data,SuTraceHeaders,SuHeader]=ReadSu(filename,varargin);

dsf=[];

if ~(exist(filename)==2'),
  SegymatVerbose([mfilename,' : ', filename,' does not exist !'])
  Data=[];SuTraceHeaders=[];SuHeader=[];HeaderInfo=[];
  return
end


% SET THE NEXT FLAG TO 1, IF YOU ARE SURE ALL TRACES HAVE THE SAME LENGT
% THIS WILL GREATLT INCREASE READING WHEN USING 'JUMP'
FixedTraceLength=1;

% NEXT TWO LINES TO ENUSRE THAT VARARGIN CAN BE PASSED TO FUNCTION
if nargin==2
    % CALL USING VARARGIN
    local_nargin=1+length(varargin{1});
    varargin=varargin{1};
else
    % DIRECT CALL
    local_nargin=length(varargin);
end


%
% HERE SHOULD BE A CALL TO A DEFAULT SEGY HEADER !!!!
%

SegyHeader.SegyFormatRevisionNumber=1;    
SegyHeader.DataSampleFormat=5;    
SegyHeader.Rev=GetSegyHeaderBasics;

% TRANSFORM VARARGING INTO PARAMETERS
cargin=1;
while (cargin<local_nargin)
   
   if strcmp(varargin{cargin},'endian')
       cargin=cargin+1;
       eval(['endian=char(varargin{cargin});'])
       if strcmp(endian,'b'),SegymatVerbose(['Reading BIG ENDIAN STYLE']);end
       if strcmp(endian,'l'),SegymatVerbose(['Reading LITTLE ENDIAN STYLE']);end
    end

    if strcmp(varargin{cargin},'DataFormat')
       cargin=cargin+1;
       eval(['DataFormat=char(varargin{cargin});'])
 
       if strcmp(DataFormat,'float32'),
 	       SegyHeader.DataSampleFormat=5; % IEEE
       end

       if strcmp(DataFormat,'int32'),
	       SegyHeader.DataSampleFormat=2; % 4 Byte, two's
                                          % complement integer
       end
       if strcmp(DataFormat,'int16'),
	   SegyHeader.DataSampleFormat=3; % 2 Byte, two's
                                          % complement integer
       end
       if strcmp(DataFormat,'int8'),
	   SegyHeader.DataSampleFormat=8; % 2 Byte, two's
                                          % complement integer
       end
       SegymatVerbose(['DataFormat : ',num2str(SegyHeader.DataSampleFormat)])
    end

    if strcmp(varargin{cargin},'dsf')
      cargin=cargin+1;
      eval(['dsf=',num2str(varargin{cargin}),';']);
      SegymatVerbose(['USING Data Sample Format : dsf=',num2str(dsf)])
    end    
       
    if strcmp(varargin{cargin},'jump')
       cargin=cargin+1;
       eval(['jump=',num2str(varargin{cargin}),';']);
       SegymatVerbose(['JUMP : Read only every ',num2str(jump),'th trace'])
    end

    if strcmp(varargin{cargin},'minmax')
       cargin=cargin+1;
       eval(['header=''',varargin{cargin},''';']);
       cargin=cargin+1;
       eval(['headermin=',num2str(varargin{cargin}),';']);
       cargin=cargin+1;
       eval(['headermax=',num2str(varargin{cargin}),';']);
       SegymatVerbose(['MIN MAX : Using header ',header,' from ',num2str(headermin),' to ',num2str(headermax)])
    end    

    if strcmp(varargin{cargin},'trange')
       cargin=cargin+1;
       eval(['tmin=',num2str(varargin{cargin}),';']);
       cargin=cargin+1;
       eval(['tmax=',num2str(varargin{cargin}),';']);
       SegymatVerbose(['TRANGE : tmin=',num2str(tmin),' tmax=',num2str(tmax)])
    end    
    
    cargin=cargin+1;
    
end



%
% MAYBE DATA FORMATS CAN VARY ?
% but for now we use float32 if nothing else is set
%
if exist('SegyHeader')==0,
  SegyHeader.DataSampleFormat=5; % IEEE
end

if isempty(dsf)==0,
    SegyHeader.DataSampleFormat=dsf; 
end

%
% SU DATA CAN BE EITHER BIG OR SMALL ENDIAN
% DEFAULT IS BIG ENDIAN
%

if exist('endian')==0, 
  segyid = fopen(filename,'r');   % USE LOCAL BUTE ORDER AS DEFAULT
else
  segyid = fopen(filename,'r',endian); 
  endian
end
				    
if ~(exist(filename)==2'),
  SegymatVerbose([mfilename,' : ', filename,' does not exist !'])
  Data=[];SegyTraceHeaders=[];SegyHeader=[];
  return
end




if exist('SkipData','var')==0,
    SkipData=0; % [0] READ ONLY HEADER VALUES, [1] READ IN ALL DATA
end
if SkipData==1, SegymatVerbose(['Not reading data - headers only']), end

% SEGY HEADER FORMAT INFO
Revision=SegyHeader.SegyFormatRevisionNumber;
if Revision>0, Revision=1; end
Format=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).name;
BPS=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).bps;
SegymatVerbose([mfilename,' : ',Format],2)

% GET SIZE OF FILE

fseek(segyid,0,'eof'); DataEnd=ftell(segyid);
fseek(segyid,0,'bof');       % Go to the beginning of the file

%txt=['SegyRevision ',sprintf('%0.4g',Revision),', ',Format,'(',num2str(SegyHeader.DataSampleFormat),')'];
txt=[Format,'(',num2str(SegyHeader.DataSampleFormat),')'];


%hw=waitbar(0,['Reading SU-file -',txt]);
traceinfile=0;
outtrace=0;
existJump=exist('jump');
existHeader=exist('header');
existtmin=exist('tmin');
existtmax=exist('tmax');

tic;
while (~(ftell(segyid)>=DataEnd))
  traceinfile=traceinfile+1;
  usetrace=1;
  
  %if exist('jump')
  if existJump==1;
    if (traceinfile/jump)~=round(traceinfile/jump), usetrace=0; end  
  end
  
  if ((usetrace==0)&(FixedTraceLength==1)&(exist('ns')==1))
    skip=240+(SegyHeader.BytesPerSample/8)*ns;
    fseek(segyid,skip,'cof');
  else
    % Read Trace Header
    SingleSuTraceHeader=GetSegyTraceHeader(segyid);

    ns=SingleSuTraceHeader.ns;

    % IF HEADER MIN MAX HAS BEEN CHOSEN, THEN CHECK THAT TRACE IS GOOD ENOUGH
    if ((existHeader==1)&(usetrace==1))
      headervalue=getfield(SingleSuTraceHeader,header);
       if ((headervalue<headermin)|(headervalue>headermax))
           usetrace=0;
       end
    end

    if usetrace==1;
      % Read Segy Trace Data
      SkipData=0;
    else
      SkipData=1;
    end
     
    SingleSuData.data=GetSegyTraceData(segyid,SingleSuTraceHeader.ns,SegyHeader,SkipData);
  
  end % END READ TRACE OR NO FIXED LENGTH
  
  if (usetrace==1)
    %% IF TIME RANGE IS SPECIFIED, THEN EXTRACT THIS
    if (existtmin)&(existtmax)
          % NEXT LINE SHOULD CONSIDER THAT ns in Trace and Segy Header could vary !!!
          origtrange=[1:1:SingleSuTraceHeader.ns].*SingleSuTraceHeader.dt./1e+6 + SingleSuTraceHeader.DelayRecordingTime./1e+3;
          gooddata=find(origtrange>tmin & origtrange<tmax);
          SingleSuData.data=SingleSuData.data(gooddata);
          % CHECK NEXT LINE TAHT DelatRec... is in micro seconds
          SingleSuTraceHeader.DelayRecordingTime=tmin;
          SingleSuTraceHeader.ns=length(gooddata);
          ns=length(gooddata); %  for use below 
    end
    outtrace=outtrace+1;
    SuTraceHeaders(outtrace)=SingleSuTraceHeader;
    SuData(outtrace).data=SingleSuData.data;
  
    %    keyboard
  end
 
  
  %  waitbar(ftell(segyid)/DataEnd,hw);
end
%close(hw)
SegymatVerbose([mfilename,' : Elapsed time ',num2str(toc)],2);

ns=SuTraceHeaders(1).ns;
dt=SuTraceHeaders(1).dt;
nt=length(SuData);

% DefaultSegyHeader
SuHeader=SegyHeader;
SuHeader.ns=ns;
SuHeader.nsOrig=ns;
SuHeader.dt=dt;
SuHeader.dtOrig=dt;
SuHeader.FixedLengthTraceFlag=1; % CHECK THAT THIS IS TRUE
SuHeader.SegyFormatRevisionNumber=1;
SuHeader.NumberOfExtTextualHeaders=0;


Data=zeros(ns,nt);
try
  for it=1:nt
    Data(:,it)=SuData(it).data;
  end
catch
  SegymatVerbose([mfilename,' Could not collect data in one matrix - check byte order'])
end

%[HeaderInfo]=TraceheaderToInfo(SuTraceHeaders);
return


