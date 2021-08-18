% ReadSegyFast : Reads a SEG Y rev 1 formatted file, without header values (faster than ReadSegy)
%
% Call :
% [Data]=ReadSegyFast(filename);
% and equivalent to :
% [Data]=ReadSegy(filename);
% 
%
% Read only the data of a SegFile - NOT Their headers.
% Much faster than ReadSegy
%
% 'minmax', 'skip'
%


% Implemented using the syntax of the SEG-Y revised format :
% SEGY-Y rev 0, SEG-Y rev 1 as described in 
% http://seg.org/publications/tech-stand/
%
% Extended Textual Header is not yet tested
% If you would like it implemented, please send me an SEGY file with
% that sort of information, as well as a description of the segy file
%
%
% (C) 2001, Thomas Mejer Hansen, thomas@cultpenguin.com
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
function [Data]=ReadSegyFast(filename,varargin);

  mfilename='ReadSegyFast';
  
  if nargin==0,
    SegymatVerbose([mfilename,' : ', filename,' does not exist !'])
    Data=[];SegyTraceHeaders=[];SegyHeader=[];HeaderInfo=[];
    return
  end
  
  % Initialize varargin values
  dsf=[];
  revision=[];
  endian_tight=[];
  tmin=[];tmax=[];
  headermin=[];headermax=[];header=[];
  jump=[];
  SegyHeader=[];

  segyid = fopen(filename,'r','b');   % ALL DISK FILES ARE IN BIG
                                    % ENDIAN FORMAT, ACCORDING TO SEG
                                    % Y rev 1

  fseek(segyid,0,'eof');
  BytesInFile = ftell(segyid);
  fseek(segyid,0,'bof');
  

ninput=nargin;
% NEXT TWO LINES TO ENUSRE THAT VARARGIN CAN BE PASSED TO FUNCTION
if ninput==2
    % CALL USING VARARGIN
    ninput=1+length(varargin{1});
    varargin=varargin{1};
else
    % DIRECT CALL
    ninput=length(varargin);
end


% TRANSFORM VARARGING INTO PARAMETERS
cargin=1;
while (cargin<ninput)
    
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

    if strcmp(varargin{cargin},'revision')
       cargin=cargin+1;
       eval(['revision=',num2str(varargin{cargin}),';']);
       SegymatVerbose(['USING REVISION : rev=',num2str(revision)])
    end    

    if strcmp(varargin{cargin},'dsf')
       cargin=cargin+1;
       eval(['dsf=',num2str(varargin{cargin}),';']);
       SegymatVerbose(['USING Data Sample Format : dsf=',num2str(dsf)])
    end    

    if strcmp(varargin{cargin},'SegyHeader')
       cargin=cargin+1;
       SegyHeader=varargin{cargin};
       SegymatVerbose(['USING LOADED SEGYHEADER'])
    end    

    cargin=cargin+1;
    
end

% SKIP DATA SEEMS OUTDATED, COMMENTING OUT 2003-10-24

%if exist('SkipData')==0,
%    SkipData=0; % [0] READ ONLY HEADER VALUES, [1] READ IN ALL DATA
%end%%
%
%if SkipData==1, SegymatVerbose(['Not reading data - headers only']), end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BINARY HEADERS
if isempty(SegyHeader)==1,
  SegyHeader=GetSegyHeader(segyid);
else
  SegymatVerbose([mfilename,' - Using supplied SegyHeader'])
end


% APPLY CHANGES TO SEGY HEADER IF NEEDE
if isempty(revision)==0, 
  SegyHeader.SegyFormatRevisionNumber=revision; 
  SegymatVerbose([mfilename,' - Manually set SEG-Y revision to ',num2str(revision)])
end
if isempty(dsf)==0, 
  SegyHeader.DataSampleFormat=dsf; 
end



% JUST SOME INFORMATION TO WRITE TO SCREEN :

Revision=SegyHeader.SegyFormatRevisionNumber;
if Revision>0, Revision=1; end
FormatName=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).name;
Format=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).format;  
BPS=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).bps; 
txt=['SegyRevision ',sprintf('%0.4g',Revision),', ',FormatName,'(',num2str(SegyHeader.DataSampleFormat),')'];


ns=SegyHeader.ns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ TEXTURAL FILE HEADER EXTENSION IF NEEDED
%if SegyHeader.NumberOfExtTextualHeaders~=0
%  SegymatVerbose(['---------------------------------------------------'])
%  SegymatVerbose(['extendeD textual file headers are not supported yet'])
%  SegymatVerbose(['They are simply skipped'])
%  SegymatVerbose(['---------------------------------------------------'])
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ DATA
%Segy=fread(segyid,4000,'float32');
fseek(segyid,0,'eof'); DataEnd=ftell(segyid);

DataStart=3600+3200*SegyHeader.NumberOfExtTextualHeaders;
fseek(segyid,DataStart,'bof');       % Go to the beginning of the file

ntraces=(DataEnd-DataStart)./(240+(SegyHeader.ns)*(BPS/8));

SegymatVerbose(['Number of Samples Per Trace=',num2str(SegyHeader.ns)])
SegymatVerbose(['Number of Traces=',num2str(ntraces)])

% Set Format :
Revision=SegyHeader.SegyFormatRevisionNumber;
if Revision>0, Revision=1; end
Format=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).format;  
bitpersample=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).bps;
bytepersample=bitpersample/8;
ntraceheader=240/bytepersample;


Data=fread(segyid,[ns+ntraceheader ntraces],Format);
Data=Data([ntraceheader+1:1:ns+ntraceheader],:);
%SegymatVerbose(['Now=',num2str(ftell(segyid)),' END=',num2str(DataEnd)])

% THE FOLLOWING DOES NOT WORK FOR OCTAVE (uint32)
if (strcmp(Format,'uint32')==1), % IBM FLOATING POINT
        % CONVERT FROM FLOATING POINT
        verbose=1;
        if verbose>1, SegymatVerbose([mfilename,'Converting from IBM, DataFormat :',SegyHeader.DataFormat]); end
        Data=ibm2num(uint32(Data));
end;

existJump=1-isempty(jump);
existHeader=1-isempty(header);
existTmin=1-isempty(tmin);
existTmax=1-isempty(tmax);

if existJump
  usetraces=[jump:jump:ntraces];
  ntraces=length(usetraces);
  Data=Data(:,usetraces);
end



if (existTmin==1)&(existTmax==1)
  SegymatVerbose('TRANGE')
  % NEXT LINE SHOULD CONSIDER THAT ns in Trace and Segy Header could vary !!!
  origtrange=[1:1:SegyHeader.ns].*SegyHeader.dt.*1e-6;
  gooddata=find(origtrange>tmin & origtrange<tmax);
  Data=Data(gooddata,:);
end

