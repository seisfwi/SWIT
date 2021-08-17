% ReadSegyConstantTraceLength : Reads a SEG Y rev 0/1 formatted file
%                               Assumes CONSTANT TRACE LENGTH
%                               which allows much faster code
%
% Call :
% [Data,SegyTraceHeaders,SegyHeader]=ReadSegyConstantTraceLength(filename);
%
% To read time slice 0.5<t<5 :
% [Data,SegyTraceHeaders,SegyHeader]=ReadSegy(filename,'trange',.5,3);
% Skip every 5th trace :
% [Data,SegyTraceHeaders,SegyHeader]=ReadSegy(filename,'jump',5);
% Read data in a CDP header range : 5000<cdp<5800 :
% (change cdp to any other valid TraceHeader value)
% [Data,SegyTraceHeaders,SegyHeader]=ReadSegy(filename,'minmax','cdp',5000,5800);
% Use several minmax entries
% [Data,SegyTraceHeaders,SegyHeader]=ReadSegy(filename,'minmax','cdp',5000,5800,'minmax','SourceX',10,20);
%
%
% Read from trace 13-18:
% [Data,SegyTraceHeaders,SegyHeader]=ReadSegy(filename,'trace',13:18);
% Read from trace 13-18 and 100-130:
% [Data,SegyTraceHeaders,SegyHeader]=ReadSegy(filename,'trace',[13:18,100:130]);
%
% SEG-Y format revision number can be '0' (1975) or 
% '100' (similar to '1') (2002).
% By default the SEG-Y format revision number is read in the 
% binary header, but this can be overruled using :
% [Data,SegyTraceHeaders,SegyHeader]=ReadSegy(filename,'revision',0);
%
% Read using a specific Data Sample Format :
% Rev 0, IBM FLOATING POINT
% [Data,SegyTraceHeaders,SegyHeader]=ReadSegy(filename,'revision',0,'dsf',1);
% Rev 1, IEEE FLOATING POINT
% [Data,SegyTraceHeaders,SegyHeader]=ReadSegy(filename,'revision',1,'dsf',5);
%
% A SegyHeader can be forced on the SEG-Y file using :
% [Data,SegyTraceHeaders,SegyHeader]=ReadSegy(filename,'SegyHeader',SegyHeader);
% The SegyHeader can be obtain by GetSegyHeader(segyfilename), and
% then edited.
%
% To read using little endian :
% [Data,SegyTraceHeaders,SegyHeader]=ReadSegy(filename,'endian','l');
%
% Combine any combination of the above
% [Data,SegyTraceHeaders,SegyHeader]=ReadSegy(filename,'jump',1,'minmax','cdp',5300,5400);
%
%
% Plot the data using e.g. 
% imagesc([SegyTraceHeaders.cdp],SegyHeader.time,Data);
% wiggle([SegyTraceHeaders.TraceNumber],SegyHeader.time,Data);
%
% (C) 2003-2004, Thomas Mejer Hansen, tmh@gfy.ku.dk
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
% 
%
% 0.1 : INitial Release
% 0.2 : Added SkipData var, to skip reading of data.
% 0.3 : May 01, 2002 
%       Added ability to read in ever 'jump' traces.
%       Added ability to read in time range.
%       Added abiliy to read in header range (ex. mincdp to maxcdp).
%

% TODO : WHEN READING ONLY PART OF DATASET MAKE SURE TO ADJUST THE SEGY
% HEADER ACCORDINGLY !!!!!!

function [Data,SegyTraceHeaders,SegyHeader]=ReadSegy(filename,varargin);


  if isoctave
    doWaitBar=0; % [1] show progress bar gui
    mfilename='ReadSegy';
  else
    doWaitBar=1;
    mfilename='ReadSegy';
  end
    
  dsf=[];
  revision=[];
  endian_tight=[];
  tmin=[];tmax=[];
  headermin=[];headermax=[];header=[];
  jump=[];
  SkipData=[];
  %tracestart=[];  traceend=[];
  trace=[];
  
  SegymatVerbose([mfilename,' : reading ',filename])
  
  if ~(exist(filename)==2'),
    SegymatVerbose([mfilename,' : ', filename,' does not exist !'])
    Data=[];SegyTraceHeaders=[];SegyHeader=[];HeaderInfo=[];
    return
  end
  
  % IF ONLY 'filename', AND one outpuet HAS BEEN
  % SPECIFIED AS IN/OUTPUT, THEN USE THE FAST
  % ALGORITHM FOR READING.
  if (nargin==1)&(nargout==1)
    [Data]=ReadSegyFast(filename);
    return
  end
  
  SegymatVerbose([mfilename,' - Checking Varargin'],90)
  
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
  
  
  nmm=0;
  
  % TRANSFORM VARARGING INTO PARAMETERS
  cargin=1;
  while (cargin<ninput)
    
    SegymatVerbose([mfilename,' - Converting varargin, ',num2str(cargin)],90)
    
    if strcmp(varargin{cargin},'jump')
      cargin=cargin+1;
      eval(['jump=',num2str(varargin{cargin}),';']);
      SegymatVerbose(['JUMP : Read only every ',num2str(jump),'th trace'])
    end

    
        
    if strcmp(varargin{cargin},'trace')
      cargin=cargin+1;
      trace=varargin{cargin};
      SegymatVerbose(sprintf('trace: reading %d from trace %d to %d',length(trace),min(trace),max(trace)));

    end

    
    if strcmp(varargin{cargin},'minmax')
      nmm=nmm+1;
      cargin=cargin+1;
      eval(sprintf('header{%d}=''%s'';',nmm,varargin{cargin}))
      cargin=cargin+1;
      eval(sprintf('headermin(%d)=%f;',nmm,varargin{cargin}));
      cargin=cargin+1;
      eval(sprintf('headermax(%d)=%f;',nmm,varargin{cargin}));

      SegymatVerbose(['MIN MAX : Using header ',header{nmm},' from ',num2str(headermin(nmm)),' to ',num2str(headermax(nmm))])
    end    
    
    if strcmp(varargin{cargin},'trange')
      cargin=cargin+1;
      eval(['tmin=',num2str(varargin{cargin}),';']);
      cargin=cargin+1;
      eval(['tmax=',num2str(varargin{cargin}),';']);
      SegymatVerbose(['TRANGE : tmin=',num2str(tmin),' tmax=',num2str(tmax)])
    end    
    
    % ENDIAN FORMAT 
    endian='ieee-be'; % Big Endian is default
    if strcmp(varargin{cargin},'endian')
      cargin=cargin+1;
      eval(['endian_tight=varargin{cargin};'])
      if endian_tight=='l',
        SegymatVerbose(['USING LITTLE ENDIAN TYPE'])
        endian='ieee-le';
      else
        SegymatVerbose(['USING BIG ENDIAN TYPE'])
      end
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
    
    if strcmp(varargin{cargin},'SkipData')
      cargin=cargin+1;
      eval(['SkipData=',num2str(varargin{cargin}),';']);
      SegymatVerbose(['SKIPPING DATA - READING ONLY HEADERS'])
    end    
    
    if strcmp(varargin{cargin},'SegyHeader')
      cargin=cargin+1;
      SegyHeader=varargin{cargin};
      SegymatVerbose(['USING LOADED SEGYHEADER'])
    end    
    
    cargin=cargin+1;
    
  end
  
  
  if isempty(SkipData)==1,
    SegymatVerbose([mfilename,' : Skip data is not set (defautls to 0)'],90)
    SkipData=0; % [0] READ ONLY HEADER VALUES, [1] READ IN ALL DATA
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % OPEN FILE HANDLE
  
  if exist('endian')==1,
    SegymatVerbose([mfilename,' : ENDIAN : ',endian],90)
    segyid = fopen(filename,'r',endian);   
  else
    endian='ieee-be';
    SegymatVerbose([mfilename,' : ENDIAN SET TO ',endian],90)
    segyid = fopen(filename,'r','ieee-be');  % ALL DISK FILES ARE IN BIG
  end                                        % ENDIAN FORMAT, ACCORDING TO 
                                             % SEGY Y rev 1
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % BINARY HEADERS
  if exist('SegyHeader')==0
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
  
  if (SegyHeader.DataSampleFormat>length(SegyHeader.Rev(Revision+1).DataSampleFormat));
    SegymatVerbose([mfilename,' : WARNING : YOU HAVE SELECTED (OR THE FILE IS FORMATTED SUCH THAT) A DATASAMPLE FORMAT THAT IS NOT DEFINED. \nREMEBER IEEE IS NOT SPECIFIED IN THE SEGY REV0 STANDARD !'])
    
    if (Revision==0)
      SegymatVerbose([mfilename,' : TRYING TO USE REVISION 1 AS OPPOSED TO REVISION 0'])
      Revision=1;
      
      if (SegyHeader.DataSampleFormat>length(SegyHeader.Rev(Revision+1).DataSampleFormat));
        SegymatVerbose([mfilename,' : FATAL ERROR : STILL THE DATASAMPLE FORMAT IS NOT SUPPRTED - EXITING (Report error to tmh@gfy.ku.dk)'])
      else 
        SegymatVerbose([mfilename,' : APPARENT SUCCES CHANING FROM Revision 0 to 1 - Continuing'])
        SegyHeader.SegyFormatRevisionNumber=1; % FORCING REVISION TO BE 1 !!!
      end
      
    end
    
  end
  
  FormatName=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).name;
  Format=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).format;  
  BPS=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).bps; 
  txt=['SegyRevision ',sprintf('%0.4g',Revision),', ',FormatName,'(',num2str(SegyHeader.DataSampleFormat),')'];
  
  
  ns=SegyHeader.ns;
  % YOU CAN FORCE FixedLengthTraceFlag=1;
  % This will make the code much faster (especially when using
  % the 'jump' option) but reading data with varying trace lengt will fail.
  % It is here since many old data sets with Constant trace length 
  % has FixedLengthTraceFlag=0;
  % 
  % As of version 1.01 this has been enable by default.
  % Change the variable below to '0' if you do not want this behaviour
  %
  SegyHeader.FixedLengthTraceFlag=1;
  
  
  SegymatVerbose([mfilename,' : Reading Data'],90);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% READ DATA
  %Segy=fread(segyid,4000,'float32');
  fseek(segyid,0,'eof'); DataEnd=ftell(segyid);
  fseek(segyid,0,'eof'); DataEnd=ftell(segyid);
  
  DataStart=3600+3200*SegyHeader.NumberOfExtTextualHeaders;
  fseek(segyid,DataStart,'bof');       % Go to the beginning of the file
  fseek(segyid,DataStart,'bof');       % Go to the beginning of the file
    
  ntraces=(DataEnd-DataStart)./(240+(SegyHeader.ns)*(BPS/8));
  
  SegymatVerbose(['Number of Samples Per Trace=',num2str(SegyHeader.ns)])
  SegymatVerbose(['Number of Traces=',num2str(ntraces)])
  if (ntraces~=round(ntraces)) 
    SegymatVerbose(['Trace lengths seems to vary. trying to read the file anyway'])
  end
  
  existJump=~isempty(jump);
  existHeader=~isempty(header);
  existTmin=~isempty(tmin);
  existTmax=~isempty(tmax);
  existTrace=~isempty(trace);
  
 

  dwaitbar=1;
  if DataEnd>1e+6, dwaitbar=10;  end
  if DataEnd>1e+7, dwaitbar=50;  end
  if DataEnd>1e+8, dwaitbar=200;  end
  traceinfile=0;  
  outtrace=0;
  tic;
  toc_old=toc;
  
  % FIND TRACES TO USE
  getTrace=1:1:ntraces;

  if existTrace
      %      getTrace=tracestart:1:traceend;
      getTrace=trace;
  end

  if nmm>0
      allTraces=ones(1,ntraces);
      for imm=1:nmm
          val=ReadSegyTraceHeaderValue(filename,'key',header{imm});
          badTraces=find((val<headermin(imm))|(val>headermax(imm)));
          allTraces(badTraces)=0;
      end
      getTrace=find(allTraces);
  end
  
  if existJump==1
      getTrace=getTrace(jump:jump:length(getTrace));
  end

  
  out_ntraces=length(getTrace);


  
  
  % LOOP OVER TRACES
  
  if doWaitBar==1;
      hw=waitbar(0,['Reading Segy - ',txt]);
  end

  
  %while (~(ftell(segyid)>=DataEnd))
  for traceinfile=1:length(getTrace);
    
    usetrace=1; % DEFAULT USING TRACE WHEN [1].
    
    ishow=500;
    if ((traceinfile/ishow)==round(traceinfile/ishow)), 
      PerTrace=(toc-toc_old)/ishow;
      TimeLeft=(out_ntraces-traceinfile)*PerTrace;
      txt=sprintf('Reading trace %d/%d, (%5.0fs left)',traceinfile,out_ntraces,TimeLeft);
      toc_old=toc;
      SegymatVerbose(txt)
    end
    TraceStart=ftell(segyid);

    TraceStart=DataStart+(getTrace(traceinfile)-1)*(SegyHeader.ns*BPS/8+240);

    
    %    if ((usetrace==0)&(SegyHeader.FixedLengthTraceFlag==1)),
    %   % SKIP FORWARD IN FILE'
    %   skip=240+(BPS/8)*SegyHeader.ns;
    %   fseek(segyid,skip,'cof');
    %   %SegymatVerbose([num2str(traceinfile),' - SKIPPING TRACE ... ',num2str(outtrace)])
    % else
      SingleSegyTraceHeaders=GetSegyTraceHeader(segyid,TraceStart,Format,SegyHeader.ns,[]);
      SingleSegyData.data=GetSegyTraceData(segyid,SingleSegyTraceHeaders.ns,SegyHeader);

      if SingleSegyTraceHeaders.TraceNumber<1
        SingleSegyTraceHeaders.TraceNumber=traceinfile;
        SegymatVerbose(sprintf('TraceNumber malformatetd. Setting TraceNumber=%d',traceinfile),10);
      end
      
      SegymatVerbose(sprintf('ns=%d, Trace in line : %d, Trace in file : %d, ns=%10.5f dt=%10.5f',SingleSegyTraceHeaders.ns,SingleSegyTraceHeaders.TraceSequenceLine,SingleSegyTraceHeaders.TraceSequenceFile,SingleSegyTraceHeaders.ns,SingleSegyTraceHeaders.dt),10)

      
      %end
    
    
    %    % IF HEADER MIN MAX HAS BEEN CHOSEN, THEN CHECK THAT TRACE IS GOOD ENOUGH
    % if ((existHeader==1)&(usetrace==1))
    %   headervalue=getfield(SingleSegyTraceHeaders,header); 
    %   if ((headervalue<headermin)|(headervalue>headermax))
    %     usetrace=0;
    %   end
    % end
    
    % USE THIS TRACE IF usetrace=1 !!
    if usetrace==1,
      %% IF TIME RANGE IS SPECIFIED, THEN EXTRACT THIS
      if (existTmin==1)&(existTmax==1)
        % NEXT LINE SHOULD CONSIDER THAT ns in Trace and Segy Header could vary !!!
        origtrange=[1:1:SegyHeader.ns].*SegyHeader.dt.*1e-6+SingleSegyTraceHeaders.DelayRecordingTime.*1e-3;
        gooddata=find(origtrange>tmin & origtrange<tmax);
        SingleSegyData.data=SingleSegyData.data(gooddata);
        % CHECK NEXT LINE TAHT DelatRec... is in micro seconds
        SingleSegyTraceHeaders.DelayRecordingTime=tmin;
        SingleSegyTraceHeaders.ns=length(gooddata);
        ns=length(gooddata); %  for use below 
      end
      
      outtrace=outtrace+1;
      
      if (outtrace==1),
        % Preallocate RAM
        SegymatVerbose(sprintf('Pre allocating RAM ntraces=%d out_traces=%d',ntraces,out_ntraces));
        SegyData=repmat(struct('data',zeros(ns,1)),1,out_ntraces);
        SegyTraceHeaders=repmat(SingleSegyTraceHeaders,1,out_ntraces);
      end
      
      SegyData(outtrace).data=SingleSegyData.data;
      SegyTraceHeaders(outtrace)=SingleSegyTraceHeaders;
      
      if doWaitBar==1,
        if ((outtrace/dwaitbar)==round(outtrace/dwaitbar))
          waitbar(ftell(segyid)/DataEnd,hw);
        end
      end
      
    end
    
  end
  if doWaitBar==1
    close(hw);
  end
  SegymatVerbose([mfilename,' : Elapsed time ',num2str(toc)]);
  t=outtrace;
  
  % CHANGE SEGY HEADER IF TIME RANGE WAS SET
  if (existTmin==1)&(existTmax==1)
    SegyHeader.ns=ns;
    SegyHeader.time=[1:1:SegyHeader.ns].*SegyHeader.dt./1e+6 + SegyTraceHeaders(1).DelayRecordingTime./1e+3;  
  end
  
  % Make sure that only read SegyTraceHEaders are returned
  if outtrace~=out_ntraces
    SegyTraceHeaders=SegyTraceHeaders(1:outtrace);
  end
  
  SegyHeader.time=[1:1:SegyHeader.ns].*SegyHeader.dt./1e+6 + SegyTraceHeaders(1).DelayRecordingTime./1e+3;  

  
  
  % MOVE DATA from SegyData.data to a regular variable
  if SkipData==1,
    Data=[];
  else
    
    try
      Data=[SegyData(1:outtrace).data];
    catch

    
      Data=zeros(ns,outtrace);
      for i=1:outtrace
        try
          Data(:,i)=SegyData(i).data;
        catch
          errmsg=lasterr;
          if isempty(SegyData(i).data)
            errmsg='Empty data in trace';
          elseif (strfind(errmsg, 'In an assignment  A(:,matrix) = B, the number of rows in A and B'))            
            nns=length(SegyData(i).data);
            if nns<ns
              errmsg='Length of trace varies - padding with zeros';
              Data(1:nns,i)=SegyData(i).data;
            else
              errmsg='Length of trace varies - truncating';
              Data(:,i)=SegyData(i).data(1:ns);
            end
          end
          SegymatVerbose(sprintf('Had a problem at trace %d : %s',i,errmsg))
        end
      end 
    end
  end
  
