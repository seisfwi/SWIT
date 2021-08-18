% Sac2Segy : Reads SAC formatted data into a SegyMAT (SGY) structure
%
% CALL :
%   [Data,SegyTraceHeader,SegyHeader]=Sac2Segy(files_in,segyfile_out,varargin)
%
%   files_in : Either a single filename or a strcture of filenames
%            files_in='d1.SAC';
%            or
%            files_in{1}='d1.SAC';
%            files_in{2}='d2.SAC';
%
% Examples :
%   [D,STH,SH]=Sac2Segy('','test.segy','FixedLengthTraceFlag',1);
%              converts all SAC files into one SEGY file (test.segy), using
%              a FixedLengthTraceFlag of 1. This is compatible with mosty
%              any SEGY reader.
%
%   [D,STH,SH]=Sac2Segy('','test.segy','FixedLengthTraceFlag',0);
%              converts all SAC files into one SEGY file (test.segy), using
%              a FixedLengthTraceFlag of 0, allowing varying trace length of SEGY files
%              This is only compatible with revision 1 of the SEGY format.
%
%   [D,STH,SH]=Sac2Segy('file.sac');
%              convert file.sac to file.segy
%
%   [D,STH,SH]=Sac2Segy('file.sac','another_file.segy');
%              convert file.sac to another_file.segy
%
%
%   Force little endian byte format for SAC file:
%   Sac2Segy('file.sac','test.sgy','endian','l');
%
% Relies on sac2mat.m 
%
% Download SAC files from : http://www.iris.edu/hq/ssn/events
% 
function [Data,STH,SegyHeader]=Sac2Segy(files_in,segyfile_out,varargin) 

endian='b'; % Default big endian SAC files

if nargin==0
    d=dir('*.sac');
    if (length(d)==0)
        SegymatVerbose(sprintf('No SAC files found.',mfilename));
        Data=[];STH=[];SegyHeader=[];
        return
    end
    for i=1:length(d); files_in{i}=d(i).name;end
end


if nargin>0
    if isempty(files_in)
        d=dir('*.sac');
        for i=1:length(d); files_in{i}=d(i).name;end
    end
    if (ischar(files_in)&(nargin==1))
        [p,f,e]=fileparts(files_in);
        segyfile_out=[f,'.segy'];
    end
end

ninput=nargin;
% NEXT TWO LINES TO ENUSRE THAT VARARGIN CAN BE PASSED TO FUNCTION
if ninput==3
    % CALL USING VARARGIN
    ninput=2+length(varargin{1});
    varargin=varargin{1};
else
    % DIRECT CALL
    ninput=length(varargin);
end


% TRANSFORM VARARGING INTO PARAMETERS
cargin=1;
while (cargin<ninput)
    SegymatVerbose([mfilename,' - Converting varargin, ',varargin{cargin}],-1)

    if strcmp(varargin{cargin},'FixedLengthTraceFlag')
        cargin=cargin+1;
        eval(['FixedLengthTraceFlag=',num2str(varargin{cargin}),';']);
        SegymatVerbose(['FixedLengthTraceFlag = ',num2str(FixedLengthTraceFlag),'.'])
    end

    if strcmp(varargin{cargin},'endian')
        cargin=cargin+1;
        eval(['endian=''',(varargin{cargin}),''';']);
        SegymatVerbose(['endian = ',endian,'.'])
    end

    cargin=cargin+1;
end

if exist('sac2mat.m','file')~=2
    SegymatVerbose(['sac2mat needs to be in your path'],0)
    SegymatVerbose(['Get it from http://mgstat.sourceforge.net/'],0)
    return
end

try
    [SACdata,SeisData] = sac2mat(files_in,endian);
catch MESS
    if (strfind(MESS.message,'Out of memory'))
        SegymatVerbose(['Out of memory error calling ''sac2mat'', suggesting endian type error'],0)
        SegymatVerbose(['  Try manyally setting the endian type'],0)
        return
    end
end
ntraces=size(SeisData,2);
ns=[SACdata.trcLen];
ns_max=max([SACdata.trcLen]);
data=size(ns_max,ntraces);

% GET dt
% ONE SHOULD MULTIPLY WITH 1e+6 USING SEGY FORMAT
% HOWEEVER SINCE DT IS WRITTEN IN
% UINT16 FORMAT, AND SAC DT IS
% USUALLY VERY HIGH WE MUST CHOOSE
% TO MULTIPLY ONLY WITH 1000.
for i=1:length(ns);
    dt(i)=SACdata(i).times.delta.*1e+3;
end

% -------------------------------------------
% SET UP SegyHeader structure.
% -------------------------------------------


% IF A SPECFIC REVISION HAS BEEN CHOSEN, USE THAT
if exist('revision')==1,
    if revision==1,
        SegyHeader.SegyFormatRevisionNumber=100;
    else
        SegyHeader.SegyFormatRevisionNumber=0;
    end
    SegymatVerbose([mfilename,' :  Using user specified SEG Y revision : ',num2str(revision)],1)
else
    SegyHeader.SegyFormatRevisionNumber=100;
end


% IF A SPECFIC DATA SAMPLING FORMAT HAS BEEN SELECTED USE THAT
if exist('dsf')==1,
    SegyHeader.DataSampleFormat=dsf;
    SegymatVerbose([mfilename,' :  Using user specified Data Sample Format : ',num2str(revision)],1)
else
    SegyHeader.DataSampleFormat=5;
end

if exist('FixedLengthTraceFlag')==1,
    SegyHeader.FixedLengthTraceFlag=FixedLengthTraceFlag;
else
    if ntraces==1,
        SegyHeader.FixedLengthTraceFlag=1;
    else
        if length(unique(ns))==1,
            SegyHeader.FixedLengthTraceFlag=1;
        else;
            SegyHeader.FixedLengthTraceFlag=0;
        end
    end
end
%SegyHeader.FixedLengthTraceFlag=1;



SegyHeader.dt=dt(1);
SegyHeader.dtOrig=dt(1);

if exist('TextualFileHeader'), SegyHeader.TextualFileHeader=TextualFileHeader; end
if exist('Job')==1, SegyHeader.Job=Job; end;
if exist('Line')==1, SegyHeader.Line=Line; end
if exist('Reel')==1, SegyHeader.Reel=Reel; end
if exist('DataTracePerEnsemble')==1, SegyHeader.DataTracePerEnsemble=DataTracePerEnsemble; end
if exist('AuxiliaryTracePerEnsemble')==1, SegyHeader.AuxiliaryTracePerEnsemble=AuxiliaryTracePerEnsemble; end
if exist('ns')==1, SegyHeader.ns=ns(1); end
if exist('nsOrig')==1, SegyHeader.nsOrig=nsOrig(1); end
if exist('EnsembleFold')==1, SegyHeader.EnsembleFold=EnsembleFold;  end
if exist('TraceSorting')==1, SegyHeader.TraceSorting=TraceSorting; end
if exist('VerticalSumCode')==1, SegyHeader.VerticalSumCode=VerticalSumCode; end

% -------------------------------------------
% SETUP SEGY TRACE HEADER
% -------------------------------------------

for i=1:ntraces;

    if i==1,
        STH(i)=InitSegyTraceHeader(ns(i),dt(i));
    else
        STH(i)=STH(1);
    end
    if SegyHeader.FixedLengthTraceFlag==1;
        STH(i).ns=max(ns(i));
    else
        STH(i).ns=ns(i);
    end
    STH(i).TraceSequenceFile=i;

    % EVENT DATA
    STH(i).YearDataRecorded=SACdata(i).event.nzyear;
    STH(i).DayOfYear=SACdata(i).event.nzjday;
    STH(i).HourOfDay=SACdata(i).event.nzhour;
    STH(i).MinuteOfHour=SACdata(i).event.nzmin;
    STH(i).SecondOfMinut=SACdata(i).event.nzsec;

    % TRIMES
    try
        STH(i).dt=dt(i);
    catch
        keyboard
    end

    % STATION DATA
    STH(i).Inline3D=SACdata(i).station.stla;
    STH(i).Crossline3D=SACdata(i).station.stlo;
    STH(i).cdpX=SACdata(i).station.stla;
    STH(i).cdpY=SACdata(i).station.stlo;
    STH(i).ReceiverGroupElevation=SACdata(i).station.stel;
    STH(i).ReceiverDatumElevation=SACdata(i).station.stel;
    %SACdata(i).station.stdp
    %SACdata(i).station.cmpaz
    %SACdata(i).station.cmpinc
    %SACdata(i).station.kstnm
    %SACdata(i).station.kcmpnm
    %SACdata(i).station.knetwk

    Data(:,i)=SeisData(:,i);
    %if ns(i)<max(ns)
    %    Data((ns(i)+1):max(ns),i)=NaN;
    %end

end

if SegyHeader.FixedLengthTraceFlag==1;
    ins=find(ns==max(ns));ins=ins(1);
    SegyHeader.ns=ns(ins);
    SegyHeader.dt=dt(ins);
    for i=1:ntraces
        STH(i).ns=ns(ins);
        STH(i).dt=dt(ins);
    end
end

% WRITE SEGY STRUCTURE IF REQUESTED
if ((nargin>1)|(exist('segyfile_out','var')))
    SegyHeader=WriteSegyStructure(segyfile_out,SegyHeader,STH,Data);
end