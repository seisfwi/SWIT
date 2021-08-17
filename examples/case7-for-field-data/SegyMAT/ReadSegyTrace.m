% ReadSegyTrace
function [Data,STH,SegyHeader]=ReadSegyTrace(filename,traces,SegyHeader);

Data=[];
if nargin<3
    SegyHeader=GetSegyHeader(filename);
end

if nargin<2
    [Data,STH,SegyHeader]=ReadSegyTrace(filename);
end

if (length(traces)>1)
    try
        Data=zeros(SegyHeader.ns,length(traces));
    catch
        keyboard
    end
    for i=1:length(traces)
        [Data(:,i),STH(i)]=ReadSegyTrace(filename,traces(i),SegyHeader);
    end
    return
end


endian='ieee-be'; % Big Endian is default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPEN FILE HANDLE
if exist('endian')==1,
    SegymatVerbose([mfilename,' : ENDIAN : ',endian],1)
    segyid = fopen(filename,'r',endian);
else
    endian='ieee-be';
    SegymatVerbose([mfilename,' : ENDIAN SET TO ',endian],0)
    segyid = fopen(filename,'r','ieee-be');  % ALL DISK FILES ARE IN BIG
end                                        % ENDIAN FORMAT, ACCORDING TO

Revision=SegyHeader.SegyFormatRevisionNumber;
if Revision>0, Revision=1; end

FormatName=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).name;
Format=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).format;
BPS=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).bps;
txt=['SegyRevision ',sprintf('%0.4g',Revision),', ',FormatName,'(',num2str(SegyHeader.DataSampleFormat),')'];


% move to proper location in segyfile
DataStart=3600+3200*SegyHeader.NumberOfExtTextualHeaders;
fseek(segyid,DataStart,'bof');       % Go to the beginning of the file


skip=240+(BPS/8)*SegyHeader.ns;
fseek(segyid,(traces-1)*skip,'cof');
%SegymatVerbose([num2str(traceinfile),' - SKIPPING TRACE ... ',num2str(outtrace)],2)

TraceStart=ftell(segyid);

STH=GetSegyTraceHeader(segyid,TraceStart,Format,SegyHeader.ns,[]);
Data=GetSegyTraceData(segyid,STH.ns,SegyHeader);


fclose(segyid);