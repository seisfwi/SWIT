% WriteSegyTraceHeaderValue : Write trace header valaue at specific location
%  
% Call:
%
%    % Update all trace header values starting at position 72, in integer32
%    % format, to the value 30
%    data=30;
%    WriteSegyTraceHeaderValue(filename,data,'pos',72,'precision','int32',);
%
%    % Update all trace header values starting at position 72, in integer32
%    % format, to the values in array 'data'
%    ntraces=311;
%    data=[1:1:311]*10;
%    WriteSegyTraceHeaderValue(filename,data,'pos',72,'precision','int32');
%    d_header=ReadSegyTraceHeaderValue(filename,'pos',72,'precision','int32');
%    
%    % Update the 'cdp' TraceHeader value: 
%    cdp=ReadSegyTraceHeaderValue(file,'key','cdp');  % READ CDP
%    cdp=cdp+10;                                      % change CDP 
%    WriteSegyTraceHeaderValue(file,cdp,'key','cdp'); % UPDATE CDP
%
%    Call 'TraceHeaderDef(1)' to see a list of TraceHeader 'key' names
% 
% See also ReadSegyTraceHeaderValue, PutSegyTraceHeader, TraceHeaderDef
%

function dread=WriteSegyTraceHeaderValue(filename,data,varargin);%pos,type)

pos=0;
precision='int32';

ninput=length(varargin);
% TRANSFORM VARARGING INTO PARAMETERS
cargin=1;
while (cargin<ninput)
    % ENDIAN FORMAT 
    endian='ieee-be'; % Big Endian is default
    if cargin==5, keyboard; end
    if strcmp(varargin{cargin},'endian')
       cargin=cargin+1;
       eval(['endian_tight=varargin{cargin};'])
       if endian_tight=='l',
         disp(['USING LITTLE ENDIAN TYPE'])
         endian='ieee-le';
       else
         disp(['USING BIG ENDIAN TYPE'])
       end
    end    

  
    if strcmp(varargin{cargin},'pos')
       cargin=cargin+1;
       eval(['pos=',num2str(varargin{cargin}),';']);
       disp(['Reading at header postision : pos=',num2str(pos)])
    end    

    if strcmp(varargin{cargin},'precision')
       cargin=cargin+1;
       eval(['precision=''',varargin{cargin},''';']);
       disp(['precision : ',precision])
    end    

    if strcmp(varargin{cargin},'key')
        cargin=cargin+1;
        eval(['key=''',varargin{cargin},''';']);
        
        STH=TraceHeaderDef;
        try
            pos=STH.(key).pos;
            precision=STH.(key).precision;
        catch
            disp(sprintf('Trace Header Value %s not defined',key))
            pos=0; precision='int32';
            hval=[];
            return
        end
        disp(sprintf('key=%s, startpos=%d, precision=%s ',key,pos+1,precision))
    end
        
        cargin=cargin+1;
end


if exist('endian')==1,
  segyid = fopen(filename,'r+',endian);   
else
  segyid = fopen(filename,'r+','ieee-be');  % ALL DISK FILES ARE IN BIG
end                                        % ENDIAN FORMAT, ACCORDING TO 



[SegyHeader]=ReadSegyHeader(filename);

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


fseek(segyid,0,'eof'); DataEnd=ftell(segyid);

DataStart=3600+3200*SegyHeader.NumberOfExtTextualHeaders;
fseek(segyid,DataStart,'bof');       % Go to the beginning of the file

ntraces=(DataEnd-DataStart)./(240+(SegyHeader.ns)*(BPS/8));

hval=zeros(1,SegyHeader.ntraces);
for itrace=1:SegyHeader.ntraces
    if ((itrace/1000)==round(itrace/1000))
	  progress_txt(itrace,ntraces,'Trace #')
 	end	

	%GOTO START OF TRACE HEADER
 	skip=DataStart+(itrace-1)*(240+(BPS/8)*SegyHeader.ns);
	fseek(segyid,skip,'bof');
    fseek(segyid,pos,'cof');
    if length(data)==1;
        d_trace=data;
    else
        d_trace=data(itrace);
    end
	fwrite(segyid,d_trace,precision);

end	

if nargout>0
    dread=ReadSegyTraceHeaderValue(filename,varargin);
end

fclose(segyid);