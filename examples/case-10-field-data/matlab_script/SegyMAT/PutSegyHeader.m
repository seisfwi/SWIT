% PutSegyHeader : Writes SEG-Y header to disk.
% PutSegyHeader(segyid,SegyHeader)
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
function SegyHeader=PutSegyHeader(segyid,SegyHeader);

NA=0;

if isfield(SegyHeader,'Rev')==0
  SegyHeader.Rev=GetSegyHeaderBasics;
end

if (isfield(SegyHeader,'DataSampleFormat')==0)&(isfield(SegyHeader,'SegyFormatRevisionNumber')==0)
  SegyHeader.DataSampleFormat=5; % '5'->4-byte IEEE floating point 
  SegyHeader.SegyFormatRevisionNumber=100;
  SegymatVerbose([mfilename,' : Using datasample format : ',SegyHeader.Rev(2).DataSampleFormat(SegyHeader.DataSampleFormat).name],1)
  SegymatVerbose([mfilename,' : Using SEG Y revision : ',num2str(floor(SegyHeader.SegyFormatRevisionNumber/100))],1)
end

% SET Revsion number. Always use revision 1, if not otherwise specified.
if (isfield(SegyHeader,'SegyFormatRevisionNumber')==0)
  SegyHeader.SegyFormatRevisionNumber=100;
  SegymatVerbose([mfilename,' : Using SEG Y revision : ',num2str(floor(SegyHeader.SegyFormatRevisionNumber/100))],1)
end

% SET DATA SAMPLE FORMAT IF NOT SPECIFIED
% Rev0->IBM FLOATING POINT
% Rev1->IEEE FLOATING POINT
if (isfield(SegyHeader,'DataSampleFormat')==0)
  if SegyHeader.SegyFormatRevisionNumber==0
    SegyHeader.DataSampleFormat=1; % '1'->IBM floating point
    SegymatVerbose([mfilename,' : Using datasample format : ',SegyHeader.Rev(1).DataSampleFormat(SegyHeader.DataSampleFormat).name],1)
  else
    SegyHeader.DataSampleFormat=5; % '5'->4-byte IEEE floating point 
    SegymatVerbose([mfilename,' : Using datasample format : ',SegyHeader.Rev(2).DataSampleFormat(SegyHeader.DataSampleFormat).name],1)
  end
end

if isfield(SegyHeader,'TextualFileHeader')
  % MAKE SURE LENGTH IF Text.. is 3200
  if length(SegyHeader.TextualFileHeader)<3200
    dummyTXT=SegyHeader.TextualFileHeader;
    SegyHeader.TextualFileHeader=zeros(3200,0);
    SegyHeader.TextualFileHeader(1:length(dummyTXT))=dummyTXT;
  end  
else
  SegyHeader.TextualFileHeader=sprintf('%3200s','SEGY READER (tmh@gfy.ku.dk)');
end


if ~isfield(SegyHeader,'Job'), SegyHeader.Job=NA; end
if ~isfield(SegyHeader,'Line'), SegyHeader.Line=NA; end
if ~isfield(SegyHeader,'Reel'), SegyHeader.Reel=NA; end
if ~isfield(SegyHeader,'DataTracePerEnsemble'), SegyHeader.DataTracePerEnsemble=NA;end
if ~isfield(SegyHeader,'AuxiliaryTracePerEnsemble'), SegyHeader.AuxiliaryTracePerEnsemble=0;end
if ~isfield(SegyHeader,'dt'), SegyHeader.dt=4000; end
if ~isfield(SegyHeader,'dtOrig'), SegyHeader.dtOrig=NA; end
if ~isfield(SegyHeader,'ns'), SegyHeader.ns=NA; end
if ~isfield(SegyHeader,'nsOrig'), SegyHeader.nsOrig=NA; end
if ~isfield(SegyHeader,'DataSampleFormat'), SegyHeader.DataSampleFormat=5; end % '5'->4-byte IEEE floating point 
if ~isfield(SegyHeader,'EnsembleFold'), SegyHeader.EnsembleFold=NA; end
if ~isfield(SegyHeader,'TraceSorting'), SegyHeader.TraceSorting=NA; end
if ~isfield(SegyHeader,'VerticalSumCode'), SegyHeader.VerticalSumCode=NA; end
if ~isfield(SegyHeader,'SweepFrequencyStart');SegyHeader.SweepFrequencyStart=NA; end
if ~isfield(SegyHeader,'SweepFrequencyEnd');SegyHeader.SweepFrequencyEnd=NA; end
if ~isfield(SegyHeader,'SweepLength');SegyHeader.SweepLength=NA; end
if ~isfield(SegyHeader,'SweepType');SegyHeader.SweepType=NA; end
if ~isfield(SegyHeader,'SweepChannel');SegyHeader.SweepChannel=NA; end
if ~isfield(SegyHeader,'SweepTaperlengthStart');SegyHeader.SweepTaperlengthStart=NA; end
if ~isfield(SegyHeader,'SweepTaperLengthEnd');SegyHeader.SweepTaperLengthEnd=NA; end
if ~isfield(SegyHeader,'TaperType');SegyHeader.TaperType=NA; end
if ~isfield(SegyHeader,'CorrelatedDataTraces');SegyHeader.CorrelatedDataTraces=NA; end
if ~isfield(SegyHeader,'BinaryGain');SegyHeader.BinaryGain=NA; end
if ~isfield(SegyHeader,'AmplitudeRecoveryMethod');SegyHeader.AmplitudeRecoveryMethod=NA; end
if ~isfield(SegyHeader,'MeasurementSystem');SegyHeader.MeasurementSystem=1; end %1-Meters, 2-Feet
if ~isfield(SegyHeader,'ImpulseSignalPolarity');SegyHeader.ImpulseSignalPolarity=NA; end
if ~isfield(SegyHeader,'VibratoryPolarityCode');SegyHeader.VibratoryPolarityCode=NA; end
% 3261-3500 UNASSIGNED
if ~isfield(SegyHeader,'Unassigned1');SegyHeader.Unassigned1=NA*ones(1,120); end

%if ~isfield(SegyHeader,'SegyFormatRevisionNumber'); SegyHeader.SegyFormatRevisionNumber=100; end % ACCORDING TO SEGY REV1 draft 6
if ~isfield(SegyHeader,'FixedLengthTraceFlag');SegyHeader.FixedLengthTraceFlag=1; end % ALL TRACES HAVE THE SAME dt AND ns
if ~isfield(SegyHeader,'NumberOfExtTextualHeaders'); SegyHeader.NumberOfExtTextualHeaders=0; end % WE DO NOT YET USE TEXTURAL HEADERS as can be done according to draft 6
% 3506-3600 UNASSIGNED
if ~isfield(SegyHeader,'Unassigned2');SegyHeader.Unassigned2=NA*ones(1,47); end

% IF EXTENDED TEXTUAL FILE HEADERS EXISTS THEN WRITE THEM TO DISK


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE TO DISK
%


fseek(segyid,0,'bof');

fwrite(segyid,SegyHeader.TextualFileHeader(1:3200),'uchar');         % 1-3200
fwrite(segyid,SegyHeader.Job,'int32');                           % 3204 

fwrite(segyid,SegyHeader.Line,'int32');                          % 3208
fwrite(segyid,SegyHeader.Reel,'int32');                          % 3212
fwrite(segyid,SegyHeader.DataTracePerEnsemble,'int16');        % 3214
fwrite(segyid,SegyHeader.AuxiliaryTracePerEnsemble,'int16');   % 3216
fwrite(segyid,SegyHeader.dt,'uint16');                          % 3218
fwrite(segyid,SegyHeader.dtOrig,'uint16');                      % 3220
fwrite(segyid,SegyHeader.ns,'uint16');                          % 3222
fwrite(segyid,SegyHeader.nsOrig,'uint16');                      % 3224
fwrite(segyid,SegyHeader.DataSampleFormat,'int16');            % 3226
fwrite(segyid,SegyHeader.EnsembleFold,'int16');
fwrite(segyid,SegyHeader.TraceSorting,'int16');               % 3228
fwrite(segyid,SegyHeader.VerticalSumCode,'int16');            % 3230

fwrite(segyid,SegyHeader.SweepFrequencyStart,'int16');        % 3232
fwrite(segyid,SegyHeader.SweepFrequencyEnd,'int16');          % 3234
fwrite(segyid,SegyHeader.SweepLength,'int16');                % 3236
fwrite(segyid,SegyHeader.SweepType,'int16');                  % 3238
fwrite(segyid,SegyHeader.SweepChannel,'int16');               % 3240
fwrite(segyid,SegyHeader.SweepTaperlengthStart,'int16');               % 3242
fwrite(segyid,SegyHeader.SweepTaperLengthEnd,'int16');               % 3244
fwrite(segyid,SegyHeader.TaperType,'int16');               % 3246

fwrite(segyid,SegyHeader.CorrelatedDataTraces,'int16');               % 3248

fwrite(segyid,SegyHeader.BinaryGain,'int16');               % 3250
fwrite(segyid,SegyHeader.AmplitudeRecoveryMethod,'int16');               % 3252

fwrite(segyid,SegyHeader.MeasurementSystem,'int16');               % 3254

fwrite(segyid,SegyHeader.ImpulseSignalPolarity,'int16');               % 3256
fwrite(segyid,SegyHeader.VibratoryPolarityCode,'int16');               % 3258

% 3261-3500 UNASSIGNED1 => (120 int32 = 240 bytes)
fwrite(segyid,SegyHeader.Unassigned1,'int16');               % 3260
%SegymatVerbose(ftell(segyid),10);
%fseek(segyid,3500,'bof');
fwrite(segyid,SegyHeader.SegyFormatRevisionNumber,'uint16');   % 3500
fwrite(segyid,SegyHeader.FixedLengthTraceFlag,'int16');        % 3502
fwrite(segyid,SegyHeader.NumberOfExtTextualHeaders,'uint16');        % 3504

% 3506-3600 UNASSIGNED2 => 94/2=47 int16
fwrite(segyid,SegyHeader.Unassigned2,'int16');               % 3506


%
% 
%
if SegyHeader.NumberOfExtTextualHeaders>0
  n=SegyHeader.NumberOfExtTextualHeaders;
  SegymatVerbose(sprintf('Writing %d Extended Textual File Headers',n));
  if isfield(SegyHeader,'ExtTextualHeader')
    % MAKE SURE LENGTH IF Text.. is 3200*n
    if length(SegyHeader.ExtTextualHeader)<(3200*n)
      dummyTXT=SegyHeader.ExtTextualHeader;
      SegyHeader.ExtTextualHeader=zeros(n*3200,0);
      SegyHeader.ExtTextualHeader(1:length(dummyTXT))=dummyTXT;
    end  
  else
    SegyHeader.ExtTextualHeader=sprintf('%3200s','SEGY READER (tmh@gfy.ku.dk)');
    for i=2:n,
      SegyHeader.ExtTextualHeader=[SegyHeader.ExtTextualHeader,sprintf('%3200s','SEGY READER (tmh@gfy.ku.dk)')];
    end
  end
  % WRITE TEXTUAL FILE HEADER
  fwrite(segyid,SegyHeader.ExtTextualHeader(1:(n*3200)),'uchar');         % 1:(n*3200)
end




Revision=SegyHeader.SegyFormatRevisionNumber;
if Revision>0, Revision=1; end
Format=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).name;
SegymatVerbose([mfilename,' - SegyRevision ',sprintf('%0.4g',Revision),', ',Format],2)

