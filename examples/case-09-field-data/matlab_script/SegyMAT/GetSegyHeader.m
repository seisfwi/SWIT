% GetSegyHeader : Reads the segyheader of a SEGY Y formatted file
%
% Call :
%  [SegyHeader]=GetSegyHeader(segyid);
%
%  segyid can be a filehandle or a filename
%  
%
% (C) 2001-2004 Thomas Mejer Hansen, tmh@gfy.ku.dk/thomas@cultpenguin.com
%
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
function [SegyHeader]=GetSegyHeader(segyid,varargin);
  mfilename='GetSegyHeader';
  SegymatVerbose([mfilename,' : Start'],90);
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
    
   
      cargin=cargin+1;
      
  end
  
  
  
if ischar(segyid)
    % if segyid is a string it is interpreted as a filename
    %
    segyid = fopen(segyid,'r','b');   % ALL DISK FILES ARE IN BIG
                                    % ENDIAN FORMAT, ACCORDING TO SEG
                                    % Y rev 1
end
% Basic Segy Header Information.
SegyHeader.Rev=GetSegyHeaderBasics;
fseek(segyid,0,'bof');
SegyHeader.TextualFileHeader=fread(segyid,3200,'uchar');         % 3200
SegyHeader.Job=fread(segyid,1,'int32');                           % 3204 
SegyHeader.Line=fread(segyid,1,'int32');                          % 3208
SegyHeader.Reel=fread(segyid,1,'int32');                          % 3212
SegyHeader.DataTracePerEnsemble=fread(segyid,1,'int16');        % 3214
SegyHeader.AuxiliaryTracePerEnsemble=fread(segyid,1,'uint16');   % 3216
SegyHeader.dt=fread(segyid,1,'uint16');                          % 3218
SegyHeader.dtOrig=fread(segyid,1,'uint16');                      % 3220
SegyHeader.ns=fread(segyid,1,'uint16');                          % 3222
SegyHeader.nsOrig=fread(segyid,1,'uint16');                      % 3224
SegyHeader.DataSampleFormat=fread(segyid,1,'int16');            % 3226
SegyHeader.EnsembleFold=fread(segyid,1,'int16');                
SegyHeader.TraceSorting=fread(segyid,1,'int16');               % 3228
SegyHeader.VerticalSumCode=fread(segyid,1,'int16');            % 3230
SegyHeader.SweepFrequencyStart=fread(segyid,1,'int16');        % 3232
SegyHeader.SweepFrequencyEnd=fread(segyid,1,'int16');          % 3234
SegyHeader.SweepLength=fread(segyid,1,'int16');                % 3236
SegyHeader.SweepType=fread(segyid,1,'int16');                  % 3238
SegyHeader.SweepChannel=fread(segyid,1,'int16');               % 3240
SegyHeader.SweepTaperlengthStart=fread(segyid,1,'int16');               % 3242
SegyHeader.SweepTaperLengthEnd=fread(segyid,1,'int16');               % 3244
SegyHeader.TaperType=fread(segyid,1,'int16');               % 3246
SegyHeader.CorrelatedDataTraces=fread(segyid,1,'int16');               % 3248
SegyHeader.BinaryGain=fread(segyid,1,'int16');               % 3250
SegyHeader.AmplitudeRecoveryMethod=fread(segyid,1,'int16');               % 3252
SegyHeader.MeasurementSystem=fread(segyid,1,'int16');               % 3254
SegyHeader.ImpulseSignalPolarity=fread(segyid,1,'int16');               % 3256
SegyHeader.VibratoryPolarityCode=fread(segyid,1,'int16');               % 3258
% 3261-3500 UNASSIGNED (as 120*2byte integer)
SegyHeader.Unassigned1=fread(segyid,120,'int16');               % 3260
% fseek(segyid,3500,'bof');
SegyHeader.SegyFormatRevisionNumber=fread(segyid,1,'uint16');   % 3500
SegyHeader.FixedLengthTraceFlag=fread(segyid,1,'integer*2');        % 3502
SegyHeader.NumberOfExtTextualHeaders=fread(segyid,1,'uint16');        % 3504
% 3506-3600 UNASSIGNED (as 47*2byte integer = 94 byte)
SegyHeader.Unassigned2=fread(segyid,47,'int16');               % 3506
% READ TEXTURAL FILE HEADER EXTENSION IF NEEDED
%fseek(segyid,3600,'bof');
if (SegyHeader.NumberOfExtTextualHeaders>0),
  SegymatVerbose(['---------------------------------------------------'])
  SegymatVerbose(['extended textual file headers are implemented      '])
  SegymatVerbose(['but have not been tested, since I(tmh@gfy.ku.dk)   '])
  SegymatVerbose(['have had no access to SEGY REV-1 files formatted   '])
  SegymatVerbose(['like this. Please contact me if you have such a    '])
  SegymatVerbose(['file to share   '])
  SegymatVerbose(['---------------------------------------------------'])
  
  txt=sprintf('%d Extended Textual File Headers',SegyHeader.NumberOfExtTextualHeaders);
  SegymatVerbose(txt);
  nChars=3200*SegyHeader.NumberOfExtTextualHeaders;
  SegyHeader.ExtTextualHeaders=fread(segyid,nChars,'schar');        % 3600
else
  SegymatVerbose('NO extended textual file headers');
end
SegyHeader.time=[1:1:SegyHeader.ns].*SegyHeader.dt./1e+6;
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
try
  Format=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).name;
  SegymatVerbose([mfilename,' : SegyRevision ',sprintf('%0.4g',Revision),', ',Format,'(',num2str(SegyHeader.DataSampleFormat),')'])
catch
end
  
if ischar(segyid)
    fclose(segyid);
end
  SegymatVerbose([mfilename,' : End'],90);
  