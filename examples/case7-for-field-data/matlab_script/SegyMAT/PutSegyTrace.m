% PutSegyTrace(segyid,tracedata,SegyTraceHeader,TraceStart);
%  
% Write a SegyTrace to a filehandle 'segyid'
%
% (C) 2001-2004, Thomas Mejer Hansen, tmh@gfy.ku.dk/thomas@cultpenguin.com
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
function PutSegyTrace(segyid,tracedata,SegyTraceHeader,SegyHeader);

% *** DUE TO A BUG IN MATLAB 6.5 : Technical Solution Number: 31977
  
% Enable next line if error messages show that not all
% trace header values have been set.
% Enabling this check will considerably slow down writing 
% with  factor of about 6!
%  SegyTraceHeader=CheckSegyTraceHeader(SegyTraceHeader);
  

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WRITE SegyTraceHeader
    
  if ~exist('TraceStart')==1, TraceStart=ftell(segyid);end
  
  fseek(segyid,0,'bof');fseek(segyid,TraceStart,'bof'); % ***
  fwrite(segyid,ones(1,60),'int32');
  fseek(segyid,0,'bof');fseek(segyid,TraceStart,'bof'); % ***

  WRITETRACEHEADER=1;
  if WRITETRACEHEADER==1,
    fseek(segyid,0,'bof');
    fseek(segyid,TraceStart,'bof');
    fwrite(segyid,SegyTraceHeader.TraceSequenceLine,'int32');    % 0
    fwrite(segyid,SegyTraceHeader.TraceSequenceFile,'int32');    % 4 
    fwrite(segyid,SegyTraceHeader.FieldRecord,'int32');          % 8
    fwrite(segyid,SegyTraceHeader.TraceNumber,'int32');          % 12
    fwrite(segyid,SegyTraceHeader.EnergySourcePoint,'int32');    % 16
    fwrite(segyid,SegyTraceHeader.cdp,'int32');                  % 20
    fwrite(segyid,SegyTraceHeader.cdpTrace,'int32');             % 24
    fwrite(segyid,SegyTraceHeader.TraceIdenitifactionCode,'int16'); % 28
    fwrite(segyid,SegyTraceHeader.NSummedTraces,'int16'); % 30
    fwrite(segyid,SegyTraceHeader.NStackedTraces,'int16'); % 32
    fwrite(segyid,SegyTraceHeader.DataUse,'int16'); % 34
    fwrite(segyid,SegyTraceHeader.offset,'int32');             %36
    fwrite(segyid,SegyTraceHeader.ReceiverGroupElevation,'int32');             %40
    fwrite(segyid,SegyTraceHeader.SourceSurfaceElevation,'int32');             %44
    fwrite(segyid,SegyTraceHeader.SourceDepth,'int32');             %48
    fwrite(segyid,SegyTraceHeader.ReceiverDatumElevation,'int32');             %52
    fwrite(segyid,SegyTraceHeader.SourceDatumElevation,'int32');             %56
    fwrite(segyid,SegyTraceHeader.SourceWaterDepth,'int32');  %60
    fwrite(segyid,SegyTraceHeader.GroupWaterDepth,'int32');  %64
    fwrite(segyid,SegyTraceHeader.ElevationScalar,'int16');  %68
                                                             % Multiply/divide next number for following 4 values
    fwrite(segyid,SegyTraceHeader.SourceGroupScalar,'int16');  %70
    fwrite(segyid,SegyTraceHeader.SourceX,'int32');  %72
    fwrite(segyid,SegyTraceHeader.SourceY,'int32');  %76
    fwrite(segyid,SegyTraceHeader.GroupX,'int32');  %80
    fwrite(segyid,SegyTraceHeader.GroupY,'int32');  %84
    fwrite(segyid,SegyTraceHeader.CoordinateUnits,'int16');  %88
    fwrite(segyid,SegyTraceHeader.WeatheringVelocity,'int16');  %90
    fwrite(segyid,SegyTraceHeader.SubWeatheringVelocity,'int16');  %92
    fwrite(segyid,SegyTraceHeader.SourceUpholeTime,'int16');  %94
    fwrite(segyid,SegyTraceHeader.GroupUpholeTime,'int16');  %96
    fwrite(segyid,SegyTraceHeader.SourceStaticCorrection,'int16');  %98
    fwrite(segyid,SegyTraceHeader.GroupStaticCorrection,'int16');  %100
    fwrite(segyid,SegyTraceHeader.TotalStaticApplied,'int16');  %102
    fwrite(segyid,SegyTraceHeader.LagTimeA,'int16');  %104
    fwrite(segyid,SegyTraceHeader.LagTimeB,'int16');  %106
    fwrite(segyid,SegyTraceHeader.DelayRecordingTime,'int16');  %108
    fwrite(segyid,SegyTraceHeader.MuteTimeStart,'int16');  %110
    fwrite(segyid,SegyTraceHeader.MuteTimeEND,'int16');  %112
    fwrite(segyid,SegyTraceHeader.ns,'uint16');  %114
    fwrite(segyid,SegyTraceHeader.dt,'uint16');  %116
    fwrite(segyid,SegyTraceHeader.GainType,'int16');  %118
    fwrite(segyid,SegyTraceHeader.InstrumentGainConstant,'int16');  %120
    fwrite(segyid,SegyTraceHeader.InstrumentInitialGain,'int16');  %%122
    fwrite(segyid,SegyTraceHeader.Correlated,'int16');  %124
    fwrite(segyid,SegyTraceHeader.SweepFrequenceStart,'int16');  %126
    fwrite(segyid,SegyTraceHeader.SweepFrequenceEnd,'int16');  %128
    fwrite(segyid,SegyTraceHeader.SweepLength,'int16');  %130
    fwrite(segyid,SegyTraceHeader.SweepType,'int16');  %132
    fwrite(segyid,SegyTraceHeader.SweepTraceTaperLengthStart,'int16');  %134
    fwrite(segyid,SegyTraceHeader.SweepTraceTaperLengthEnd,'int16');  %136
    fwrite(segyid,SegyTraceHeader.TaperType,'int16');  %138
    fwrite(segyid,SegyTraceHeader.AliasFilterFrequency,'int16');  %140
    fwrite(segyid,SegyTraceHeader.AliasFilterSlope,'int16');  %142
    fwrite(segyid,SegyTraceHeader.NotchFilterFrequency,'int16');  %144
    fwrite(segyid,SegyTraceHeader.NotchFilterSlope,'int16');  %146
    fwrite(segyid,SegyTraceHeader.LowCutFrequency,'int16');  %148
    fwrite(segyid,SegyTraceHeader.HighCutFrequency,'int16');  %150
    fwrite(segyid,SegyTraceHeader.LowCutSlope,'int16');  %152
    fwrite(segyid,SegyTraceHeader.HighCutSlope,'int16');  %154
    fwrite(segyid,SegyTraceHeader.YearDataRecorded,'int16');  %156
    fwrite(segyid,SegyTraceHeader.DayOfYear,'int16');  %158
    fwrite(segyid,SegyTraceHeader.HourOfDay,'int16');  %160
    fwrite(segyid,SegyTraceHeader.MinuteOfHour,'int16');  %162
    fwrite(segyid,SegyTraceHeader.SecondOfMinute,'int16');  %164
    fwrite(segyid,SegyTraceHeader.TimeBaseCode,'int16');  %166
    fwrite(segyid,SegyTraceHeader.TraceWeightningFactor,'int16');  %170
    fwrite(segyid,SegyTraceHeader.GeophoneGroupNumberRoll1,'int16');  %172
    fwrite(segyid,SegyTraceHeader.GeophoneGroupNumberFirstTraceOrigField,'int16');  %174
    fwrite(segyid,SegyTraceHeader.GeophoneGroupNumberLastTraceOrigField,'int16');  %176
    fwrite(segyid,SegyTraceHeader.GapSize,'int16');  %178
    fwrite(segyid,SegyTraceHeader.OverTravel,'int16');  %178
    fwrite(segyid,SegyTraceHeader.cdpX,'int32');  %180
    fwrite(segyid,SegyTraceHeader.cdpY,'int32');  %184
    fwrite(segyid,SegyTraceHeader.Inline3D,'int32');  %188
    fwrite(segyid,SegyTraceHeader.Crossline3D,'int32');  %192
    fwrite(segyid,SegyTraceHeader.ShotPoint,'int32');  %196
    fwrite(segyid,SegyTraceHeader.ShotPointScalar,'int16');  %200
    fwrite(segyid,SegyTraceHeader.TraceValueMeasurementUnit,'int16');  %202
    fwrite(segyid,SegyTraceHeader.TransductionConstantMantissa,'int32');  %204
    fwrite(segyid,SegyTraceHeader.TransductionConstantPower,'int16'); %208
    fwrite(segyid,SegyTraceHeader.TransductionUnit,'int16');  %210
    fwrite(segyid,SegyTraceHeader.TraceIdentifier,'int16');  %212
    fwrite(segyid,SegyTraceHeader.ScalarTraceHeader,'int16');  %214
    fwrite(segyid,SegyTraceHeader.SourceType,'int16');  %216
    fwrite(segyid,SegyTraceHeader.SourceEnergyDirectionMantissa,'int32');  %218
    fwrite(segyid,SegyTraceHeader.SourceEnergyDirectionExponent,'int16');  %222
    fwrite(segyid,SegyTraceHeader.SourceMeasurementMantissa,'int32');  %224
    fwrite(segyid,SegyTraceHeader.SourceMeasurementExponent,'int16');  %228
    fwrite(segyid,SegyTraceHeader.SourceMeasurementUnit,'int16');  %230
                                                                   % WRITE UNASSIGNED CHARACTERS FOR THE REST
    fwrite(segyid,SegyTraceHeader.UnassignedInt1,'int32');  %232
    fwrite(segyid,SegyTraceHeader.UnassignedInt2,'int32');  %236
                                                            %fwrite(segyid,zeros(1,4),'int16');  %230
  end
  % 217-240 Unassigned
  
  % go to end of header 
  % Any of the nex to lines should work, with the first line being the fastest  
  fseek(segyid,0,'cof'); fseek(segyid,240-216,'cof');
  % fseek(segyid,0,'bof'); fseek(segyid,TraceStart+240,'bof');
  
  
  Revision=SegyHeader.SegyFormatRevisionNumber;
  if Revision>0, Revision=1; end
  Format=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).format;  
  BPS=SegyHeader.Rev(Revision+1).DataSampleFormat(SegyHeader.DataSampleFormat).bps; 
  
  SegymatVerbose([mfilename,' SegyRevision ',sprintf('%0.4g',Revision),', ',Format],3);
  
  %% WRITE TRACE DATA;
  if (strcmp(Format,'uint32')==1)|(strcmp(Format,'uint16')==1), % IBM FLOATING POINT
                                                                % CONVERT FROM FLOATING POINT
    SegymatVerbose([mfilename,'Converting from IBM, DataFormat :',SegyHeader.DataSampleFormat],2); 
    tracedata=double(num2ibm(tracedata));
  end;
  if SegyHeader.FixedLengthTraceFlag
      fwrite(segyid,tracedata,Format);
  else
      fwrite(segyid,tracedata(1:SegyTraceHeader.ns),Format);
  end
  