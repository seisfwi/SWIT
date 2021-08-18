% GetSegyTraceHeader : Reads a seg y trace, data and header
%
% [SegyTraceHeader]=GetSegyTraceHeader(segyid,TraceStart,DataFormat,ns);
%
%
% (C) 2001-2004 Thomas Mejer Hansen, thomas.mejer.hansen@cultpenguin.com
% 
% Revisions:
% 07/2008 Kristian Stormark (<kristian.stormark@gmail.com>) : Reduce the 
%         number of discoperations causing a significant speed up 
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

function [SegyTraceHeader]=GetSegyTraceHeader(segyid,TraceStart,DataFormat,ns,SegyTraceHeader);

if exist('DataFormat')==0, DataFormat='float32'; end
if exist('TraceStart')==0, TraceStart=ftell(segyid); end

if exist('SegyTraceHeader')
    if isempty('SegyTraceHeader');
        clear SegyTraceHeader;
    end
end


fseek(segyid,TraceStart,'bof');

% GET POSITION FOR EASY LATER LOCALIZATION
SegyTraceHeader.SegyMAT_TraceStart = ftell(segyid);

chunk = fread(segyid,7,'int32');
SegyTraceHeader.TraceSequenceLine = chunk(1);    % 0
SegyTraceHeader.TraceSequenceFile = chunk(2);    % 4 
SegyTraceHeader.FieldRecord       = chunk(3);    % 8
SegyTraceHeader.TraceNumber       = chunk(4);    % 12
SegyTraceHeader.EnergySourcePoint = chunk(5);    % 16
SegyTraceHeader.cdp               = chunk(6);    % 20
SegyTraceHeader.cdpTrace          = chunk(7);    % 24

chunk = fread(segyid,4,'int16');
SegyTraceHeader.TraceIdenitifactionCode = chunk(1); % 28
SegyTraceHeader.NSummedTraces           = chunk(2); % 30
SegyTraceHeader.NStackedTraces          = chunk(3); % 32
SegyTraceHeader.DataUse                 = chunk(4); % 34

chunk = fread(segyid,8,'int32');
SegyTraceHeader.offset                  = chunk(1);  %36
SegyTraceHeader.ReceiverGroupElevation  = chunk(2);  %40
SegyTraceHeader.SourceSurfaceElevation  = chunk(3);  %44
SegyTraceHeader.SourceDepth             = chunk(4);  %48
SegyTraceHeader.ReceiverDatumElevation  = chunk(5);  %52
SegyTraceHeader.SourceDatumElevation    = chunk(6);  %56
SegyTraceHeader.SourceWaterDepth        = chunk(7);  %60
SegyTraceHeader.GroupWaterDepth         = chunk(8);  %64


chunk = fread(segyid,2,'int16');
SegyTraceHeader.ElevationScalar   = chunk(1);  %68
% Multiply/divide next number for following 4 values
SegyTraceHeader.SourceGroupScalar = chunk(2);  %70

chunk = fread(segyid,4,'int32');
SegyTraceHeader.SourceX = chunk(1); %72
SegyTraceHeader.SourceY = chunk(2); %76
SegyTraceHeader.GroupX  = chunk(3); %80
SegyTraceHeader.GroupY  = chunk(4); %84

chunk = fread(segyid,13,'int16');
SegyTraceHeader.CoordinateUnits       	= chunk(1);  %88
SegyTraceHeader.WeatheringVelocity    	= chunk(2);  %90
SegyTraceHeader.SubWeatheringVelocity 	= chunk(3);  %92
SegyTraceHeader.SourceUpholeTime		= chunk(4);  %94
SegyTraceHeader.GroupUpholeTime			= chunk(5);  %96
SegyTraceHeader.SourceStaticCorrection	= chunk(6);  %98
SegyTraceHeader.GroupStaticCorrection	= chunk(7);  %100
SegyTraceHeader.TotalStaticApplied		= chunk(8);  %102
SegyTraceHeader.LagTimeA				= chunk(9);  %104
SegyTraceHeader.LagTimeB				= chunk(10);  %106
SegyTraceHeader.DelayRecordingTime		= chunk(11);  %108
SegyTraceHeader.MuteTimeStart			= chunk(12);  %110
SegyTraceHeader.MuteTimeEND				= chunk(13);  %112

chunk = fread(segyid,2,'uint16');
SegyTraceHeader.ns = chunk(1);  %114
SegyTraceHeader.dt = chunk(2);  %116

chunk = fread(segyid,31,'int16');
SegyTraceHeader.GainType					= chunk(1);  %118
SegyTraceHeader.InstrumentGainConstant		= chunk(2);  %120
SegyTraceHeader.InstrumentInitialGain		= chunk(3);  %%122
SegyTraceHeader.Correlated					= chunk(4);  %124
SegyTraceHeader.SweepFrequenceStart			= chunk(5);  %126
SegyTraceHeader.SweepFrequenceEnd			= chunk(6);  %128
SegyTraceHeader.SweepLength					= chunk(7);  %130
SegyTraceHeader.SweepType					= chunk(8);  %132
SegyTraceHeader.SweepTraceTaperLengthStart	= chunk(9);  %134
SegyTraceHeader.SweepTraceTaperLengthEnd	= chunk(10);  %136
SegyTraceHeader.TaperType					= chunk(11);  %138
SegyTraceHeader.AliasFilterFrequency		= chunk(12);  %140
SegyTraceHeader.AliasFilterSlope			= chunk(13);  %142
SegyTraceHeader.NotchFilterFrequency		= chunk(14);  %144
SegyTraceHeader.NotchFilterSlope			= chunk(15);  %146
SegyTraceHeader.LowCutFrequency				= chunk(16);  %148
SegyTraceHeader.HighCutFrequency			= chunk(17);  %150
SegyTraceHeader.LowCutSlope					= chunk(18);  %152
SegyTraceHeader.HighCutSlope				= chunk(19);  %154
SegyTraceHeader.YearDataRecorded			= chunk(20);  %156
SegyTraceHeader.DayOfYear					= chunk(21);  %158
SegyTraceHeader.HourOfDay					= chunk(22);  %160
SegyTraceHeader.MinuteOfHour				= chunk(23);  %162
SegyTraceHeader.SecondOfMinute				= chunk(24);  %164
SegyTraceHeader.TimeBaseCode				= chunk(25);  %166


if SegyTraceHeader.TimeBaseCode==1, SegyTraceHeader.TimeBaseCodeText='Local'; 
elseif SegyTraceHeader.TimeBaseCode==2, SegyTraceHeader.TimeBaseCodeText='GMT';
elseif SegyTraceHeader.TimeBaseCode==3, SegyTraceHeader.TimeBaseCodeText='Other';
elseif SegyTraceHeader.TimeBaseCode==4, SegyTraceHeader.TimeBaseCodeText='UTC';
else SegyTraceHeader.TimeBaseCodeText=''; end


SegyTraceHeader.TraceWeightningFactor					= chunk(26);  %168
SegyTraceHeader.GeophoneGroupNumberRoll1				= chunk(27);  %170
SegyTraceHeader.GeophoneGroupNumberFirstTraceOrigField	= chunk(28);  %172
SegyTraceHeader.GeophoneGroupNumberLastTraceOrigField	= chunk(29);  %174
SegyTraceHeader.GapSize									= chunk(30);  %176
SegyTraceHeader.OverTravel								= chunk(31);  %178

chunk = fread(segyid,5,'int32');
SegyTraceHeader.cdpX		= chunk(1); %180
SegyTraceHeader.cdpY		= chunk(2); %184
SegyTraceHeader.Inline3D	= chunk(3); %188
SegyTraceHeader.Crossline3D	= chunk(4); %192
SegyTraceHeader.ShotPoint	= chunk(5); %196

SegyTraceHeader.ShotPointScalar=fread(segyid,1,'int16');  %200

SegyTraceHeader.TraceValueMeasurementUnit=fread(segyid,1,'int16');  %202
if SegyTraceHeader.TraceValueMeasurementUnit==-1, SegyTraceHeader.TraceValueMeasurementUnitText='Other';
elseif SegyTraceHeader.TraceValueMeasurementUnit==0, SegyTraceHeader.TraceValueMeasurementUnitText='Unknown';
elseif SegyTraceHeader.TraceValueMeasurementUnit==1, SegyTraceHeader.TraceValueMeasurementUnitText='Pascal (Pa)';
elseif SegyTraceHeader.TraceValueMeasurementUnit==2, SegyTraceHeader.TraceValueMeasurementUnitText='Volts (v)';
elseif SegyTraceHeader.TraceValueMeasurementUnit==3, SegyTraceHeader.TraceValueMeasurementUnitText='Millivolts (v)';
elseif SegyTraceHeader.TraceValueMeasurementUnit==4, SegyTraceHeader.TraceValueMeasurementUnitText='Amperes (A)';  
elseif SegyTraceHeader.TraceValueMeasurementUnit==5, SegyTraceHeader.TraceValueMeasurementUnitText='Meters (m)';  
elseif SegyTraceHeader.TraceValueMeasurementUnit==6, SegyTraceHeader.TraceValueMeasurementUnitText='Meters Per Second (m/s)';  
elseif SegyTraceHeader.TraceValueMeasurementUnit==7, SegyTraceHeader.TraceValueMeasurementUnitText='Meters Per Second squared (m/&s2)Other';  
elseif SegyTraceHeader.TraceValueMeasurementUnit==8, SegyTraceHeader.TraceValueMeasurementUnitText='Newton (N)';  
elseif SegyTraceHeader.TraceValueMeasurementUnit==8, SegyTraceHeader.TraceValueMeasurementUnitText='Watt (W)';  
else SegyTraceHeader.TraceValueMeasurementUnitText='Undefined'; end

SegyTraceHeader.TransductionConstantMantissa=fread(segyid,1,'int32');  %204


chunk = fread(segyid,5,'int16');
SegyTraceHeader.TransductionConstantPower	= chunk(1);%208
SegyTraceHeader.TransductionUnit			= chunk(2); %210
SegyTraceHeader.TraceIdentifier				= chunk(3); %212
SegyTraceHeader.ScalarTraceHeader			= chunk(4); %214
SegyTraceHeader.SourceType					= chunk(5); %216


SegyTraceHeader.SourceEnergyDirectionMantissa=fread(segyid,1,'int32');  %218
SegyTraceHeader.SourceEnergyDirectionExponent=fread(segyid,1,'int16');  %222
SegyTraceHeader.SourceMeasurementMantissa=fread(segyid,1,'int32');  %224

chunk = fread(segyid,2,'int16');
SegyTraceHeader.SourceMeasurementExponent	= chunk(1);  %228
SegyTraceHeader.SourceMeasurementUnit		= chunk(2);  %230

chunk = fread(segyid,2,'int32');
SegyTraceHeader.UnassignedInt1 = chunk(1);  %232
SegyTraceHeader.UnassignedInt2 = chunk(2);  %236


if exist('ns')==0,
  ns=SegyTraceHeader.ns;
end


% GO TO POSITION OF DATA
fseek(segyid,TraceStart+240,'bof');
SegyTraceHeader.SegyMAT_TraceDataStart = ftell(segyid);

