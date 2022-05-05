% SegyTraceHeader=CheckSegyTraceHeader(SegyTraceHeader);
%
% Checks that all fields of the SegyTraceHeader is set. 
% If not, they are initialized.
%

%
% (C) 2004, Thomas Mejer Hansen, tmh@gfy.ku.dk/thomas@cultpenguin.com
%
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
function SegyTraceHeader=CheckSegyTraceHeader(SegyTraceHeader);
%% SET UP SegyTraceHeader
  NA=0;
  
  if ~isfield(SegyTraceHeader,'TraceSequenceLine'), SegyTraceHeader.TraceSequenceLine=NA;end
  if ~isfield(SegyTraceHeader,'TraceSequenceFile'), SegyTraceHeader.TraceSequenceFile=NA;end
  if ~isfield(SegyTraceHeader,'FieldRecord'), SegyTraceHeader.FieldRecord=NA;end
  if ~isfield(SegyTraceHeader,'TraceNumber'), SegyTraceHeader.TraceNumber=NA;end
  if ~isfield(SegyTraceHeader,'EnergySourcePoint'), SegyTraceHeader.EnergySourcePoint=NA;end
  if ~isfield(SegyTraceHeader,'cdp'), SegyTraceHeader.cdp=NA;end
  if ~isfield(SegyTraceHeader,'cdpTrace'), SegyTraceHeader.cdpTrace=NA;end
  if ~isfield(SegyTraceHeader,'TraceIdenitifactionCode'), SegyTraceHeader.TraceIdenitifactionCode=NA;end
  if ~isfield(SegyTraceHeader,'NSummedTraces'), SegyTraceHeader.NSummedTraces=NA;end
  if ~isfield(SegyTraceHeader,'NStackedTraces'), SegyTraceHeader.NStackedTraces=NA;end
  if ~isfield(SegyTraceHeader,'DataUse'), SegyTraceHeader.DataUse=NA;end
  if ~isfield(SegyTraceHeader,'offset'), SegyTraceHeader.offset=NA;end
  if ~isfield(SegyTraceHeader,'ReceiverGroupElevation'), SegyTraceHeader.ReceiverGroupElevation=NA;end
  if ~isfield(SegyTraceHeader,'SourceSurfaceElevation'), SegyTraceHeader.SourceSurfaceElevation=NA;end
  if ~isfield(SegyTraceHeader,'SourceDepth'), SegyTraceHeader.SourceDepth=NA;end
  if ~isfield(SegyTraceHeader,'ReceiverDatumElevation'), SegyTraceHeader.ReceiverDatumElevation=NA;end
  if ~isfield(SegyTraceHeader,'SourceDatumElevation'), SegyTraceHeader.SourceDatumElevation=NA;end
  if ~isfield(SegyTraceHeader,'SourceWaterDepth'), SegyTraceHeader.SourceWaterDepth=NA;end
  if ~isfield(SegyTraceHeader,'GroupWaterDepth'), SegyTraceHeader.GroupWaterDepth=NA;end
  if ~isfield(SegyTraceHeader,'ElevationScalar'), SegyTraceHeader.ElevationScalar=NA;end
  % Multiply/divide next number for following 4 values
  if ~isfield(SegyTraceHeader,'SourceGroupScalar'), SegyTraceHeader.SourceGroupScalar=NA;end
  if ~isfield(SegyTraceHeader,'SourceX'), SegyTraceHeader.SourceX=NA;end
  if ~isfield(SegyTraceHeader,'SourceY'), SegyTraceHeader.SourceY=NA;end
  if ~isfield(SegyTraceHeader,'GroupX'), SegyTraceHeader.GroupX=NA;end
  if ~isfield(SegyTraceHeader,'GroupY'), SegyTraceHeader.GroupY=NA;end
  if ~isfield(SegyTraceHeader,'CoordinateUnits'), SegyTraceHeader.CoordinateUnits=NA;end
  if ~isfield(SegyTraceHeader,'WeatheringVelocity'), SegyTraceHeader.WeatheringVelocity=NA;end
  if ~isfield(SegyTraceHeader,'SubWeatheringVelocity'), SegyTraceHeader.SubWeatheringVelocity=NA;end
  if ~isfield(SegyTraceHeader,'SourceUpholeTime'), SegyTraceHeader.SourceUpholeTime=NA;end
  if ~isfield(SegyTraceHeader,'GroupUpholeTime'), SegyTraceHeader.GroupUpholeTime=NA;end
  if ~isfield(SegyTraceHeader,'SourceStaticCorrection'), SegyTraceHeader.SourceStaticCorrection=NA;end
  if ~isfield(SegyTraceHeader,'GroupStaticCorrection'), SegyTraceHeader.GroupStaticCorrection=NA;end
  if ~isfield(SegyTraceHeader,'TotalStaticApplied'), SegyTraceHeader.TotalStaticApplied=NA;end
  if ~isfield(SegyTraceHeader,'LagTimeA'), SegyTraceHeader.LagTimeA=NA;end
  if ~isfield(SegyTraceHeader,'LagTimeB'), SegyTraceHeader.LagTimeB=NA;end
  if ~isfield(SegyTraceHeader,'DelayRecordingTime'), SegyTraceHeader.DelayRecordingTime=NA;end
  if ~isfield(SegyTraceHeader,'MuteTimeStart'), SegyTraceHeader.MuteTimeStart=NA;end
  if ~isfield(SegyTraceHeader,'MuteTimeEND'), SegyTraceHeader.MuteTimeEND=NA;end
  if ~isfield(SegyTraceHeader,'ns'), SegyTraceHeader.ns=NA;end
  if ~isfield(SegyTraceHeader,'dt'), SegyTraceHeader.dt=NA;end
  if ~isfield(SegyTraceHeader,'GainType'), SegyTraceHeader.GainType=NA;end
  if ~isfield(SegyTraceHeader,'InstrumentGainConstant'), SegyTraceHeader.InstrumentGainConstant=NA;end
  if ~isfield(SegyTraceHeader,'InstrumentInitialGain'), SegyTraceHeader.InstrumentInitialGain=NA;end
  if ~isfield(SegyTraceHeader,'Correlated'), SegyTraceHeader.Correlated=NA;end
  if ~isfield(SegyTraceHeader,'SweepFrequenceStart'), SegyTraceHeader.SweepFrequenceStart=NA;end
  if ~isfield(SegyTraceHeader,'SweepFrequenceEnd'), SegyTraceHeader.SweepFrequenceEnd=NA;end
  if ~isfield(SegyTraceHeader,'SweepLength'), SegyTraceHeader.SweepLength=NA;end
  if ~isfield(SegyTraceHeader,'SweepType'), SegyTraceHeader.SweepType=NA;end
  if ~isfield(SegyTraceHeader,'SweepTraceTaperLengthStart'), SegyTraceHeader.SweepTraceTaperLengthStart=NA;end
  if ~isfield(SegyTraceHeader,'SweepTraceTaperLengthEnd'), SegyTraceHeader.SweepTraceTaperLengthEnd=NA;end
  if ~isfield(SegyTraceHeader,'TaperType'), SegyTraceHeader.TaperType=NA;end
  if ~isfield(SegyTraceHeader,'AliasFilterFrequency'), SegyTraceHeader.AliasFilterFrequency=NA;end
  if ~isfield(SegyTraceHeader,'AliasFilterSlope'), SegyTraceHeader.AliasFilterSlope=NA;end
  if ~isfield(SegyTraceHeader,'NotchFilterFrequency'), SegyTraceHeader.NotchFilterFrequency=NA;end
  if ~isfield(SegyTraceHeader,'NotchFilterSlope'), SegyTraceHeader.NotchFilterSlope=NA;end
  if ~isfield(SegyTraceHeader,'LowCutFrequency'), SegyTraceHeader.LowCutFrequency=NA;end
  if ~isfield(SegyTraceHeader,'HighCutFrequency'), SegyTraceHeader.HighCutFrequency=NA;end
  if ~isfield(SegyTraceHeader,'LowCutSlope'), SegyTraceHeader.LowCutSlope=NA;end
  if ~isfield(SegyTraceHeader,'HighCutSlope'), SegyTraceHeader.HighCutSlope=NA;end
  da=clock;
  if ~isfield(SegyTraceHeader,'YearDataRecorded'), SegyTraceHeader.YearDataRecorded=da(1);end
  DayOfYear=datenum(0,da(2),da(3)); % 22/7-2002 by Michael Toews (mwtoews@ucalgary.ca)
  if ~isfield(SegyTraceHeader,'DayOfYear'), SegyTraceHeader.DayOfYear=DayOfYear;end
  if ~isfield(SegyTraceHeader,'HourOfDay'), SegyTraceHeader.HourOfDay=da(4);end
  if ~isfield(SegyTraceHeader,'MinuteOfHour'), SegyTraceHeader.MinuteOfHour=da(5);end
  if ~isfield(SegyTraceHeader,'SecondOfMinute'), SegyTraceHeader.SecondOfMinute=round(da(6));end
  if ~isfield(SegyTraceHeader,'TimeBaseCode'), SegyTraceHeader.TimeBaseCode=NA;end
  if ~isfield(SegyTraceHeader,'TraceWeightningFactor'), SegyTraceHeader.TraceWeightningFactor=NA;end
  if ~isfield(SegyTraceHeader,'GeophoneGroupNumberRoll1'), SegyTraceHeader.GeophoneGroupNumberRoll1=NA;end
  if ~isfield(SegyTraceHeader,'GeophoneGroupNumberFirstTraceOrigField'), SegyTraceHeader.GeophoneGroupNumberFirstTraceOrigField=NA;end
  if ~isfield(SegyTraceHeader,'GeophoneGroupNumberLastTraceOrigField'), SegyTraceHeader.GeophoneGroupNumberLastTraceOrigField=NA;end
  if ~isfield(SegyTraceHeader,'GapSize'), SegyTraceHeader.GapSize=NA;end
  if ~isfield(SegyTraceHeader,'OverTravel'), SegyTraceHeader.OverTravel=NA;end
  if ~isfield(SegyTraceHeader,'cdpX'), SegyTraceHeader.cdpX=NA;end
  if ~isfield(SegyTraceHeader,'cdpY'), SegyTraceHeader.cdpY=NA;end
  if ~isfield(SegyTraceHeader,'Inline3D'), SegyTraceHeader.Inline3D=NA;end
  if ~isfield(SegyTraceHeader,'Crossline3D'), SegyTraceHeader.Crossline3D=NA;end
  if ~isfield(SegyTraceHeader,'ShotPoint'), SegyTraceHeader.ShotPoint=NA;end
  if ~isfield(SegyTraceHeader,'ShotPointScalar'), SegyTraceHeader.ShotPointScalar=NA;end
  if ~isfield(SegyTraceHeader,'TraceValueMeasurementUnit'), SegyTraceHeader.TraceValueMeasurementUnit=NA;end
  if ~isfield(SegyTraceHeader,'TransductionConstantMantissa'), SegyTraceHeader.TransductionConstantMantissa=NA;end
  if ~isfield(SegyTraceHeader,'TransductionConstantPower'), SegyTraceHeader.TransductionConstantPower=NA;end
  if ~isfield(SegyTraceHeader,'TransductionUnit'), SegyTraceHeader.TransductionUnit=NA;end
  if ~isfield(SegyTraceHeader,'TraceIdentifier'), SegyTraceHeader.TraceIdentifier=NA;end
  if ~isfield(SegyTraceHeader,'ScalarTraceHeader'), SegyTraceHeader.ScalarTraceHeader=NA;end
  
  
  if ~isfield(SegyTraceHeader,'SourceType'), SegyTraceHeader.SourceType=NA;end
  if ~isfield(SegyTraceHeader,'SourceEnergyDirectionMantissa'), SegyTraceHeader.SourceEnergyDirectionMantissa=NA;end
  if ~isfield(SegyTraceHeader,'SourceEnergyDirectionExponent'), SegyTraceHeader.SourceEnergyDirectionExponent=NA;end
  if ~isfield(SegyTraceHeader,'SourceMeasurementMantissa'), SegyTraceHeader.SourceMeasurementMantissa=NA;end
  if ~isfield(SegyTraceHeader,'SourceMeasurementExponent'), SegyTraceHeader.SourceMeasurementExponent=NA;end
  if ~isfield(SegyTraceHeader,'SourceMeasurementUnit'), SegyTraceHeader.SourceMeasurementUnit=NA;end
  
  if ~isfield(SegyTraceHeader,'UnassignedInt1'), SegyTraceHeader.UnassignedInt1=NA;end
  if ~isfield(SegyTraceHeader,'UnassignedInt2'), SegyTraceHeader.UnassignedInt2=NA;end
  
