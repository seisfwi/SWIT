% InitSegyTraceHeaders : returns an empty SegyTraceHeader structure
%
% EX:
% SegyTraceHeader=InitSegyTraceHeader(ns,dt);
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
function SegyTraceHeader=InitSegyTraceHeader(ns,dt);
  
  SegyTraceHeader=[];
  
  if exist('ns')==0, 
    SegymatVerbose([mfilename,' : ns not set !! returning'],0)
    return
  end
  
  if exist('dt')==0, 
    dt=4000;
    SegymatVerbose([mfilename,' : dt not set. using dt=',num2str(dt)],0)    
  end

  NA=0;
  da=clock;
  DayOfYear=datenum(0,da(2),da(3));

  
  SegyTraceHeader.TraceSequenceLine=NA;
  SegyTraceHeader.TraceSequenceFile=NA;
  SegyTraceHeader.FieldRecord=NA;
  SegyTraceHeader.TraceNumber=NA;;
  SegyTraceHeader.EnergySourcePoint=NA;
  SegyTraceHeader.cdp=NA;
  SegyTraceHeader.cdpTrace=NA;
  SegyTraceHeader.TraceIdenitifactionCode=NA;
  SegyTraceHeader.NSummedTraces=NA;
  SegyTraceHeader.NStackedTraces=NA;
  SegyTraceHeader.DataUse=NA;
  SegyTraceHeader.offset=NA;
  SegyTraceHeader.ReceiverGroupElevation=NA;
  SegyTraceHeader.SourceSurfaceElevation=NA;
  SegyTraceHeader.SourceDepth=NA;
  SegyTraceHeader.ReceiverDatumElevation=NA;
  SegyTraceHeader.SourceDatumElevation=NA;
  SegyTraceHeader.SourceWaterDepth=NA;
  SegyTraceHeader.GroupWaterDepth=NA;
  SegyTraceHeader.ElevationScalar=NA;
  SegyTraceHeader.SourceGroupScalar=NA;
  SegyTraceHeader.SourceX=NA;
  SegyTraceHeader.SourceY=NA;
  SegyTraceHeader.GroupX=NA;
  SegyTraceHeader.GroupY=NA;
  SegyTraceHeader.CoordinateUnits=NA;
  SegyTraceHeader.WeatheringVelocity=NA;
  SegyTraceHeader.SubWeatheringVelocity=NA;
  SegyTraceHeader.SourceUpholeTime=NA;
  SegyTraceHeader.GroupUpholeTime=NA;
  SegyTraceHeader.SourceStaticCorrection=NA;
  SegyTraceHeader.GroupStaticCorrection=NA;
  SegyTraceHeader.TotalStaticApplied=NA;
  SegyTraceHeader.LagTimeA=NA;
  SegyTraceHeader.LagTimeB=NA;
  SegyTraceHeader.DelayRecordingTime=NA;
  SegyTraceHeader.MuteTimeStart=NA;
  SegyTraceHeader.MuteTimeEND=NA;
  SegyTraceHeader.ns=ns;
  SegyTraceHeader.dt=dt;
  SegyTraceHeader.GainType=NA;
  SegyTraceHeader.InstrumentGainConstant=NA;
  SegyTraceHeader.InstrumentInitialGain=NA;
  SegyTraceHeader.Correlated=NA;
  SegyTraceHeader.SweepFrequenceStart=NA;
  SegyTraceHeader.SweepFrequenceEnd=NA;
  SegyTraceHeader.SweepLength=NA;
  SegyTraceHeader.SweepType=NA;
  SegyTraceHeader.SweepTraceTaperLengthStart=NA;
  SegyTraceHeader.SweepTraceTaperLengthEnd=NA;
  SegyTraceHeader.TaperType=NA;
  SegyTraceHeader.AliasFilterFrequency=NA;
  SegyTraceHeader.AliasFilterSlope=NA;
  SegyTraceHeader.NotchFilterFrequency=NA;
  SegyTraceHeader.NotchFilterSlope=NA;
  SegyTraceHeader.LowCutFrequency=NA;
  SegyTraceHeader.HighCutFrequency=NA;
  SegyTraceHeader.LowCutSlope=NA;
  SegyTraceHeader.HighCutSlope=NA;
  
  SegyTraceHeader.YearDataRecorded=da(1);
  SegyTraceHeader.DayOfYear=DayOfYear;
  SegyTraceHeader.HourOfDay=da(4);
  SegyTraceHeader.MinuteOfHour=da(5);
  SegyTraceHeader.SecondOfMinute=round(da(6));
  SegyTraceHeader.TimeBaseCode=NA;
  SegyTraceHeader.TraceWeightningFactor=NA;
  SegyTraceHeader.GeophoneGroupNumberRoll1=NA;
  SegyTraceHeader.GeophoneGroupNumberFirstTraceOrigField=NA;
  SegyTraceHeader.GeophoneGroupNumberLastTraceOrigField=NA;
  SegyTraceHeader.GapSize=NA;
  SegyTraceHeader.OverTravel=NA;
  SegyTraceHeader.cdpX=NA;
  SegyTraceHeader.cdpY=NA;
  SegyTraceHeader.Inline3D=NA;
  SegyTraceHeader.Crossline3D=NA;
  SegyTraceHeader.ShotPoint=NA;
  SegyTraceHeader.ShotPointScalar=NA;
  SegyTraceHeader.TraceValueMeasurementUnit=NA;
  SegyTraceHeader.TransductionConstantMantissa=NA;
  SegyTraceHeader.TransductionConstantPower=NA;
  SegyTraceHeader.TransductionUnit=NA;
  SegyTraceHeader.TraceIdentifier=NA;
  SegyTraceHeader.ScalarTraceHeader=NA;
  
  
  SegyTraceHeader.SourceType=NA;
  SegyTraceHeader.SourceEnergyDirectionMantissa=NA;
  SegyTraceHeader.SourceEnergyDirectionExponent=NA;
  SegyTraceHeader.SourceMeasurementMantissa=NA;
  SegyTraceHeader.SourceMeasurementExponent=NA;
  SegyTraceHeader.SourceMeasurementUnit=NA;
  
  SegyTraceHeader.UnassignedInt1=NA;
  SegyTraceHeader.UnassignedInt2=NA;
  
  
  