% TraceHeaderDef : Defines names, position, and precision for Trace Headers
%
% % To get a Matlab structure with trace header definitions call:
% STH==TraceHeaderDef;
% % To get a list fo trace header definision listed on the screen call:
% STH==TraceHeaderDef(1)
%
% See also: ReadSegyTraceHeaderValue, WriteSegyTraceHeaderValue
%
function STH=TraceHeaderDef(print);

if nargin>0
    STH=TraceHeaderDef;
    fn=fieldnames(STH);
    disp(sprintf('%4s %6s %s','POS','PREC','Traece Header Name'))
    for i=1:length(fn)
        disp(sprintf('%4d %6s %s',STH.(fn{i}).pos,STH.(fn{i}).precision,fn{i}))
    end
    return
end

%NEW
STH.TraceSequenceLine.pos=0; STH.TraceSequenceLine.precision='int32';
STH.TraceSequenceFile.pos=4; STH.TraceSequenceFile.precision='int32';
STH.FieldRecord.pos=8;       STH.FieldRecord.precision='int32';
STH.TraceNumber.pos=12;      STH.TraceNumber.precision='int32';
STH.EnergySourcePoint.pos=16;STH.EnergySourcePoint.precision='int32';
STH.cdp.pos=20;              STH.cdp.precision='int32';
STH.cdpTrace.pos=24;         STH.cdpTrace.precision='int32';

STH.TraceIdenitifactionCode.pos=28;STH.TraceIdenitifactionCode.precision='int16';
STH.NSummedTraces.pos=30;           STH.NSummedTraces.precision='int16';
STH.NStackedTraces.pos=32;          STH.NStackedTraces.precision='int16';
STH.DataUse.pos=34;                 STH.DataUse.precision='int16';

STH.offset.pos=36;STH.offset.precision='int32';
STH.ReceiverGroupElevation.pos=40;STH.ReceiverGroupElevation.precision='int32';
STH.SourceSurfaceElevation.pos=44;STH.SourceSurfaceElevation.precision='int32';
STH.SourceDepth.pos=48;STH.SourceDepth.precision='int32';
STH.ReceiverDatumElevation.pos=52;STH.ReceiverDatumElevation.precision='int32';
STH.SourceDatumElevation.pos=56;STH.SourceDatumElevation.precision='int32';
STH.SourceWaterDepth.pos=60;STH.SourceWaterDepth.precision='int32';
STH.GroupWaterDepth.pos=64;STH.GroupWaterDepth.precision='int32';

STH.ElevationScalar.pos=68;STH.ElevationScalar.precision='int16';
STH.SourceGroupScalar.pos=70;STH.SourceGroupScalar.precision='int16';

STH.SourceX.pos=72;STH.SourceX.precision='int32';
STH.SourceY.pos=76;STH.SourceY.precision='int32';
STH.GroupX.pos=80;STH.GroupX.precision='int32';
STH.GroupY.pos=84;STH.GroupY.precision='int32';

STH.CoordinateUnits.pos=88;STH.CoordinateUnits.precision='int16';
STH.WeatheringVelocity.pos=90;STH.WeatheringVelocity.precision='int16';
STH.SubWeatheringVelocity.pos=92;STH.SubWeatheringVelocity.precision='int16';
STH.SourceUpholeTime.pos=94;STH.SourceUpholeTime.precision='int16';
STH.GroupUpholeTime.pos=96;STH.GroupUpholeTime.precision='int16';
STH.SourceStaticCorrection.pos=98;STH.SourceStaticCorrection.precision='int16';
STH.GroupStaticCorrection.pos=100;STH.GroupStaticCorrection.precision='int16';
STH.TotalStaticApplied.pos=102;STH.TotalStaticApplied.precision='int16';
STH.LagTimeA.pos=104;STH.LagTimeA.precision='int16';
STH.LagTimeB.pos=106;STH.LagTimeB.precision='int16';
STH.DelayRecordingTime.pos=108;STH.DelayRecordingTime.precision='int16';
STH.MuteTimeStart.pos=110;STH.MuteTimeStart.precision='int16';
STH.MuteTimeEND.pos=112;STH.MuteTimeEND.precision='int16';

STH.ns.pos=114;STH.ns.precision='uint16';
STH.dt.pos=116;STH.dt.precision='uint16';

STH.GainType.pos=118;STH.GainType.precision='int16';
STH.InstrumentGainConstant.pos=120;STH.InstrumentGainConstant.precision='int16';
STH.InstrumentInitialGain.pos=122;STH.InstrumentInitialGain.precision='int16';
STH.Correlated.pos=124;STH.Correlated.precision='int16';
STH.SweepFrequenceStart.pos=126;STH.SweepFrequenceStart.precision='int16';
STH.SweepFrequenceEnd.pos=128;STH.SweepFrequenceEnd.precision='int16';
STH.SweepLength.pos=130;STH.SweepLength.precision='int16';
STH.SweepType.pos=132;STH.SweepType.precision='int16';
STH.SweepTraceTaperLengthStart.pos=134;STH.SweepTraceTaperLengthStart.precision='int16';
STH.SweepTraceTaperLengthEnd.pos=136;STH.SweepTraceTaperLengthEnd.precision='int16';
STH.TaperType.pos=138;STH.TaperType.precision='int16';
STH.AliasFilterFrequency.pos=140;STH.AliasFilterFrequency.precision='int16';
STH.AliasFilterSlope.pos=142;STH.AliasFilterSlope.precision='int16';
STH.NotchFilterFrequency.pos=144;STH.NotchFilterFrequency.precision='int16';
STH.NotchFilterSlope.pos=146;STH.NotchFilterSlope.precision='int16';
STH.LowCutFrequency.pos=148;STH.LowCutFrequency.precision='int16';
STH.HighCutFrequency.pos=150;STH.HighCutFrequency.precision='int16';
STH.LowCutSlope.pos=152;STH.LowCutSlope.precision='int16';
STH.HighCutSlope.pos=154;STH.HighCutSlope.precision='int16';
STH.YearDataRecorded.pos=156;STH.YearDataRecorded.precision='int16';
STH.DayOfYear.pos=158;STH.DayOfYear.precision='int16';
STH.HourOfDay.pos=160;STH.HourOfDay.precision='int16';
STH.MinuteOfHour.pos=162;STH.MinuteOfHour.precision='int16';
STH.SecondOfMinute.pos=164;STH.SecondOfMinute.precision='int16';
STH.TimeBaseCode.pos=166;STH.TimeBaseCode.precision='int16';
%STH.TimeBaseCodeText.pos=168;STH.TimeBaseCodeText.precision='int16';
STH.TraceWeightningFactor.pos=168;STH.TraceWeightningFactor.precision='int16';
STH.GeophoneGroupNumberRoll1.pos=170;STH.GeophoneGroupNumberRoll1.precision='int16';
STH.GeophoneGroupNumberFirstTraceOrigField.pos=172;STH.GeophoneGroupNumberFirstTraceOrigField.precision='int16';
STH.GeophoneGroupNumberLastTraceOrigField.pos=174;STH.GeophoneGroupNumberLastTraceOrigField.precision='int16';
STH.GapSize.pos=176;STH.GapSize.precision='int16';
STH.OverTravel.pos=178;STH.OverTravel.precision='int16';

STH.cdpX.pos=180;STH.cdpX.precision='int32';
STH.cdpY.pos=184;STH.cdpY.precision='int32';
STH.Inline3D.pos=188;STH.Inline3D.precision='int32';
STH.Crossline3D.pos=192;STH.Crossline3D.precision='int32';
STH.ShotPoint.pos=196;STH.ShotPoint.precision='int32';

STH.ShotPointScalar.pos=200;STH.ShotPointScalar.precision='int16';

STH.TraceValueMeasurementUnit.pos=202;STH.TraceValueMeasurementUnit.precision='int16';

STH.TransductionConstantMantissa.pos=204;STH.TransductionConstantMantissa.precision='int32';

STH.TransductionConstantPower.pos=208;STH.TransductionConstantPower.precision='int16';
STH.TransductionUnit.pos=210;STH.TransductionUnit.precision='int16';
STH.TraceIdentifier.pos=212;STH.TraceIdentifier.precision='int16';
STH.ScalarTraceHeader.pos=214;STH.ScalarTraceHeader.precision='int16';
STH.SourceType.pos=216;STH.SourceType.precision='int16';

STH.SourceEnergyDirectionMantissa.pos=218;STH.SourceEnergyDirectionMantissa.precision='int32';

STH.SourceEnergyDirectionExponent.pos=222;STH.SourceEnergyDirectionExponent.precision='int16';

STH.SourceMeasurementMantissa.pos=224;STH.SourceMeasurementMantissa.precision='in32';

STH.SourceMeasurementExponent.pos=228;STH.SourceMeasurementExponent.precision='int16';
STH.SourceMeasurementUnit.pos=230;STH.SourceMeasurementUnit.precision='int16';

STH.UnassignedInt1.pos=232;STH.UnassignedInt1.precision='int32';
STH.UnassignedInt2.pos=236;STH.UnassignedInt2.precision='int32';
