% GetSegyHeaderBasics : Default Segy Header Header settings
% 
% Call :
% Rev=GetSegyHeaderBasics
%
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
function Rev=GetSegyHeaderBasics;

%
%  TODO
%  Fixed Point ???
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revision 0
Rev(1).name='Revision 0 (1975)';
Rev(1).SegyFormatRevisionNumber=0;
% DataSampleFormat 
Rev(1).DataSampleFormat(1).name='4-byte IBM Floating Point';
Rev(1).DataSampleFormat(2).name='4-byte Fixed Point - (UNSUPPORTED)';
Rev(1).DataSampleFormat(3).name='2-byte Fixed Point - (UNSUPPORTED)';
Rev(1).DataSampleFormat(4).name='4-byte Fixed Point with gain - (UNSUPPORTED)';

Rev(1).DataSampleFormat(1).format='uint32'; %
Rev(1).DataSampleFormat(2).format='int32'; % 
Rev(1).DataSampleFormat(3).format='int16'; %
%Rev(1).DataSampleFormat(2).format=''; % 'uint32'; % 
%Rev(1).DataSampleFormat(3).format=''; % 'uint16'; %
Rev(1).DataSampleFormat(4).format=''; % 'uint32'; % 

Rev(1).DataSampleFormat(1).bps=32;
Rev(1).DataSampleFormat(2).bps=32;
Rev(1).DataSampleFormat(3).bps=16;
%Rev(1).DataSampleFormat(4).bps=32;
%Rev(1).DataSampleFormat(3).bps=0; % NOT IMPLEMENETD
Rev(1).DataSampleFormat(4).bps=0; % NOT IMPLEMENETD

% TraceSorting
Rev(1).TraceSorting(1).name='as recorded (no sorting)';
Rev(1).TraceSorting(2).name='CDP ensemble';
Rev(1).TraceSorting(3).name='single fold continous profile';
Rev(1).TraceSorting(4).name='horizontally stacked';
Rev(1).TraceSorting(1).value=1;
Rev(1).TraceSorting(2).value=2;
Rev(1).TraceSorting(3).value=3;
Rev(1).TraceSorting(4).value=4;

% SweepType
Rev(1).SweepType(1).value=1;
Rev(1).SweepType(1).name='Linear';
Rev(1).SweepType(2).value=2;
Rev(1).SweepType(2).name= 'Parabolic';
Rev(1).SweepType(3).value=3;
Rev(1).SweepType(3).name='Exponential';
Rev(1).SweepType(4).value=4;
Rev(1).SweepType(4).name= 'Other';

% TaperType
Rev(1).TaperType(1).value=1;
Rev(1).TaperType(1).name='Linear';
Rev(1).TaperType(2).value=2;
Rev(1).TaperType(2).name= 'Cos2';
Rev(1).TaperType(3).value=3;
Rev(1).TaperType(3).name='Other';

% CorrelatedDataTraces
Rev(1).CorrelatedDataTraces(1).value=1;
Rev(1).CorrelatedDataTraces(1).name='No';
Rev(1).CorrelatedDataTraces(2).value=2;
Rev(1).CorrelatedDataTraces(2).name= 'Yes';

% BinaryGain
Rev(1).BinaryGain(1).value=1;
Rev(1).BinaryGain(1).name='Yes';
Rev(1).BinaryGain(2).value=2;
Rev(1).BinaryGain(2).name= 'No';

% AmplitudeRecoveryMethod
Rev(1).AmplitudeRecoveryMethod(1).value=1;
Rev(1).AmplitudeRecoveryMethod(1).name='none';
Rev(1).AmplitudeRecoveryMethod(2).value=2;
Rev(1).AmplitudeRecoveryMethod(2).name= 'Spherical Divergence';
Rev(1).AmplitudeRecoveryMethod(3).value=3;
Rev(1).AmplitudeRecoveryMethod(3).name= 'AGC';
Rev(1).AmplitudeRecoveryMethod(4).value=4;
Rev(1).AmplitudeRecoveryMethod(4).name= 'Other';

% MeasurementSystem
Rev(1).MeasurementSystem(1).value=1;
Rev(1).MeasurementSystem(1).name='Meters';
Rev(1).MeasurementSystem(2).value=2;
Rev(1).MeasurementSystem(2).name= 'Feets';

% ImpulseSignalPolarity
Rev(1).ImpulseSignalPolarity(1).value=1;
Rev(1).ImpulseSignalPolarity(1).name='Negative number';
Rev(1).ImpulseSignalPolarity(2).value=2;
Rev(1).ImpulseSignalPolarity(2).name= 'Positive Number';

% VibratoryPolarityCode
Rev(1).VibratoryPolarityCode(1).value=1;
Rev(1).VibratoryPolarityCode(1).name='337.5 - 22.5';
Rev(1).VibratoryPolarityCode(2).value=2;
Rev(1).VibratoryPolarityCode(2).name= '22.5-67.5';
Rev(1).VibratoryPolarityCode(3).value=3;
Rev(1).VibratoryPolarityCode(3).name= '67.5-112.5';
Rev(1).VibratoryPolarityCode(4).value=4;
Rev(1).VibratoryPolarityCode(4).name= '112.5-157.5';
Rev(1).VibratoryPolarityCode(5).value=5;
Rev(1).VibratoryPolarityCode(5).name= '157.5-202.5';
Rev(1).VibratoryPolarityCode(6).value=6;
Rev(1).VibratoryPolarityCode(6).name= '202.5-247.5';
Rev(1).VibratoryPolarityCode(7).value=7;
Rev(1).VibratoryPolarityCode(7).name= '247.5-292.5';
Rev(1).VibratoryPolarityCode(8).value=8;
Rev(1).VibratoryPolarityCode(8).name= '292.5-337.5';

% FixedLengthTraceFlag
Rev(1).FixedLengthTraceFlag(1).value=0;
Rev(1).FixedLengthTraceFlag(1).name='Varying Trace Length';
Rev(1).FixedLengthTraceFlag(2).value=1;
Rev(1).FixedLengthTraceFlag(2).name= 'Fixed Trace Length';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revision 1
Rev(2).name='Revision 1 (2002)';
Rev(2).SegyFormatRevisionNumber=100;
% DataSampleFormat
Rev(2).DataSampleFormat(1).format='uint32'; %
Rev(2).DataSampleFormat(2).format='int32'; % 
Rev(2).DataSampleFormat(3).format='int16'; %
Rev(2).DataSampleFormat(4).format=''; % 
Rev(2).DataSampleFormat(5).format='float32'; %
Rev(2).DataSampleFormat(6).format=''; % 
Rev(2).DataSampleFormat(7).format=''; % 
Rev(2).DataSampleFormat(8).format='int8'; % 

Rev(2).DataSampleFormat(1).bps=32;
Rev(2).DataSampleFormat(2).bps=32;
Rev(2).DataSampleFormat(3).bps=16;
Rev(2).DataSampleFormat(4).bps=0;
Rev(2).DataSampleFormat(5).bps=32;
Rev(2).DataSampleFormat(6).bps=0;
Rev(2).DataSampleFormat(7).bps=0;
Rev(2).DataSampleFormat(8).bps=8;

Rev(2).DataSampleFormat(1).name='4-byte IBM Floating Point';
Rev(2).DataSampleFormat(2).name='4-byte two''s complement';
Rev(2).DataSampleFormat(3).name='2-byte two''s complement';
Rev(2).DataSampleFormat(4).name='4-byte fixed-point with gain (obsolete) - NOT IMPLEMENTED';
Rev(2).DataSampleFormat(5).name='4-byte IEEE floating-point';
Rev(2).DataSampleFormat(6).name='DSF Not currently used';
Rev(2).DataSampleFormat(7).name='DSF Not currently used';
Rev(2).DataSampleFormat(8).name='1-byte, two''s complement integer';

% TraceSorting
Rev(2).TraceSorting(1).value=-1;
Rev(2).TraceSorting(1).name='Other (should be explained in user Extended Textual File Header stanza)';
Rev(2).TraceSorting(2).value=0;
Rev(2).TraceSorting(2).name= 'Unknown';
Rev(2).TraceSorting(3).value=1;
Rev(2).TraceSorting(3).name= 'As recorded (no sorting)';
Rev(2).TraceSorting(4).value=2;
Rev(2).TraceSorting(4).name= 'CDP ensemble';
Rev(2).TraceSorting(5).value=3;
Rev(2).TraceSorting(5).name= 'Single fold continuous profile';
Rev(2).TraceSorting(6).value=4;
Rev(2).TraceSorting(6).name= 'Horizontally stacked';
Rev(2).TraceSorting(7).value=5;
Rev(2).TraceSorting(7).name= 'Common source';
Rev(2).TraceSorting(8).value=6;
Rev(2).TraceSorting(8).name= 'Common receiver';
Rev(2).TraceSorting(9).value=7;
Rev(2).TraceSorting(9).name= 'Common offset';
Rev(2).TraceSorting(10).value=8;
Rev(2).TraceSorting(10).name= 'Common mid-point';

% SweepType
Rev(2).SweepType(1).value=1;
Rev(2).SweepType(1).name='Linear';
Rev(2).SweepType(2).value=2;
Rev(2).SweepType(2).name= 'Parabolic';
Rev(2).SweepType(3).value=3;
Rev(2).SweepType(3).name='Exponential';
Rev(2).SweepType(4).value=4;
Rev(2).SweepType(4).name= 'Other';

% TaperType
Rev(2).TaperType(1).value=1;
Rev(2).TaperType(1).name='Linear';
Rev(2).TaperType(2).value=2;
Rev(2).TaperType(2).name= 'Cos2';
Rev(2).TaperType(3).value=3;
Rev(2).TaperType(3).name='Other';

% CorrelatedDataTraces
Rev(2).CorrelatedDataTraces(1).value=1;
Rev(2).CorrelatedDataTraces(1).name='No';
Rev(2).CorrelatedDataTraces(2).value=2;
Rev(2).CorrelatedDataTraces(2).name= 'Yes';

% BinaryGain
Rev(2).BinaryGain(1).value=1;
Rev(2).BinaryGain(1).name='Yes';
Rev(2).BinaryGain(2).value=2;
Rev(2).BinaryGain(2).name= 'No';

% AmplitudeRecoveryMethod
Rev(2).AmplitudeRecoveryMethod(1).value=1;
Rev(2).AmplitudeRecoveryMethod(1).name='none';
Rev(2).AmplitudeRecoveryMethod(2).value=2;
Rev(2).AmplitudeRecoveryMethod(2).name= 'Spherical Divergence';
Rev(2).AmplitudeRecoveryMethod(3).value=3;
Rev(2).AmplitudeRecoveryMethod(3).name= 'AGC';
Rev(2).AmplitudeRecoveryMethod(4).value=4;
Rev(2).AmplitudeRecoveryMethod(4).name= 'Other';

% MeasurementSystem
Rev(2).MeasurementSystem(1).value=1;
Rev(2).MeasurementSystem(1).name='Meters';
Rev(2).MeasurementSystem(2).value=2;
Rev(2).MeasurementSystem(2).name= 'Feets';

% ImpulseSignalPolarity
Rev(2).ImpulseSignalPolarity(1).value=1;
Rev(2).ImpulseSignalPolarity(1).name='Negative number';
Rev(2).ImpulseSignalPolarity(2).value=2;
Rev(2).ImpulseSignalPolarity(2).name= 'Positive Number';

% VibratoryPolarityCode
Rev(2).VibratoryPolarityCode(1).value=1;
Rev(2).VibratoryPolarityCode(1).name='337.5 - 22.5';
Rev(2).VibratoryPolarityCode(2).value=2;
Rev(2).VibratoryPolarityCode(2).name= '22.5-67.5';
Rev(2).VibratoryPolarityCode(3).value=3;
Rev(2).VibratoryPolarityCode(3).name= '67.5-112.5';
Rev(2).VibratoryPolarityCode(4).value=4;
Rev(2).VibratoryPolarityCode(4).name= '112.5-157.5';
Rev(2).VibratoryPolarityCode(5).value=5;
Rev(2).VibratoryPolarityCode(5).name= '157.5-202.5';
Rev(2).VibratoryPolarityCode(6).value=6;
Rev(2).VibratoryPolarityCode(6).name= '202.5-247.5';
Rev(2).VibratoryPolarityCode(7).value=7;
Rev(2).VibratoryPolarityCode(7).name= '247.5-292.5';
Rev(2).VibratoryPolarityCode(8).value=8;
Rev(2).VibratoryPolarityCode(8).name= '292.5-337.5';

% FixedLengthTraceFlag
Rev(2).FixedLengthTraceFlag(1).value=0;
Rev(2).FixedLengthTraceFlag(1).name='Varying Trace Length';
Rev(2).FixedLengthTraceFlag(2).value=1;
Rev(2).FixedLengthTraceFlag(2).name= 'Fixed Trace Length';

