function varargout = GUIEditSegyHeader(varargin)
% GUIEditSegyHeader : GUI for editing the SGY header
%
% Call : 
%   SH=GUIEditSegyHeader(SH);
%
%
% Example : 
%   [Data,STH,SH]=ReadSegy('841_m.sgy');
%   SH=GUIEditSegyHeader(SH);
%
%

% GUIEDITSEGYHEADER Application M-file for GUIEditSegyTraceHeader.fig
%    FIG = GUIEDITSEGYHEADER launch GUIEditSegyTraceHeader GUI.
%    GUIEDITSEGYHEADER('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 03-Mar-2003 15:33:07

if (nargin == 0)|(isstruct(varargin{1}))  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

    if nargin>0,
       data = guihandles(fig); % initialize it to contain handles
       data.SegyHeader=varargin{1};
       data.SegyHeaderOrig=varargin{1};
       guidata(fig,data);
       
       set(fig,'HandleVisibility','On')
       GUIEditSegyHeader('actionSetHeader',fig,handles)
       GUIEditSegyHeader('actionUpdateHeader',fig,handles)
       set(fig,'HandleVisibility','CallBack')

    end
    
	% Wait for callbacks to run and window to be dismissed:
	uiwait(fig);

	if nargout > 0
	  data=guidata(fig);
	  if isstruct(data.SegyHeader)
	    varargout{1}=data.SegyHeader;
	  else
	    varargout{1} = fig;
	  end
    end
    
    
    close(handles.figure1);
	set(handles.figure1,'Visible','off')
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', eSweepFrequencyStart.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, eSweepFrequencyStart.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.

% --------------------------------------------------------------------

% --------------------------------------------------------------------
function varargout = figure1_ResizeFcn(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = pbDone_Callback(h, eventdata, handles, varargin)
uiresume;

% --------------------------------------------------------------------
function varargout = pbCancel_Callback(h, eventdata, handles, varargin)
data=guidata(h);
data.SegyHeader=data.SegyHeaderOrig;
guidata(h,data);
uiresume;


clear name;
function actionSetHeader(h,handles)
data=guidata(h);
SH=data.SegyHeader;

% Revision
revision=SH.SegyFormatRevisionNumber;
r=1; % Default old revisoin
if revision>50, r=2;end
if revision==0, r=1;end
for i=1:length(SH.Rev), Revision{i}=SH.Rev(i).name;end
set(handles.popSegyFormatRevision,'String',Revision);

clear name;
% DataSampleFormat
DSF=SH.Rev(r).DataSampleFormat;
for i=1:length(DSF), name{i}=DSF(i).name;end
set(handles.popDataSampleFormat,'String',name);

clear name;
% TraceSorting
TS=SH.Rev(r).TraceSorting;
for i=1:length(TS), name{i}=TS(i).name;end
set(handles.popTraceSorting,'String',name);

clear name;
% Sweep Type
ST=SH.Rev(r).SweepType;
for i=1:length(ST), name{i}=ST(i).name;end
set(handles.popSweepType,'String',name);

clear name;
% BinaryGain
BG=SH.Rev(r).BinaryGain;
for i=1:length(BG), name{i}=BG(i).name;end
set(handles.popBinaryGain,'String',name);

clear name;
% TaperType
TT=SH.Rev(r).TaperType;
for i=1:length(TT), name{i}=TT(i).name;end
set(handles.popTaperType,'String',name);

clear name;
% CorrelatedDataTraces
CD=SH.Rev(r).CorrelatedDataTraces;
for i=1:length(CD), name{i}=CD(i).name;end
set(handles.popCorrelatedDataTraces,'String',name);

clear name;
% AmplitudeRecoveryMethod
ARM=SH.Rev(r).AmplitudeRecoveryMethod;
for i=1:length(ARM), name{i}=ARM(i).name;end
set(handles.popAmplitudeRecoveryMethod,'String',name);

clear name;
% MeasurementSystem
M=SH.Rev(r).MeasurementSystem;
for i=1:length(M), name{i}=M(i).name;end
set(handles.popMeasurementSystem,'String',name);

clear name;
% ImpulseSignalPolarity
ISR=SH.Rev(r).ImpulseSignalPolarity;
for i=1:length(ISR), name{i}=ISR(i).name;end
set(handles.popImpulseSignalPolarity,'String',name);

clear name;
% popVibratoryPolarityCode
VPC=SH.Rev(r).VibratoryPolarityCode;
for i=1:length(VPC), name{i}=VPC(i).name;end
set(handles.popVibratoryPolarityCode,'String',name);

clear name;
% popFixedLengthTraceFlag
FLT=SH.Rev(r).FixedLengthTraceFlag;
for i=1:length(FLT), name{i}=FLT(i).name;end
set(handles.popFixedLengthTraceFlag,'String',name);


% EDIT BOXES
set(handles.eJob,'String',SH.Job);
set(handles.eLine,'String',SH.Line);
set(handles.eReel,'String',SH.Reel);
set(handles.eDataTracePerEnsemble,'String',SH.DataTracePerEnsemble);
set(handles.eAuxiliaryTracePerEnsemble,'String',SH.AuxiliaryTracePerEnsemble);
set(handles.edt,'String',SH.dt);
set(handles.edtOrig,'String',SH.dtOrig);
set(handles.ens,'String',SH.ns);
set(handles.ensOrig,'String',SH.nsOrig);
set(handles.eEnsembleFold,'String',SH.EnsembleFold);
set(handles.eVerticalSumCode,'String',SH.VerticalSumCode);

set(handles.eSweepFrequencyStart,'String',SH.SweepFrequencyStart);
set(handles.eSweepFrequencyEnd,'String',SH.SweepFrequencyEnd);
set(handles.eSweepLength,'String',SH.SweepLength);
set(handles.eSweepChannel,'String',SH.SweepChannel);
set(handles.eSweepTaperlengthStart,'String',SH.SweepTaperlengthStart);
set(handles.eSweepTaperLengthEnd,'String',SH.SweepTaperLengthEnd);




%---------------------------------------------------------
function r=actionGetrevision(h);
data=guidata(h);
SH=data.SegyHeader;
% Revision
revision=SH.SegyFormatRevisionNumber;
r=1; % Default old revisoin
if revision>50, r=2;end
if revision==0, r=1;end


%---------------------------------------------------------
function actionUpdateHeader(h,handles)
data=guidata(h);

SH=data.SegyHeader;

% Revision
revision=SH.SegyFormatRevisionNumber;
r=1; % Default old revisoin
if revision>50, r=2;end
if revision==0, r=1;end
set(handles.popSegyFormatRevision,'value',r);

% DataSampleFormat
DSF=SH.Rev(r).DataSampleFormat;
if length(get(handles.popDataSampleFormat,'String'))<SH.DataSampleFormat, 
    SH.DataSampleFormat=1;
    disp([mfilename,' : Using DataSampleFormat ',num2str(SH.DataSampleFormat)])
end
set(handles.popDataSampleFormat,'value',SH.DataSampleFormat);

% MAYBE THE NEXT LINES COULD BE USED FOR MORE HEADERS
% TraceSorting
TraceSortingSelect=1;
for i=1:length(SH.Rev(r).TraceSorting)
    if SH.TraceSorting==SH.Rev(r).TraceSorting(i).value;
        TraceSortingSelect=i;
        SH.TraceSorting;
    end
end

if TraceSortingSelect>length(get(handles.popTraceSorting,'String'));
    TraceSortingSelect=1;
end
set(handles.popTraceSorting,'value',TraceSortingSelect);

SweepTypeSelect=1;
for i=1:length(SH.Rev(r).SweepType)
    if SH.SweepType==SH.Rev(r).SweepType(i).value;
        SweepTypeSelect=i;
        SH.SweepType;
    end
end
if SweepTypeSelect>length(get(handles.popSweepType,'String'));
    SweepTypeSelect=1;
end
set(handles.popSweepType,'value',SweepTypeSelect);





data.SegyHeader=SH;
guidata(h,data)


% --------------------------------------------------------------------
function varargout = popSegyFormatRevision_Callback(h, eventdata, handles, varargin)
data=guidata(h);
sel=get(h,'value');
data.SegyHeader.SegyFormatRevisionNumber =  data.SegyHeader.Rev(sel).SegyFormatRevisionNumber;

guidata(h,data);
actionSetHeader(h,handles)
actionUpdateHeader(h,handles)


% --------------------------------------------------------------------
function varargout = popDataSampleFormat_Callback(h, eventdata, handles, varargin)
data=guidata(h);
sel=get(h,'value');
data.SegyHeader.DataSampleFormat=get(h,'value');
guidata(h,data);


% --------------------------------------------------------------------
function varargout = popTraceSorting_Callback(h, eventdata, handles, varargin)
data=guidata(h);
TraceSortingSelect=get(h,'value');
r=actionGetrevision(h);
data.SegyHeader.TraceSorting=data.SegyHeader.Rev(r).TraceSorting(TraceSortingSelect).value;
guidata(h,data);

function popSweepType_Callback(h, eventdata, handles)
data=guidata(h);
SweepTypeSelect=get(h,'value');
r=actionGetrevision(h);
data.SegyHeader.SweepType=data.SegyHeader.Rev(r).SweepType(SweepTypeSelect).value;
guidata(h,data);


function popTaperType_Callback(h, eventdata, handles)
data=guidata(h);
TaperTypeSelect=get(h,'value');
r=actionGetrevision(h);
data.SegyHeader.TaperType=data.SegyHeader.Rev(r).TaperType(TaperTypeSelect).value;
guidata(h,data);
actionUpdateHeader(h,handles)

function popCorrelatedDataTraces_Callback(h, eventdata, handles)
data=guidata(h);
CorrelatedDataTracesSelect=get(h,'value');
r=actionGetrevision(h);
data.SegyHeader.CorrelatedDataTraces=data.SegyHeader.Rev(r).CorrelatedDataTraces(CorrelatedDataTracesSelect).value;
guidata(h,data);

function popBinaryGain_Callback(h, eventdata, handles)
data=guidata(h);
BinaryGainSelect=get(h,'value');
r=actionGetrevision(h);
data.SegyHeader.BinaryGain=data.SegyHeader.Rev(r).BinaryGain(BinaryGainSelect).value;
guidata(h,data);



function popAmplitudeRecoveryMethod_Callback(h, eventdata, handles)
data=guidata(h);
AmplitudeRecoveryMethodSelect=get(h,'value');
r=actionGetrevision(h);
data.SegyHeader.AmplitudeRecoveryMethod=data.SegyHeader.Rev(r).AmplitudeRecoveryMethod(AmplitudeRecoveryMethodSelect).value;
guidata(h,data);

function varargout = popMeasurementSystem_Callback(h, eventdata, handles, varargin)
data=guidata(h);
MSSelect=get(h,'value');
r=actionGetrevision(h);
data.SegyHeader.MeasurementSystem=data.SegyHeader.Rev(r).MeasurementSystem(MSSelect).value;
guidata(h,data);

function popImpulseSignalPolarity_Callback(h, eventdata, handles)
data=guidata(h);
ImpulseSignalPolaritySelect=get(h,'value');
r=actionGetrevision(h);
data.SegyHeader.ImpulseSignalPolarity=data.SegyHeader.Rev(r).ImpulseSignalPolarity(ImpulseSignalPolaritySelect).value;
guidata(h,data);

function popVibratoryPolarityCode_Callback(h, eventdata, handles)
data=guidata(h);
VibratoryPolarityCodeSelect=get(h,'value');
r=actionGetrevision(h);
data.SegyHeader.VibratoryPolarityCode=data.SegyHeader.Rev(r).VibratoryPolarityCode(VibratoryPolarityCodeSelect).value;
guidata(h,data);

% --- Executes on selection change in popFixedLengthTraceFlag.
function popFixedLengthTraceFlag_Callback(h, eventdata, handles)
data=guidata(h);
FixedLengthTraceFlagSelect=get(h,'value');
r=actionGetrevision(h);
data.SegyHeader.FixedLengthTraceFlag=data.SegyHeader.Rev(r).FixedLengthTraceFlag(FixedLengthTraceFlagSelect).value;
guidata(h,data);



% --- Executes on button press in pbEditTextualFileHeader.
function pbEditTextualFileHeader_Callback(h, eventdata, handles)

data=guidata(h);
data.SegyHeader=GUIEditTextualFileHeader(data.SegyHeader);
guidata(h,data);




function eJob_Callback(hObject, eventdata, handles)


function eLine_Callback(hObject, eventdata, handles)


function eReel_Callback(hObject, eventdata, handles)

function eDataTracePerEnsemble_Callback(hObject, eventdata, handles)

function eAuxiliaryTracePerEnsemble_Callback(hObject, eventdata, handles)

function edt_Callback(hObject, eventdata, handles)

function edtOrig_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ens_CreateFcn(hObject, eventdata, handles)

% --- Executes on button press in ens.
function ens_Callback(hObject, eventdata, handles)

% --- Executes on button press in ensOrig.
function ensOrig_Callback(hObject, eventdata, handles)


% --- Executes on button press in eEnsembleFold.
function eEnsembleFold_Callback(hObject, eventdata, handles)

% --- Executes on button press in eVerticalSumCode.
function eVerticalSumCode_Callback(hObject, eventdata, handles)

% --- Executes on button press in eSweepFrequencyStart.
function eSweepFrequencyStart_Callback(hObject, eventdata, handles)

function eSweepFrequencyEnd_Callback(hObject, eventdata, handles)

function eSweepLength_Callback(hObject, eventdata, handles)

function eSweepChannel_Callback(hObject, eventdata, handles)

function eSweepTaperlengthStart_Callback(hObject, eventdata, handles)

function eSweepTaperLengthEnd_Callback(hObject, eventdata, handles)
