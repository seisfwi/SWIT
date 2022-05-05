function varargout = DIAOpenSEGY(varargin)
% DIAOPENSEGY Application M-file for DIAOpenSEGY.fig
%    FIG = DIAOPENSEGY launch DIAOpenSEGY GUI.
%    DIAOPENSEGY('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 12-Jun-2002 10:25:58

if ((nargin == 0)|isstruct(varargin{1}))  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

    data = guihandles(fig); % initialize it to contain handles
    data.JUMP.jump=5;
    data.JUMP.enable=0;
    data.TRACEHEADER.name='offset';
    data.TRACEHEADER.min=1;
    data.TRACEHEADER.max=1e+9;
    data.TRACEHEADER.enable=0;
    data.TIMERANGE.min=1;
    data.TIMERANGE.max=5;
    data.TIMERANGE.enable=0;

    if nargin==0;
      data.pathname='';
      data.filename='';
      guidata(fig,data);
    else
      try;data.JUMP=varargin{1}.JUMP;end
      try;data.TRACEHEADER=varargin{1}.TRACEHEADER;end
      try;data.TIMERANGE=varargin{1}.TIMERANGE;end
      if isfield(varargin{1},'pathname');data.pathname=varargin{1}.pathname;end;
      if isfield(varargin{1},'filename');data.filename=varargin{1}.filename;end
    end
    guidata(fig,data);    
      
    set(fig,'HandleVisibility','On')
    DIAOpenSEGY('actionUpdateGUI',fig,handles)
    set(fig,'HandleVisibility','CallBack')
	
    % Wait for callbacks to run and window to be dismissed:
	uiwait(fig);

    if nargout > 0
       data=guidata(fig);
       SEGYFILE.JUMP=data.JUMP;
       SEGYFILE.TRACEHEADER=data.TRACEHEADER;
       SEGYFILE.TIMERANGE=data.TIMERANGE;
       if isfield(data,'SegyHeader')
         SEGYFILE.SegyHeader=data.SegyHeader;
       end
       if isfield(data,'filename'),SEGYFILE.filename=data.filename;end
       if isfield(data,'pathname'),SEGYFILE.pathname=data.pathname;end
       varargout{1} = SEGYFILE;
    end
    

    close(handles.figure1);
    try;set(handles.figure1,'Visible','off');end

    
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
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
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

function actionUpdateGUI(h,handles)
data=guidata(h);

% RESET FILEMAME,PATHNAME if FILE DOES NOT EXIST
if isfield(data,'filename')  
  if exist([data.pathname,data.filename])~=2
      data=rmfield(data,'filename');
      data=rmfield(data,'pathname');
  end
end
    
% CHECK THAT A FILENAME I SET, and disbale if not.
if isfield(data,'filename'),    
    cbEnable='on';
    set(handles.eFilename,'String',[data.pathname,data.filename]);
else    
    cbEnable='off';
    set(handles.eFilename,'String','');
end
set(findobj('style','checkbox'),'enable',cbEnable)


% JUMP
set(handles.cbJump,'Value',data.JUMP.enable);
data.Jump.value=get(handles.cbJump,'Value');
if ((get(handles.cbJump,'Value')==1)&(isfield(data,'filename'))), 
    enableJump='on'; 
    data.JUMP.enable=1;
else, 
    enableJump='off'; 
    data.JUMP.enable=0;
end
set(handles.eJump,'Enable',enableJump)    
set(handles.eJump,'String',data.JUMP.jump) 

% TraceHeader
set(handles.cbTraceHeader,'Value',data.TRACEHEADER.enable);
data.TraceHeader.value=get(handles.cbTraceHeader,'Value');
if ((get(handles.cbTraceHeader,'Value')==1)&(isfield(data,'filename'))), 
    enableTraceHeader='on';
    data.TRACEHEADER.enable=1;
else, 
    enableTraceHeader='off'; 
    data.TRACEHEADER.enable=0;
end
hnames=get(handles.popTraceHeader,'string');
selected=get(handles.popTraceHeader,'value');
for i=1:length(hnames)
    if strcmp(hnames{i},data.TRACEHEADER.name),
        set(handles.popTraceHeader,'Value',i);
    end
end
set(handles.popTraceHeader,'Enable',enableTraceHeader)  
set(handles.eTraceHeaderMin,'Enable',enableTraceHeader) 
set(handles.eTraceHeaderMin,'String',data.TRACEHEADER.min)
set(handles.eTraceHeaderMax,'Enable',enableTraceHeader)    
set(handles.eTraceHeaderMax,'String',data.TRACEHEADER.max)
    
% TimeRange
set(handles.cbTimeRange,'Value',data.TIMERANGE.enable);
data.TimeRange.value=get(handles.cbTimeRange,'Value');
if ((get(handles.cbTimeRange,'Value')==1)&(isfield(data,'filename'))), 
    enableTimeRange='on'; 
    data.TIMERANGE.enable=1;
else, 
    data.TIMERANGE.enable=0;
    enableTimeRange='off'; 
end
set(handles.eTimeRangeMin,'Enable',enableTimeRange)    
set(handles.eTimeRangeMin,'String',data.TIMERANGE.min)
set(handles.eTimeRangeMax,'Enable',enableTimeRange)    
set(handles.eTimeRangeMax,'String',data.TIMERANGE.max)

guidata(h,data)

% --------------------------------------------------------------------
function varargout = eFilename_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = pbSelectFile_Callback(h, eventdata, handles, varargin)
data=guidata(h);
[filename, pathname] = uigetfile('*.su;*.SU;*.sgy;*.SGY;*.seg*;*.SEGY', 'Pick a SEGY or SU file');

    if isequal(filename,0)|isequal(pathname,0)
       disp('File not found')
       return
    else
       disp(['File ', pathname, filename, ' found'])
       set(handles.eFilename,'String',[pathname,filename])
       data.pathname=pathname;
       data.filename=filename;
    end
guidata(h,data)
actionUpdateGUI(h,handles)



% --------------------------------------------------------------------
function varargout = cbJump_Callback(h, eventdata, handles, varargin)
data=guidata(h);
data.JUMP.enable=get(h,'Value');
guidata(h,data);
actionUpdateGUI(h,handles)



% --------------------------------------------------------------------
function varargout = eJump_Callback(h, eventdata, handles, varargin)
data=guidata(h);
if ~isempty(str2num(get(h,'String')));data.JUMP.jump=str2num(get(h,'String'));end
guidata(h,data);
actionUpdateGUI(h,handles)




% --------------------------------------------------------------------
function varargout = cbTraceHeader_Callback(h, eventdata, handles, varargin)
data=guidata(h);
data.TRACEHEADER.enable=get(h,'Value');
guidata(h,data);
actionUpdateGUI(h,handles)



% --------------------------------------------------------------------
function varargout = popTraceHeader_Callback(h, eventdata, handles, varargin)
data=guidata(h);
hnames=get(h,'string');
selected=get(h,'value');
data.TRACEHEADER.name=hnames{selected};
guidata(h,data);
actionUpdateGUI(h,handles)



% --------------------------------------------------------------------
function varargout = eTraceHeaderMin_Callback(h, eventdata, handles, varargin)
data=guidata(h);
if ~isempty(str2num(get(h,'String')));data.TRACEHEADER.min=str2num(get(h,'String'));end
if data.TRACEHEADER.min>data.TRACEHEADER.max
    data.TRACEHEADER.min=data.TRACEHEADER.max;
end
guidata(h,data);
actionUpdateGUI(h,handles)



% --------------------------------------------------------------------
function varargout = eTraceHeaderMax_Callback(h, eventdata, handles, varargin)
data=guidata(h);
if ~isempty(str2num(get(h,'String')));data.TRACEHEADER.max=str2num(get(h,'String'));end
if data.TRACEHEADER.max<data.TRACEHEADER.min
    data.TRACEHEADER.max=data.TRACEHEADER.min;
end
guidata(h,data);
actionUpdateGUI(h,handles)




% --------------------------------------------------------------------
function varargout = cbTimeRange_Callback(h, eventdata, handles, varargin)
data=guidata(h);
data.TIMERANGE.enable=get(h,'Value');
guidata(h,data);
actionUpdateGUI(h,handles)



% --------------------------------------------------------------------
function varargout = eTimeRangeMin_Callback(h, eventdata, handles, varargin)
data=guidata(h);
if ~isempty(str2num(get(h,'String')));data.TIMERANGE.min=str2num(get(h,'String'));end
if data.TIMERANGE.min>data.TIMERANGE.max;   
    data.TIMERANGE.min=data.TIMERANGE.max;
end
guidata(h,data);
actionUpdateGUI(h,handles)




% --------------------------------------------------------------------
function varargout = eTimeRangeMax_Callback(h, eventdata, handles, varargin)
data=guidata(h);
if ~isempty(str2num(get(h,'String')));data.TIMERANGE.max=str2num(get(h,'String'));end
if data.TIMERANGE.max<data.TIMERANGE.min
    data.TIMERANGE.max=data.TIMERANGE.min;
end
guidata(h,data);
actionUpdateGUI(h,handles)



% --------------------------------------------------------------------
function varargout = pbReadSegy_Callback(h, eventdata, handles, varargin)
uiresume




% --------------------------------------------------------------------
function varargout = pbCancel_Callback(h, eventdata, handles, varargin)
data=guidata(h);
if isfield(data,'filename'), data=rmfield(data,'filename'); end
if isfield(data,'pathname'), data=rmfield(data,'pathname'); end
guidata(h,data);
uiresume



% --------------------------------------------------------------------
function varargout = pbGUIEditSegyHeader_Callback(h, eventdata, handles, varargin)
data=guidata(h);

% GET SEGY HEDER FROM FILE, UNLESS WE ALLREADY GOT IT
if isfield(data,'SegyHeader')==0
  data.SegyHeader=GetSegyHeader([data.pathname,data.filename]);
end

% EDIT THE SEGY HEADER
data.SegyHeader=GUIEditSegyHeader(data.SegyHeader);

guidata(h,data);
