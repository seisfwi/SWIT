function varargout = GUIEditSegyTraceHeader(varargin)
% GUIPlotXY Application M-file for GUIEditSegyTraceHeader.fig
%    FIG = GUIPLOTXY launch GUIEditSegyTraceHeader GUI.
%    GUIPLOTXY('callback_name', ...) invoke the named callback.
%
% OBSOLETE ??
%
  
  
% Last Modified by GUIDE v2.0 16-Jun-2002 11:41:23

if (nargin == 0)|(isstruct(varargin{1}))  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));


	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

    if nargin>0,
       data = guihandles(fig); % initialize it to contain handles
       data.H=varargin{1};
       data.Hname=fieldnames(data.H(1));
 
       data.NShowHeaders=7; % NUMBER OF HEADER VALUES TO LIST
       data.ShowNTraces=10; % NUMBER OF TRACES TO LIST 
       data.StartTrace=1;
       
       set(handles.popX,'String',data.Hname);
       set(handles.popY,'String',data.Hname);
       
       set(handles.popX,'Value',2);
       set(handles.popY,'Value',4);

       guidata(fig,data);
       
       set(fig,'HandleVisibility','On')
       GUIPlotXY('actionPlotData',fig,handles)
       set(fig,'HandleVisibility','CallBack')

    end
    
	% Wait for callbacks to run and window to be dismissed:
	uiwait(fig);

	if nargout > 0
	  data=guidata(fig);
	  if isstruct(data.H)
	    varargout{1}=data.H
	  else
	    varargout{1} = fig;
	  end
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

% --------------------------------------------------------------------

% --------------------------------------------------------------------
function varargout = popX_Callback(h, eventdata, handles, varargin)
actionPlotData(h,handles)

% --------------------------------------------------------------------
function varargout = popY_Callback(h, eventdata, handles, varargin)
actionPlotData(h,handles)

% --------------------------------------------------------------------
function varargout = figure1_ResizeFcn(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = pbDone_Callback(h, eventdata, handles, varargin)
uiresume;


function actionPlotData(h,handles)
data=guidata(h);
Xv=get(handles.popX,'value');
Yv=get(handles.popY,'value');

Xheader=data.Hname(Xv);
Yheader=data.Hname(Yv);

X=[data.H.(char(Xheader))];
Y=[data.H.(char(Yheader))];

axes(handles.axMain);
plot(X,Y,'k.');
xlabel(char(Xheader))
ylabel(char(Yheader))
zoom on
