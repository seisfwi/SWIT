function varargout = GUIEditSegyTraceHeader(varargin)
% GUIEditSegyTraceHeader Application M-file for GUIEditSegyTraceHeader.fig
%    FIG = GUIEDITSEGYTRACEHEADER launch GUIEditSegyTraceHeader GUI.
%    GUIEDITSEGYTRACEHEADER('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 05-Jun-2002 09:24:38

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
       data.Horig=varargin{1};
       data.Hname=fieldnames(data.H(1));
 
       data.NShowHeaders=7; % NUMBER OF HEADER VALUES TO LIST
       data.ShowNTraces=10; % NUMBER OF TRACES TO LIST 
       data.StartTrace=1;
       
       set(handles.popHV1,'String',data.Hname);
       set(handles.popHV2,'String',data.Hname);
       set(handles.popHV3,'String',data.Hname);
       set(handles.popHV4,'String',data.Hname);
       set(handles.popHV5,'String',data.Hname);
       set(handles.popHV6,'String',data.Hname);
       set(handles.popHV7,'String',data.Hname);
       
       set(handles.popHV1,'Value',2);
       set(handles.popHV2,'Value',4);
       set(handles.popHV3,'Value',5);
       set(handles.popHV4,'Value',7);
       set(handles.popHV5,'Value',13);
       set(handles.popHV6,'Value',23);
       set(handles.popHV7,'Value',24);

       guidata(fig,data);
       
       set(fig,'HandleVisibility','On')
       GUIEditSegyTraceHeader('actionCreateHandles',fig,handles)
       set(fig,'HandleVisibility','CallBack')

    end
    
	% Wait for callbacks to run and window to be dismissed:
	uiwait(fig);

	if nargout > 0
	  data=guidata(fig);
	  if isstruct(data.H)
	    varargout{1}=data.H;
	  else
	    varargout{1} = fig;
	  end
    end

    try
        close(handles.figure1);
        set(handles.figure1,'Visible','off')
    catch
        SegymatVerbose(sprintf('%s: Could not close WINDOW',mfilename))
    end
	
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
function actionCreateHandles(h,handles)
data=guidata(h);
ShowNTraces=data.ShowNTraces;

MainPos=get(handles.fMain,'Position');
for ih=1:data.NShowHeaders
   poph=findobj('tag',['popHV',num2str(ih)]);
   HandlePos=get(poph,'Position');
   X=HandlePos(1)+HandlePos(3);
   Xdist=(MainPos(3)+MainPos(1))-(HandlePos(1)+HandlePos(3));
   dx=Xdist/ShowNTraces;
   ddx=.1;
   % CREATE PLOT BUTTONS
   Xwidth=HandlePos(1)-MainPos(1);
   TagString=['PH',num2str(ih+10)];
   h = uicontrol('units','normalized','Style', 'pushbutton','String','P',...
        'Position', [MainPos(1)+ddx*Xwidth HandlePos(2) Xwidth-2*ddx*Xwidth HandlePos(4)],...
        'Tag',TagString,...
        'ToolTipstring',['Plot these header values'],...
        'HorizontalAlignment','center',...
        'Callback','GUIEditSegyTraceHeader(''actionPlotHeader'',gcbo,guidata(gcbo))');
   
   % CREATE EDIT BOXES
   for it=1:ShowNTraces;
      Xpos=X+(it-1)*Xdist/(ShowNTraces); 
      
      TagString=['T',num2str(it+10),'H',num2str(ih+10)];
      h = uicontrol('units','normalized','Style', 'edit','String','P',...
        'Position', [Xpos+ddx*dx HandlePos(2) dx-2*ddx*dx HandlePos(4)],...
        'Tag',TagString,...
        'HorizontalAlignment','Right',...
        'Callback','GUIEditSegyTraceHeader(''actionChangeValueAction'',gcbo,guidata(gcbo))');
    
  end
end
for it=1:ShowNTraces;
      Xpos=X+(it-1)*Xdist/(ShowNTraces); 
      HeaderValuePos=get(handles.fHeaderValue,'Position');
      TagString=['T',num2str(it+10)];
      h = uicontrol('units','normalized','Style', 'edit','String','',...
        'Position', [Xpos+ddx*dx HeaderValuePos(2) dx-2*ddx*dx HeaderValuePos(4)],...
        'Tag',TagString,...
        'HorizontalAlignment','Right',...
        'Callback','GUIEditSegyTraceHeader(''actionHandleAction'',gcbo,guidata(gcbo))');
    
end
actionUpdateHandles(h,handles)

% --------------------------------------------------------------------
function actionChangeValueAction(h,handles)
data=guidata(h);


Tag=get(h,'Tag');
trace=str2num(Tag(2:3))-10+(data.StartTrace-1);
ih=str2num(Tag(5:6))-10;
poph=findobj('tag',['popHV',num2str(ih)]);
headername=char(data.Hname(get(poph,'value')));


% ONLY USE THE TYPES IN VALUE IF IT IS A NUMBER
origvalue=getfield(data.H(trace),headername);
value=str2num(get(h,'string'));
if isempty(value), value=origvalue; end

% UPDATE THE HEADER VALUE
data.H(trace)=setfield(data.H(trace),headername,value);

guidata(h,data);
actionUpdateHandles(h,handles)


% --------------------------------------------------------------------
function actionUpdateHandles(h,handles)
data=guidata(h);

% MAKE SURE WE ARE NOT OUT OF BOUNDS
if data.StartTrace<1, data.StartTrace=1; end
if data.StartTrace>(length(data.H)-data.ShowNTraces+1), data.StartTrace=length(data.H)-data.ShowNTraces+1;end

for ih=1:data.NShowHeaders
  poph=findobj('tag',['popHV',num2str(ih)]);
  hname=char(data.Hname(get(poph,'value')));
  for it=[1:data.ShowNTraces];
    hvalue=getfield(data.H(it+data.StartTrace-1),hname);
    TagString=['T',num2str(it+10),'H',num2str(ih+10)];
    set(findobj('Tag',TagString),'String',hvalue);
  end
end
for it=1:data.ShowNTraces;
  TagString=['T',num2str(it+10)];
  set(findobj('Tag',TagString),'String',it+data.StartTrace-1);    
end
guidata(h,data);

% --------------------------------------------------------------------
function actionPlotHeader(h,handles,varargin)
if length(varargin)==1
  headers=varargin{1};
else
  Tag=get(h,'Tag');
  headers=str2num(Tag(3:length(Tag)))-10;
end
data=guidata(h);

traces=[1:1:length(data.H)];

figure
set(gcf,'name','Plot of Trace Header Values') 
nheaders=length(headers);
ch=0;
for ih=headers
 
  ch=ch+1;
  subplot(nheaders,1,ch);

  
  poph=findobj('tag',['popHV',num2str(ih)]);
  hname=char(data.Hname(get(poph,'value')));

  % GET HEADER VALUES FROM STRUCTURE
  for it=1:length(traces), hv(it)=getfield(data.H(it),hname);  end
  
  %bar(traces,hv);
  area(traces,hv,min(hv(:)));
  %hold on;bar(traces,hv);  hold off
  if ch==nheaders, 
    xlabel(['TraceNumber']);
  else
    set(gca,'XtickLabel','');
  end
  ylabel(hname);
  grid on
  zoom on
  
  
end


% --------------------------------------------------------------------
function varargout = pbPlotAll_Callback(h, eventdata, handles, varargin)
data=guidata(h);
actionPlotHeader(h,handles,[1:1:data.NShowHeaders]);


% --------------------------------------------------------------------
function varargout = popHV_Callback(h, eventdata, handles, varargin)
actionUpdateHandles(h,handles)


% --------------------------------------------------------------------
function varargout = popHV1_Callback(h, eventdata, handles, varargin)
data=guidata(h);
data.StartTrace=data.StartTrace+1;
guidata(h,data);
actionUpdateHandles(h,handles)


% --------------------------------------------------------------------
function varargout = popHV2_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = popHV3_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = popHV4_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = popHV5_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = popHV6_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = popHV7_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
% -------------- JUMP IN THE TRACES ----------------------------------
% --------------------------------------------------------------------
%


% --------------------------------------------------------------------
function varargout = pbLeftStart_Callback(h, eventdata, handles, varargin)
data=guidata(h);
data.StartTrace=1;
guidata(h,data);
actionUpdateHandles(h,handles)


% --------------------------------------------------------------------
function varargout = pbLeftJump_Callback(h, eventdata, handles, varargin)
data=guidata(h);
data.StartTrace=data.StartTrace-data.ShowNTraces;
guidata(h,data);
actionUpdateHandles(h,handles)

% --------------------------------------------------------------------
function varargout = pbLeft_Callback(h, eventdata, handles, varargin)
data=guidata(h);
data.StartTrace=data.StartTrace-1;
guidata(h,data);
actionUpdateHandles(h,handles)


% --------------------------------------------------------------------
function varargout = pbRight_Callback(h, eventdata, handles, varargin)
data=guidata(h);
data.StartTrace=data.StartTrace+1;
guidata(h,data);
actionUpdateHandles(h,handles)



% --------------------------------------------------------------------
function varargout = pbRightJump_Callback(h, eventdata, handles, varargin)
data=guidata(h);
data.StartTrace=data.StartTrace+data.ShowNTraces;
StartTrace=data.StartTrace;
guidata(h,data);
actionUpdateHandles(h,handles)


% --------------------------------------------------------------------
function varargout = pbRightEnd_Callback(h, eventdata, handles, varargin)
data=guidata(h);
data.StartTrace=length(data.H)-data.ShowNTraces+1;
guidata(h,data);
actionUpdateHandles(h,handles)



% --------------------------------------------------------------------
function varargout = figure1_ResizeFcn(h, eventdata, handles, varargin)








% --------------------------------------------------------------------
function varargout = pbDone_Callback(h, eventdata, handles, varargin)
uiresume;



% --------------------------------------------------------------------
function varargout = pbClose_Callback(h, eventdata, handles, varargin)
data=guidata(h);
data.H=data.Horig;
guidata(h,data);
uiresume;






% --------------------------------------------------------------------
function varargout = pbReset_Callback(h, eventdata, handles, varargin)
data=guidata(h);
data.H=data.Horig;
guidata(h,data);
actionUpdateHandles(h,handles)
