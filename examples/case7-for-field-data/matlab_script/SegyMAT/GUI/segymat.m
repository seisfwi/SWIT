function varargout = segymat(varargin)
% segymat : Garphical User Interface for SegyMAT
%
%   (C) 2001-2004 Thomas Mejer Hansen, tmh@gfy.ku.dk/thomas@cultpenguin.com
%
 

% SEGYMAT M-file for segymat.fig
%      SEGYMAT, by itself, creates a new SEGYMAT or raises the existing
%      singleton*.
%
%      H = SEGYMAT returns the handle to a new SEGYMAT or the handle to
%      the existing singleton*.
%
%      SEGYMAT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEGYMAT.M with the given input arguments.
%
%      SEGYMAT('Property','Value',...) creates a new SEGYMAT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before segymat_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to segymat_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

  
% Edit the above text to modify the response to help segymat

% Last Modified by GUIDE v2.5 05-Jan-2009 17:29:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @segymat_OpeningFcn, ...
                   'gui_OutputFcn',  @segymat_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);

% Next line is edited               
% ORIGINAL : if nargin & isstr(varargin{1})
if nargin>1 & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

if nargin==1 & isstr(varargin{1})
    [path,file,suffix]=fileparts(varargin{1});
    segyfile.pathname=path;
    if length(suffix)>1,
        segyfile.filename=[file,'.',suffix];
    else    
        segyfile.filename=file;
    end
    % FIGURE OUT HOW TO GET THE hObject ??
    % fReadSegy_Callback(hObject, eventdata, handles,segyfile);
end

% End initialization code - DO NOT EDIT


% --- Executes just before segymat is made visible.
function segymat_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to segymat (see VARARGIN)

  SegymatVerbose(['GUI : OPENING GUI'],20)
% Choose default command line output for segymat
handles.output = hObject;

% Initialize PlotPref 
handles.PlotPref.Show=0;

% Update handles structure
guidata(hObject, handles);

% Resize Figure
fMain_ResizeFcn(hObject, eventdata, handles)

%  read in segy file
segyfile.pathname='c:\MaTLAB6p5\work\';
segyfile.filename='test.sgy';
%segyfile.filename='peaksh.segy';
%fReadSegy_Callback(hObject, eventdata, handles,segyfile);

% disable edit menu
UpdateMenus(hObject, eventdata, handles)

% UIWAIT makes segymat wait for user response (see UIRESUME)
% uiwait(handles.fMain);


% --- Outputs from this function are returned to the command line.
function varargout = segymat_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when fMain window is resized.
function fMain_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to fMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  
if (isfield(handles,'fMain')==0)
    return
end

SegymatVerbose(['GUI : Resizing'],20)

data=guidata(hObject);

% Get Size of Main Window
[mainPos]=get(handles.fMain,'Position');

d=4; % Distance to edges
d_ax=.2; % distance between axes
if isfield(data,'PlotPref')==0,
    data.PlotPref.Show=1;
    guidata(hObject,data);
end    
if isfield(data.PlotPref,'Show')==0,
    data.PlotPref.Show=1;
    guidata(hObject,data);
end

% fHeader Location
if data.PlotPref.Show==0,
    fHeaderW=d;
    set(handles.fHeader,'Visible','Off');
    set(handles.text7,'Visible','Off'); 
    set(handles.txtGain,'Visible','Off'); 
    set(handles.pbGainUp,'Visible','Off'); 
    set(handles.pbGainDown,'Visible','Off'); 
    set(handles.eGainMin,'Visible','Off'); 
    set(handles.eGainMax,'Visible','Off'); 
    set(handles.popTop,'Visible','Off'); 
    set(handles.popBot,'Visible','Off'); 
    set(handles.txtTop,'Visible','Off'); 
    set(handles.txtBot,'Visible','Off'); 
    set(handles.txtColormap,'Visible','Off'); 
    set(handles.popColormap,'Visible','Off'); 
    set(handles.txtStyle,'Visible','Off'); 
    set(handles.popStyle,'Visible','Off'); 
else
    fHeaderW=30; % WIDTH OF INFO PAGE
    x0=d;y0=d;
    w=fHeaderW-d;
    h=mainPos(4)-2*d;
    x1=x0+w;
    y1=y0+h;
    
    Wuic=8; % WIDTH OF UICONTROL
    Huic=1.4; %    HEIGHT OF UICONTROL
    
    set(handles.fHeader,'Position',[x0 y0 w h])
    d_frame=1;
    x0=x0+d_frame;
    wframe=w-2*d_frame;
    % TXT HEADER
    ytop=mainPos(4)-d;
%    set(handles.text7,'Position',[d+d ytop fHeaderW-3*d 2]);
    set(handles.text7,'Position',[x0 ytop wframe 2]);

    % GAIN
    ytop=ytop-d_ax-Huic;
    set(handles.txtGain,'Position',[x0,ytop,wframe,Huic])
    ytop=ytop-d_ax-Huic;
    w_sign=4;
    set(handles.eGainMin,'Position',[x0,ytop,wframe-w_sign,Huic])
    set(handles.pbGainDown,'Position',[x0+wframe-w_sign,ytop,w_sign,Huic])
    ytop=ytop-d_ax-Huic;
    set(handles.eGainMax,'Position',[x0,ytop,Wuic*1.5,Huic])
    set(handles.pbGainUp,'Position',[x0+Wuic*1.5,ytop,Wuic/2,Huic])
    set(handles.eGainMax,'Position',[x0,ytop,wframe-w_sign,Huic])
    set(handles.pbGainUp,'Position',[x0+wframe-w_sign,ytop,w_sign,Huic])

    % TopSub
    ytop=ytop-2*d_ax-Huic;
    set(handles.txtTop,'Position',[x0,ytop,wframe,Huic])
    ytop=ytop-d_ax-Huic;
    set(handles.popTop,'Position',[x0,ytop,wframe,Huic])
    
    % BaseSub
    ytop=ytop-2*d_ax-Huic;
    set(handles.txtBot,'Position',[x0,ytop,wframe,Huic])
    ytop=ytop-d_ax-Huic;
    set(handles.popBot,'Position',[x0,ytop,wframe,Huic])

    % Style
    ytop=ytop-2*d_ax-Huic;
    set(handles.txtStyle,'Position',[x0,ytop,wframe,Huic])
    ytop=ytop-d_ax-Huic;
    set(handles.popStyle ,'Position',[x0,ytop,wframe,Huic])
    
    % Colormap
    ytop=ytop-2*d_ax-Huic;
    set(handles.txtColormap,'Position',[x0,ytop,wframe,Huic])
    ytop=ytop-d_ax-Huic;
    set(handles.popColormap ,'Position',[x0,ytop,wframe,Huic])
    
    %
    
    % MAKE ALL VISIBLE
    set(handles.fHeader,'Visible','On');
    set(handles.text7,'Visible','On');
    set(handles.txtGain,'Visible','On'); 
    set(handles.pbGainUp,'Visible','On'); 
    set(handles.pbGainDown,'Visible','On'); 
    set(handles.eGainMin,'Visible','On'); 
    set(handles.eGainMax,'Visible','On'); 
    set(handles.popTop,'Visible','On'); 
    set(handles.popBot,'Visible','On'); 
    set(handles.txtTop,'Visible','On'); 
    set(handles.txtBot,'Visible','On'); 
    set(handles.txtColormap,'Visible','On'); 
    set(handles.popColormap,'Visible','On'); 
    set(handles.txtStyle,'Visible','On'); 
    set(handles.popStyle,'Visible','On'); 
end

% AXIS
axLeft=fHeaderW+2*d;
Hsmallax=4;

% axBot
axBotX=axLeft;
axBotY=d;
axBotW=mainPos(3)-axLeft-2*d;
axBotH=Hsmallax;
axB=get(handles.axBot,'Position');
set(handles.axBot,'Position',[axBotX axBotY axBotW axBotH])

% axTop
axTopX=axLeft;
axTopY=mainPos(4)-Hsmallax-d;
axTopW=axBotW;
axTopH=Hsmallax;
axT=get(handles.axTop,'Position');
set(handles.axTop,'Position',[axTopX axTopY axTopW axTopH])

% axMain
axMainX=axLeft;
axMainY=axBotY+Hsmallax+d_ax;
axMainW=axBotW;
axMainH=mainPos(4)-axMainY-Hsmallax-d-d_ax;
set(handles.axMain,'Position',[axMainX axMainY axMainW axMainH])


% --------------------------------------------------------------------
function UpdatePrefs(hObject, eventdata, handles)
  SegymatVerbose(['GUI : UpdatePrefs'],20)
  data=guidata(hObject);

if isfield(data,'SegyData')==0,
    set(handles.eGainMin,'Enable','off');
    set(handles.eGainMin,'String','');
    set(handles.eGainMax,'Enable','off');
    set(handles.eGainMax,'String','');
   return
else
    set(handles.eGainMin,'Enable','on');
    set(handles.eGainMax,'Enable','on');
end

if isfield(data.PlotPref,'caxis')==0
    data.PlotPref.caxis=[min(data.SegyData(:)) max(data.SegyData(:))];
end 
set(handles.eGainMin,'String',data.PlotPref.caxis(1));
set(handles.eGainMax,'String',data.PlotPref.caxis(2));

if isfield(data.PlotPref,'TraceHeaders')==0,
    data.PlotPref.TraceHeaders=fieldnames(data.SegyTraceHeaders);
end 
if isfield(data.PlotPref,'TopPlot')==0, data.PlotPref.TopPlot=7; end
if isfield(data.PlotPref,'BotPlot')==0, data.PlotPref.BotPlot=13; end
set(handles.popTop,'String',data.PlotPref.TraceHeaders);
set(handles.popTop,'value',data.PlotPref.TopPlot);
set(handles.popTop,'Enable','On');
set(handles.popBot,'String',data.PlotPref.TraceHeaders);
set(handles.popBot,'value',data.PlotPref.BotPlot);
set(handles.popBot,'Enable','On');
guidata(hObject,data);

% --------------------------------------------------------------------
function UpdatePlots(hObject, eventdata, handles)
  SegymatVerbose(['GUI : Update Plots'],20)
  %data=guidata(hObject);

  UpdateMainPlot(hObject, eventdata, handles)
  
  UpdateTopPlot(hObject, eventdata, handles)
  UpdateBotPlot(hObject, eventdata, handles)
  
axes(handles.axMain);%%%zoom on;

% --------------------------------------------------------------------
function UpdateMainPlot(hObject, eventdata, handles)
  SegymatVerbose(['GUI : Update Main Plot'],20)
  data=guidata(hObject);
  Style=get(handles.popStyle,'Value'); % GET PLOTTING STYLE
  axes(handles.axMain);
  
  if ((Style==1)|(Style==4)|(Style==5))
    imagesc(data.SegyTrace,data.SegyTime,data.SegyData);   
    caxis(data.PlotPref.caxis)
    hold on
  end 

  gain=max(abs(data.PlotPref.caxis));
  if ((Style==2)|(Style==4))
    wiggle(data.SegyTrace,data.SegyTime,data.SegyData,'wiggle',gain);
  end
  if ((Style==3)|(Style==5))
    wiggle(data.SegyTrace,data.SegyTime,data.SegyData,'VA',gain);
  end
  hold off
  set(gca,'YaxisLocation','Right');
  set(gca,'XtickLabel',[]);
  ylabel('TWT [ms]')

  %% RESET ZOOM SUCH THAT THIS STATE IS THE 'ZOOM OUT' STATE.
  zoom reset;
  
  %%%zoom on;



% --------------------------------------------------------------------
function UpdateTopPlot(hObject, eventdata, handles)
  SegymatVerbose(['GUI : Update TOP Plots'],60)
  data=guidata(hObject);
  axes(handles.axMain);axMain=axis;
  
  pop=get(handles.popTop,'String');
  ipop=get(handles.popTop,'value');
  hname=pop{ipop};
  
  axes(handles.axTop)
  if isfield(data,'dTopPlot')
		if ipop~=data.dTopPlot_ipop
			data.dTopPlot=bar([data.SegyTraceHeaders.(hname)]);
			data.dTopPlot_ipop=ipop; 
			guidata(hObject,data)				
		end
	else
	  data.dTopPlot=bar([data.SegyTraceHeaders.(hname)]);
		data.dTopPlot_ipop=ipop; 
	  guidata(hObject,data)	
	end
  ax=axis;axis([axMain(1) axMain(2) ax(3) ax(4)])
  set(gca,'XaxisLocation','Top');
  set(gca,'YaxisLocation','Right');
  xlabel(hname)
  axes(handles.axMain);%%%zoom on;
  
  
function UpdateBotPlot(hObject, eventdata, handles)
  SegymatVerbose(['GUI : Update BOT Plots'],60)
  data=guidata(hObject);
  axes(handles.axMain);axMain=axis;
  
  pop=get(handles.popBot,'String');
  ipop=get(handles.popBot,'value');
  hname=pop{ipop};
  
  axes(handles.axBot)
	if isfield(data,'dBotPlot')
		if ipop~=data.dBotPlot_ipop
			data.dBotPlot=plot([data.SegyTraceHeaders.(hname)]);
			data.dBotPlot_ipop=ipop; 
			guidata(hObject,data)				
		end
	else
		data.dBotPlot=plot([data.SegyTraceHeaders.(hname)]);
	  data.dBotPlot_ipop=ipop; 
		guidata(hObject,data)	
	end
 	
  % bar([data.SegyTraceHeaders.(hname)]);
  ax=axis;axis([axMain(1) axMain(2) ax(3) ax(4)])
  set(gca,'XaxisLocation','Bot');
  set(gca,'YaxisLocation','Right');
  xlabel(hname)
  axes(handles.axMain);%%%zoom on;

% --------------------------------------------------------------------
function UpdateGain(hObject, eventdata, handles)
  SegymatVerbose(['GUI : UpdateGain'],20)
  data=guidata(hObject);
  set(handles.eGainMin,'String',data.PlotPref.caxis(1));
  set(handles.eGainMax,'String',data.PlotPref.caxis(2));
  axes(handles.axMain);
  caxis(data.PlotPref.caxis);
  drawnow;
  
% --------------------------------------------------------------------
function KeyPressFcn_Callback(hObject, eventdata, handles)
data=guidata(hObject);
Key=get(gcf,'CurrentCharacter');
%SegymatVerbose(['GUI : KeyPressFcn :  ',num2str(double(Key)),' ',char(Key)],20)
%disp(char(Key));
%disp(double(Key));

%%%%%%%%%%%
% GAIN
%if (Key=='+'),
if ((double(Key)==30)|(double(Key)==43)|(double(Key)==115))
    data.PlotPref.caxis=data.PlotPref.caxis./(1.2);
    guidata(hObject,data);
    UpdateGain(hObject, eventdata, handles) ;
end
%if (Key=='-');
if ((double(Key)==31)|(double(Key)==45)|(double(Key)==120))
    data.PlotPref.caxis=data.PlotPref.caxis.*(1.2);
    guidata(hObject,data);
    UpdateGain(hObject, eventdata, handles) ;
end

%% SHOW PREFS
if (lower(Key)=='h')
    if isfield(data.PlotPref,'Show')==0, data.PlotPref.Show=1;end
    data.PlotPref.Show=1-data.PlotPref.Show;
    guidata(hObject,data);
    fMain_ResizeFcn(hObject,[],handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ZOOM ON MAIN AX
if isfield(data,'zoomfac')==0
    data.zoomfac=1.5;
    guidata(hObject,data);
end 
zoom_in_keys=[29,97]; % ARROW RIGHT OR 'a'
zoom_out_keys=[28,122]; % ARROW LEFT OR 'z'
if isempty(Key)==1
  Key=0;
end

try
if find(double(Key)==zoom_in_keys);
	zoom(data.zoomfac)
	UpdateTopPlot(hObject, eventdata, handles)
	UpdateBotPlot(hObject, eventdata, handles)
end
catch
SegymatVerbose(sprintf('%s : failed to zoom',mfilename))
%keyboard;
end
if find(double(Key)==zoom_out_keys);
    zoom(-data.zoomfac)
	UpdateTopPlot(hObject, eventdata, handles)
	UpdateBotPlot(hObject, eventdata, handles)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MOVE
if isfield(data,'movefac')==0
    % Movefac=1; move in stepos of size of window
    % Movefac=0.5; move in steps of half-size of window
    data.movefac=.02;
    guidata(hObject,data);
end 
if (double(Key)>=49)&(double(Key)<=117), % NUM PAD
	
  go_left_keys=[52,117]; % 2,j
	go_right_keys=[54,111]; % 2,j
	go_down_keys=[50,107]; % 2,j
	go_up_keys=[56]; % 2,j
	
  xlim=get(handles.axMain,'Xlim');
  ylim=get(handles.axMain,'Ylim');
  wx=(xlim(2)-xlim(1)).*data.movefac;
  wy=(ylim(2)-ylim(1)).*data.movefac;

  if (find(go_up_keys==double(Key))), set(handles.axMain,'Ylim',ylim-wy); end % 8
  if (find(go_down_keys==double(Key))), set(handles.axMain,'Ylim',ylim+wy); end % 2
  if (find(go_left_keys==double(Key))), set(handles.axMain,'Xlim',xlim-wx); end % 4 
  if (find(go_right_keys==double(Key))), set(handles.axMain,'Xlim',xlim+wx); end % 6
	
  if (double(Key)==49), % 1
    set(handles.axMain,'Ylim',ylim+wy); 
    set(handles.axMain,'Xlim',xlim-wx);
  end % 1
  if (double(Key)==51),  % 3
    set(handles.axMain,'Ylim',ylim+wy); 
    set(handles.axMain,'Xlim',xlim+wx);
  end % 3
  if (double(Key)==55),  % 7
    set(handles.axMain,'Ylim',ylim-wy); 
    set(handles.axMain,'Xlim',xlim-wx);
  end % 7
  if (double(Key)==57), % 9
    set(handles.axMain,'Ylim',ylim-wy); 
    set(handles.axMain,'Xlim',xlim+wx);
  end % 9
  if (double(Key)==53),  % 5
    axes(handles.axMain)
    zoom out;
    
  end % 5
  xlim_after=get(handles.axMain,'Xlim');
  if (sum(xlim==xlim_after)<2)
    UpdateTopPlot(hObject, eventdata, handles)
    UpdateBotPlot(hObject, eventdata, handles)
  end
end 


%
% IO FUNCTIONS
%

% --------------------------------------------------------------------
function mFileOpen_Callback(hObject, eventdata, handles)
% hObject    handle to mFileOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    data=guidata(hObject);
    if isfield(data,'SegyFile');
        segyfile=DIAOpenSEGY(data.SegyFile);
    else
        segyfile=DIAOpenSEGY;
    end
	fReadSegy_Callback(hObject, eventdata, handles,segyfile);
catch
end
% --------------------------------------------------------------------
function mOpenFileFast_Callback(hObject, eventdata, handles)
% hObject    handle to mOpenFileFast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
	[segyfile.filename,segyfile.pathname]=uigetfile( ...
    {'*.segy;*.SEGY;*.SEGY;*.sgy','All Segy files'; ...
     '*.su;*.SU;*.sU;*.Su','All SU files'; ...
     '*','All Files'},...
    'Pick A SEGY file');
	fReadSegy_Callback(hObject, eventdata, handles,segyfile);
catch
end

% --------------------------------------------------------------------
function fReadSegy_Callback(hObject, eventdata, handles,segyfile);
data=guidata(hObject);
try
	[dpath,dfile,dsuffix]=fileparts(segyfile.filename);
	if (strcmp(lower(dsuffix),'.su'))
		[data.SegyData,data.SegyTraceHeaders,data.SegyHeader]=ReadSu(fullfile(segyfile.pathname,segyfile.filename));	
    else	
        if isfield(segyfile,'JUMP');
            k=0;
            try
            if segyfile.JUMP.enable==1
                k=k+1;v{k}='jump';
                k=k+1;v{k}=segyfile.JUMP.jump;
            end
            if segyfile.TIMERANGE.enable==1
                k=k+1;v{k}='trange';
                k=k+1;v{k}=segyfile.TIMERANGE.min;
                k=k+1;v{k}=segyfile.TIMERANGE.max;
            end
            if segyfile.TRACEHEADER.enable==1
                k=k+1;v{k}='minmax';
                k=k+1;v{k}=segyfile.TRACEHEADER.name;
                k=k+1;v{k}=segyfile.TRACEHEADER.min;
                k=k+1;v{k}=segyfile.TRACEHEADER.max;
            end
            catch
                SegymatVerbose(sprintf('%s : failed to make use of headers for reading',mfilename))
                %keyboard
            end
            if k>0                
                [data.SegyData,data.SegyTraceHeaders,data.SegyHeader]=ReadSegy(fullfile(segyfile.pathname,segyfile.filename),v);
            else
                [data.SegyData,data.SegyTraceHeaders,data.SegyHeader]=ReadSegy(fullfile(segyfile.pathname,segyfile.filename));
            end
        else
            [data.SegyData,data.SegyTraceHeaders,data.SegyHeader]=ReadSegy(fullfile(segyfile.pathname,segyfile.filename));
        end 
    end
    if isempty(data.SegyData)
          f = warndlg('No data was read !', 'ReadSegy warning', 'modal');
    else
        set(handles.fMain,'name',['SegyMAT : ',segyfile.filename])
        data.SegyFile=segyfile;
        data.SegyTime=[1:1:data.SegyHeader.ns].*data.SegyHeader.dt/1e+6;
        data.SegyTrace=[1:length(data.SegyTraceHeaders)];
        data.PlotPref.Show=1;
        data.PlotPref.caxis=[min(data.SegyData(:)) max(data.SegyData(:))];
        guidata(hObject,data);
        fMain_ResizeFcn(hObject, eventdata, handles);
        UpdatePrefs(hObject, eventdata, handles);
        UpdatePlots(hObject, eventdata, handles);
    end
catch
    errordlg('An error occured while reading the SEGY file','Error reading SGY file','modal')
    %keyboard
end 
UpdateMenus(hObject, eventdata, handles);

% --------------------------------------------------------------------
function mFileSave_Callback(hObject, eventdata, handles)
% hObject    handle to mFileSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
	data=guidata(hObject);
	file=fullfile(data.SegyFile.pathname,data.SegyFile.filename);
	ButtonName=questdlg(['Are you sure you want to override :',file], ...
                       'Warning !!!', ...
                       'Yes','No','No');
  if strcmp(ButtonName,'Yes')
		[dpath,dfile,dsuffix]=fileparts(file);
		if (strcmp(lower(dsuffix),'.su'))
			WriteSuStructure(file,data.SegyHeader,data.SegyTraceHeaders,data.SegyData);
		else
			WriteSegyStructure(file,data.SegyHeader,data.SegyTraceHeaders,data.SegyData);
		end
	end
catch
	errordlg('An error occured while writing SEGY file','Error writing SGY file','modal')
end


% --------------------------------------------------------------------
function mFileSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to mFileSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
	data=guidata(hObject);
	[segyfile.filename,segyfile.pathname]=uiputfile( ...
    {'*.segy;*.SEGY;*.SEGY;*.sgy','All Segy files'; ...
     '*.su;*.SU;*.sU;*.Su','All SU files'; ...
     '*','All Files'},...
    'Save as ');
	file=fullfile(segyfile.pathname,segyfile.filename);
	[dpath,dfile,dsuffix]=fileparts(file);
	if (strcmp(lower(dsuffix),'.su'))
		WriteSuStructure(file,data.SegyHeader,data.SegyTraceHeaders,data.SegyData);
	else
		WriteSegyStructure(file,data.SegyHeader,data.SegyTraceHeaders,data.SegyData);
	end
catch
	errordlg('An error occured while writing SEGY file','Error writing SGY file','modal')
end


% --------------------------------------------------------------------
function mHelp_Callback(hObject, eventdata, handles)
% hObject    handle to mHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mHelpHelp_Callback(hObject, eventdata, handles)
% hObject    handle to mHelpHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SegymatHelp('index')


% --------------------------------------------------------------------
function mHelpAbout_Callback(hObject, eventdata, handles)
% hObject    handle to mHelpAbout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DIAAboutSegymat;


% --- Executes on button press in pbGainUp.
function pbGainUp_Callback(hObject, eventdata, handles)
% hObject    handle to pbGainUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  data=guidata(hObject);
  data.PlotPref.caxis=data.PlotPref.caxis./(1.2);
  guidata(hObject,data);
  UpdateGain(hObject, eventdata, handles) ;

	Style=get(handles.popStyle,'Value'); % GET PLOTTING STYLE
  if ((Style==2)|(Style==4|(Style==3)|(Style==5)))
		UpdateMainPlot(hObject, eventdata, handles)	
	end
	

% --- Executes on button press in pbGainDown.
function pbGainDown_Callback(hObject, eventdata, handles)
% hObject    handle to pbGainDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  data=guidata(hObject);
  data.PlotPref.caxis=data.PlotPref.caxis.*(1.2);
  guidata(hObject,data);
  UpdateGain(hObject, eventdata, handles) ;
  
	Style=get(handles.popStyle,'Value'); % GET PLOTTING STYLE
  if ((Style==2)|(Style==4|(Style==3)|(Style==5)))
		UpdateMainPlot(hObject, eventdata, handles)	
	end

% --- Executes during object creation, after setting all properties.
function eGainMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eGainMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function eGainMin_Callback(hObject, eventdata, handles)
% hObject    handle to eGainMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eGainMin as text
%        str2double(get(hObject,'String')) returns contents of eGainMin as a double
data=guidata(hObject);
data.PlotPref.caxis(1)=str2double(get(hObject,'String'));
guidata(hObject,data);
UpdateGain(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function eGainMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eGainMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function eGainMax_Callback(hObject, eventdata, handles)
% hObject    handle to eGainMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eGainMax as text
%        str2double(get(hObject,'String')) returns contents of eGainMax as a double
data=guidata(hObject);
data.PlotPref.caxis(2)=str2double(get(hObject,'String'));
guidata(hObject,data);
UpdateGain(hObject, eventdata, handles);




% --- Executes during object creation, after setting all properties.
function popTop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popTop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popTop.
function popTop_Callback(hObject, eventdata, handles)
% hObject    handle to popTop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popTop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popTop
  SegymatVerbose(['GUI : Update popTOP'],20)
  UpdateTopPlot(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popBot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popBot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popBot.
function popBot_Callback(hObject, eventdata, handles)
% hObject    handle to popBot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popBot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popBot
  SegymatVerbose(['GUI : Update popBOT'],20)
  UpdateBotPlot(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popColormap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popColormap.
function popColormap_Callback(hObject, eventdata, handles)
% hObject    handle to popColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popColormap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popColormap
  SegymatVerbose(['GUI : Update popColormap'],20)
  data=guidata(hObject);
  cmap=get(hObject,'String');
  icmap=get(hObject,'value');
  axes(handles.axMain);
  colormap(cmap{icmap});
  %%% zoom on;
  



% --- Executes during object creation, after setting all properties.
function popStyle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popStyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popStyle.
function popStyle_Callback(hObject, eventdata, handles)
% hObject    handle to popStyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popStyle contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popStyle
  SegymatVerbose(['GUI : popStyle'],20)
  UpdateMainPlot(hObject, eventdata, handles);
  UpdateTopPlot(hObject, eventdata, handles)
  UpdateBotPlot(hObject, eventdata, handles)


%
% MENUS
%
% --------------------------------------------------------------------
function UpdateMenus(hObject, eventdata, handles)
  SegymatVerbose(['GUI : Update Menus'],20)
	data=guidata(hObject);
	
	if isfield(data,'SegyData')==0
		set(handles.mEdit,'Visible','Off')
	else
		set(handles.mEdit,'Visible','On')
	end
	
% --------------------------------------------------------------------
function mFile_Callback(hObject, eventdata, handles)
% hObject    handle to mFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mEdit_Callback(hObject, eventdata, handles)
% hObject    handle to mEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mEditSH_Callback(hObject, eventdata, handles)
% hObject    handle to mEditSH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=guidata(hObject);
try
	data.SegyHeader=GUIEditSegyHeader(data.SegyHeader);
	guidata(hObject,data)
catch
end
% --------------------------------------------------------------------
function mEditSTH_Callback(hObject, eventdata, handles)
% hObject    handle to mEditSTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=guidata(hObject);
try
	data.SegyTraceHeaders=GUIEditSegyTraceHeader(data.SegyTraceHeaders);
	guidata(hObject,data)
	UpdateTopPlot(hObject, eventdata, handles)
	UpdateBotPlot(hObject, eventdata, handles)
catch
  SegymatVerbose('Something Went wrong calling GUISegyMat')
end






% --------------------------------------------------------------------
function mPlotXY_Callback(hObject, eventdata, handles)
% hObject    handle to mPlotXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=guidata(hObject);
try
    GUIPlotXY(data.SegyTraceHeaders);    
catch
  SegymatVerbose('Something Went wrong calling GUISegyMat')
end



% --------------------------------------------------------------------
function mplot_Callback(hObject, eventdata, handles)
% hObject    handle to mplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function mEditTextualHeader_Callback(hObject, eventdata, handles)
% hObject    handle to mEditTextualHeader (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=guidata(hObject);
try
	data.SegyHeader=GUIEditTextualFileHeader(data.SegyHeader);
	guidata(hObject,data)
catch
end

