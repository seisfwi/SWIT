function varargout = GUIEditTextualFileHeader(varargin)
% GUIEditTextualFileHeader : GUI for editing the SGY header
%
% Call : 
%   SH=GUIEditTextualFileHeader(SH);
%
%
% Example : 
%   [Data,STH,SH]=ReadSegy('841_m.sgy');
%   SH=GUIEditTextualFileHeader(SH);
%
%

% GUIEDITTEXTUALFILEHEADER M-file for GUIEditTextualFileHeader.fig
%      GUIEDITTEXTUALFILEHEADER, by itself, creates a new GUIEDITTEXTUALFILEHEADER or raises the existing
%      singleton*.
%
%      H = GUIEDITTEXTUALFILEHEADER returns the handle to a new GUIEDITTEXTUALFILEHEADER or the handle to
%      the existing singleton*.
%
%      GUIEDITTEXTUALFILEHEADER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIEDITTEXTUALFILEHEADER.M with the given input arguments.
%
%      GUIEDITTEXTUALFILEHEADER('Property','Value',...) creates a new GUIEDITTEXTUALFILEHEADER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIEditTextualFileHeader_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIEditTextualFileHeader_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIEditTextualFileHeader

% Last Modified by GUIDE v2.5 05-Jan-2009 17:27:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIEditTextualFileHeader_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIEditTextualFileHeader_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);


if nargin & isstr(varargin{1})
  gui_State.gui_Callback = str2func(varargin{1});
end



if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end


% End initialization code - DO NOT EDIT

% --- Executes just before GUIEditTextualFileHeader is made visible.
function GUIEditTextualFileHeader_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIEditTextualFileHeader (see VARARGIN)

% Choose default command line output for GUIEditTextualFileHeader
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Load TextualFileHeader if given
if length(varargin)>0
  if nargin & isstruct(varargin{1})
    if isfield(varargin{1},'TextualFileHeader')
    data=guidata(gcf);    
    data.SegyHeader=varargin{1};
    data.SegyHeaderOrg=data.SegyHeader;
    
    if isempty(find(data.SegyHeader.TextualFileHeader>200))
        HeaderType=2; % ASCII
    else
        HeaderType=1; % EBCDIC
    end
    set(handles.popType,'Value',HeaderType);


    guidata(gcf, data);
    UpdateText(hObject, handles)  
    end 
  end
end


% UIWAIT makes GUIEditTextualFileHeader wait for user response (see UIRESUME)
uiwait(handles.figTextualHeader);

% --- Outputs from this function are returned to the command line.
function varargout = GUIEditTextualFileHeader_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
data=guidata(gcf);    
varargout{1} =     data.SegyHeader;
close(handles.figTextualHeader);


% --- Executes during object creation, after setting all properties.
function eTextHeader_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eTextHeader (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function eTextHeader_Callback(hObject, eventdata, handles)
% hObject    handle to eTextHeader (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eTextHeader as text
%        str2double(get(hObject,'String')) returns contents of
%        eTextHeader as a double
%eTextHeader_Callback(handles.eTextHeader, eventdata, handles)
  HS=get(hObject,'String');
  VAL=get(hObject,'Value');

  data=guidata(hObject);
  data.HS=HS;
  data.VAL=VAL;
  %disp(sprintf('LINE=''%s''',HS(VAL,:)))
  guidata(hObject, data);
  
  % UPDATE TEXT STRING
  set(handles.eEditLine,'String',HS(VAL,:))

  
  

% --- Executes during object creation, after setting all properties.
function popType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popType.
function popType_Callback(hObject, eventdata, handles)
% hObject    handle to popType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        popType
  UpdateText(hObject,handles)
  
  
  
function UpdateText(hObject, handles)  
  data=guidata(hObject);
  
  HeaderType=get(handles.popType,'Value');

  if HeaderType==1,
    txt=char(ebcdic2ascii(data.SegyHeader.TextualFileHeader));
  else
    txt=char(data.SegyHeader.TextualFileHeader');
  end
  %disp('--')
  %reshape(txt,80,40)'
  %disp('--')
  %disp(size(txt))
  try
    s=[];
    for i=1:40;
      s=[s,txt((i-1)*80+1:(i-1)*80+80)];
      try
          if i~=40, s=[s,'|']; end
      catch
          keyboard
      end
    end
  catch
    if (HeaderType==1)
      warndlg(['The Textual Header is most probably NOT EBCDIC formatted'],'Textual Header Info')
    else
      warndlg(['The Textual Header is most probably NOT ASCII formatted'],'Textual Header Info')
    end
  end
  
  try
    set(handles.eTextHeader,'String',s)
  catch
    disp('WHUUPS')
  end
  

  
  


% --- Executes during object creation, after setting all properties.
function eEditLine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eEditLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function eEditLine_Callback(hObject, eventdata, handles)
% hObject    handle to eEditLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eEditLine as text
%        str2double(get(hObject,'String')) returns contents of eEditLine as a double

data=guidata(hObject);
if ~isfield(data,'HS')
    eTextHeader_Callback(handles.eTextHeader, eventdata, handles)
    data=guidata(hObject);
end

try
   
    NEWSTR=get(hObject,'String');
    n=length(NEWSTR);
    if n<80
        NEWSTR(n+1:80)=' ';
    end
    
    n=length(NEWSTR);
    if n>=80
        data.HS(data.VAL,:)=NEWSTR(1:80);
    else
        data.HS(data.VAL,1:n)=NEWSTR(1:n);
    end
    
        
    HStransp=data.HS';
   
    
    if get(handles.popType,'Value')==1;
        % EBCDIC
        data.SegyHeader.TextualFileHeader=double(ascii2ebcdic(HStransp(:)));
    else
        % ASCII
        data.SegyHeader.TextualFileHeader=double(HStransp(:));
    end
    guidata(gcf, data);
    
    UpdateText(hObject, handles)

    eTextHeader_Callback(handles.eTextHeader, eventdata, handles)
catch
    keyboard


end


% --- Executes on button press in pbClose.
function pbClose_Callback(hObject, eventdata, handles)
% hObject    handle to pbClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume;

% --- Executes on button press in pbReset.
function pbReset_Callback(hObject, eventdata, handles)
% hObject    handle to pbReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=guidata(gcf);    
data.SegyHeader=data.SegyHeaderOrg;
guidata(gcf, data);
UpdateText(hObject, handles)



% --- Executes on button press in pbFont.
function pbFont_Callback(hObject, eventdata, handles)
% hObject    handle to pbFont (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
d = uisetfont(handles.eTextHeader);
% Apply those settings to c2
set(handles.eEditLine,d);

