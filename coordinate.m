function varargout = coordinate(varargin)
% COORDINATE MATLAB code for coordinate.fig
%      COORDINATE, by itself, creates a new COORDINATE or raises the existing
%      singleton*.
%
%      H = COORDINATE returns the handle to a new COORDINATE or the handle to
%      the existing singleton*.
%
%      COORDINATE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COORDINATE.M with the given input arguments.
%
%      COORDINATE('Property','Value',...) creates a new COORDINATE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before coordinate_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to coordinate_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help coordinate

% Last Modified by GUIDE v2.5 03-Sep-2018 02:20:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @coordinate_OpeningFcn, ...
                   'gui_OutputFcn',  @coordinate_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before coordinate is made visible.
function coordinate_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to coordinate (see VARARGIN)

serialPorts = instrhwinfo('serial');
nPorts = length(serialPorts.SerialPorts);
set(handles.listbox_com, 'String', ...
    [{'Select a port'} ; serialPorts.SerialPorts ]);
set(handles.listbox_com, 'Value', 2);   
set(handles.history_box, 'String', cell(1));

% Choose default command line output for coordinate
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes coordinate wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = coordinate_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(hObject,'String'),'Connect') % currently disconnected
    serPortn = get(handles.listbox_com, 'Value');
    if serPortn == 1
        errordlg('Select valid COM port');
    else
        serList = get(handles.listbox_com,'String');
        serPort = serList{serPortn};
        serConn = serial(serPort, 'TimeOut', 1, ...
            'BaudRate', str2num(get(handles.Txt_baudrate, 'String')));
        serConn.InputBufferSize = 10000; % set input buffer size 
        
        try
            fopen(serConn);
            handles.serConn = serConn;
            
            % enable Tx text field and Rx button
            set(handles.Tx_send, 'Enable', 'On');
            set(handles.Rx_button, 'Enable', 'On');
            
            set(hObject, 'String','Disconnect')
        catch e
            errordlg(e.message);
        end
        
    end
else
    set(handles.Tx_send, 'Enable', 'Off');
    set(handles.Rx_button, 'Enable', 'Off');
    
    set(hObject, 'String','Connect')
    fclose(handles.serConn);
end
guidata(hObject, handles);

% --- Executes on selection change in listbox_com.
function listbox_com_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_com (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_com contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_com


% --- Executes during object creation, after setting all properties.
function listbox_com_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_com (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Txt_baudrate_Callback(hObject, eventdata, handles)
% hObject    handle to Txt_baudrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Txt_baudrate as text
%        str2double(get(hObject,'String')) returns contents of Txt_baudrate as a double


% --- Executes during object creation, after setting all properties.
function Txt_baudrate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Txt_baudrate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an1_x_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an1_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an1_x as text
%        str2double(get(hObject,'String')) returns contents of edit_an1_x as a double


% --- Executes during object creation, after setting all properties.
function edit_an1_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an1_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an1_y_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an1_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an1_y as text
%        str2double(get(hObject,'String')) returns contents of edit_an1_y as a double


% --- Executes during object creation, after setting all properties.
function edit_an1_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an1_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an1_z_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an1_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an1_z as text
%        str2double(get(hObject,'String')) returns contents of edit_an1_z as a double


% --- Executes during object creation, after setting all properties.
function edit_an1_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an1_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an2_x_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an2_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an2_x as text
%        str2double(get(hObject,'String')) returns contents of edit_an2_x as a double


% --- Executes during object creation, after setting all properties.
function edit_an2_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an2_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an2_y_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an2_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an2_y as text
%        str2double(get(hObject,'String')) returns contents of edit_an2_y as a double


% --- Executes during object creation, after setting all properties.
function edit_an2_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an2_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an2_z_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an2_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an2_z as text
%        str2double(get(hObject,'String')) returns contents of edit_an2_z as a double


% --- Executes during object creation, after setting all properties.
function edit_an2_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an2_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an3_x_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an3_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an3_x as text
%        str2double(get(hObject,'String')) returns contents of edit_an3_x as a double


% --- Executes during object creation, after setting all properties.
function edit_an3_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an3_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an3_y_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an3_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an3_y as text
%        str2double(get(hObject,'String')) returns contents of edit_an3_y as a double


% --- Executes during object creation, after setting all properties.
function edit_an3_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an3_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an3_z_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an3_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an3_z as text
%        str2double(get(hObject,'String')) returns contents of edit_an3_z as a double


% --- Executes during object creation, after setting all properties.
function edit_an3_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an3_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an4_x_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an4_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an4_x as text
%        str2double(get(hObject,'String')) returns contents of edit_an4_x as a double


% --- Executes during object creation, after setting all properties.
function edit_an4_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an4_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an4_y_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an4_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an4_y as text
%        str2double(get(hObject,'String')) returns contents of edit_an4_y as a double


% --- Executes during object creation, after setting all properties.
function edit_an4_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an4_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an4_z_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an4_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an4_z as text
%        str2double(get(hObject,'String')) returns contents of edit_an4_z as a double


% --- Executes during object creation, after setting all properties.
function edit_an4_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an4_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an5_x_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an5_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an5_x as text
%        str2double(get(hObject,'String')) returns contents of edit_an5_x as a double


% --- Executes during object creation, after setting all properties.
function edit_an5_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an5_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an5_y_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an5_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an5_y as text
%        str2double(get(hObject,'String')) returns contents of edit_an5_y as a double


% --- Executes during object creation, after setting all properties.
function edit_an5_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an5_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an5_z_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an5_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an5_z as text
%        str2double(get(hObject,'String')) returns contents of edit_an5_z as a double


% --- Executes during object creation, after setting all properties.
function edit_an5_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an5_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an6_x_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an6_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an6_x as text
%        str2double(get(hObject,'String')) returns contents of edit_an6_x as a double


% --- Executes during object creation, after setting all properties.
function edit_an6_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an6_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an6_y_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an6_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an6_y as text
%        str2double(get(hObject,'String')) returns contents of edit_an6_y as a double


% --- Executes during object creation, after setting all properties.
function edit_an6_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an6_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_an6_z_Callback(hObject, eventdata, handles)
% hObject    handle to edit_an6_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_an6_z as text
%        str2double(get(hObject,'String')) returns contents of edit_an6_z as a double


% --- Executes during object creation, after setting all properties.
function edit_an6_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_an6_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_setup_anchor.
function pushbutton_setup_anchor_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setup_anchor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in history_box.
function history_box_Callback(hObject, eventdata, handles)
% hObject    handle to history_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns history_box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from history_box

% --- Executes during object creation, after setting all properties.
function history_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to history_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tx_send_Callback(hObject, eventdata, handles)
% hObject    handle to Tx_send (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tx_send as text
%        str2double(get(hObject,'String')) returns contents of Tx_send as a double
TxText = get(handles.Tx_send, 'String');
fprintf(handles.serConn, TxText );

currList = get(handles.history_box, 'String');

set(handles.history_box, 'String', ...
    [currList ; ['Send: ' TxText ] ]);
set(handles.history_box, 'Value', length(currList) + 1 );

set(hObject, 'String', '');

% --- Executes during object creation, after setting all properties.
function Tx_send_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tx_send (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Rx_button.
function Rx_button_Callback(hObject, eventdata, handles)
% hObject    handle to Rx_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try 
    RxText = fscanf(handles.serConn);
    currList = get(handles.history_box, 'String');
    if length(RxText) < 1
        RxText = 'Timeout.';
        set(handles.history_box, 'String', ...
            [currList ; RxText ]);
    else
        set(handles.history_box, 'String', ...
            [currList ; ['Received: ' RxText ] ]);
    end
    set(handles.history_box, 'Value', length(currList) + 1 );
catch e
    disp(e)
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'serConn')
    fclose(handles.serConn);
end
% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in Rcv_button.
function Rcv_button_Callback(hObject, eventdata, handles)
% hObject    handle to Rcv_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for k = 1:str2num(get(handles.sample_cnt, 'String'))
    try 
        RxText = fscanf(handles.serConn);
        currList = get(handles.history_box, 'String');
        if length(RxText) < 1
            RxText = 'Timeout.';
            set(handles.history_box, 'String', ...
                [currList ; RxText ]);
        else
            set(handles.history_box, 'String', ...
                [currList ; ['Received: ' RxText ] ]);
        end
        set(handles.history_box, 'Value', length(currList) + 1 );
    catch e
        disp(e)
    end    
end


function sample_cnt_Callback(hObject, eventdata, handles)
% hObject    handle to sample_cnt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sample_cnt as text
%        str2double(get(hObject,'String')) returns contents of sample_cnt as a double


% --- Executes during object creation, after setting all properties.
function sample_cnt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample_cnt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
