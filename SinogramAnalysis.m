function varargout = SinogramAnalysis(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SinogramAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @SinogramAnalysis_OutputFcn, ...
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

% --- Executes just before untitled is made visible.
function SinogramAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for untitled
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = SinogramAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% --- Executes on button press in Obzor.
function Obzor_Callback(hObject, eventdata, handles)
% hObject    handle to Obzor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

v = get(handles.checkbox4,'Value');
if v == 1
  
[FileName,PathName] = uigetfile('*.*' ,'Pick a file','MultiSelect','on');
x = FileName;
x1 = char(FileName(1,1))
x2 = char(FileName(1,2));
file12 = strcat(PathName,x1);
file22 = strcat(PathName,x2);
set(handles.text16, 'String', x1); 
x0 = 0;
handles.name = x1
guidata(gcbo, handles);   

f = waitbar(0,'Please, wait...');
pause(.5)
                            
M = importdata(file12);   %Data from Detectors
D2 = importdata(file22);

%Reference
l = get(handles.checkbox3, 'Value');
    str1=get(handles.popupmenu4,'String');
    str2=str1(get(handles.popupmenu4,'Value'));
    J1 = char(str2);
    J1=str2double(J1);
str3= ['J = ', num2str(J1)];
set(handles.text22, 'String',str3);
path = cd;
if J1 == 2.5
    Mref = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw1_ExitDetector.csv']);
    Dref2 = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw1_vars.csv']);
end
    
if J1 == 5.02
    Mref = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw2.1_ExitDetector.csv']);
    Dref2 = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw2.1_vars.csv']);
end
    
if J1 == 1.05
    Mref = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw0.35_ExitDetector.csv']);
    Dref2 = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw0.35_vars.csv']);
end

waitbar(.25,f,'Please, wait...');
pause(1)
                            
%% Data from detectors

Mraw = M.data (1:(size(M.data,1)),30:559);
Mrefraw = Mref.data (:,30:559);  %without some chanels

% Normalization by reference value
Dref = Dref2.data(:,3);
D = D2.data(:,3);
SR = mean(D);   %The average dose value of the monitor chamber
SRref = mean(Dref); %The average reference dose value of the monitor chamber
B = 5168.4; %Background value of the signal

for i = 1:size(Mraw, 1)
    for k = 1:size(Mraw, 2)
        Mnorm(i,k) = (Mraw(i,k)-B)*SRref/SR/Mrefraw(2900,k);
    end
end

waitbar(.50,f,'Please, wait...');
pause(1)
                            
% Projecting onto the X-axis
a1 = 0.064;  %Length of one detector cell
for a = 1:size(Mraw, 2)
    X(a) = 132*sin(a*a1/99)/(33/99+cos(a*a1/99));
end

% Interpolate by 800 points
d1 = 0.66;    %interpolation step
for p = 1:size(Mnorm, 1)
    for a = 1:size(Mnorm, 2)
        i = 0;
        if a ~= size(Mnorm, 2)
        k = Mnorm (p,a+1)-Mnorm(p,a);
        b = Mnorm(p,a) - k*a;
            for x1 = 0.66:d1:576
            i = i + 1;
                if x1 >= a & x1 <= a+1
                Y(p, i)=k.*x1+b;
                end
            end
        else
        end
    end
end

%Renaming matrices
for p = 1:size(Y,1)
    for i = 1:801
    	A12(p,i) = Y(p,i);
    end
end
            
waitbar(.75,f,'Please, wait...');
pause(1)

for p = 1:size(A12, 1)
    for i = 1:size(A12, 2) 
        F11 (p, i) = A12 (p, size(A12, 2) + 1 - i);
    end
end

%Creating a graf data from exit detectors
Y2 = 1:size(F11, 1);
d3 = 0.05;
X3= 0:d3:40;
Y1 = 1:size(F11, 2);   
axes(handles.axes1);
contourf(Y2,X3,F11', 'LineStyle', 'none');
caxis([0 1]);
colorbar(handles.axes1);
colormap jet;
grid on; 
xlabel('Projection number');
ylabel ('X, cm');
title('Data from Exit Detectors');
else

%% Downloading the DICOM RT file
[FileName,PathName] = uigetfile('*.*' ,'Pick a file');
file11 = strcat(PathName,FileName);
handles.file11 = file11;
guidata(gcbo, handles);   
set(handles.text16, 'String', FileName); 
handles.name = FileName;
f = waitbar(0,'Please, wait...');
pause(.5)
D=char(FileName);
D1(1:3)=D(1,(size(D,2)-2):size(D,2));
    
if isequal(char('bin'),D1)
    N = fopen(file11);
    E = fread(N);
    k=0;
    x0 = 0;
    for i = 1:4:((size(E))-3)
        k=k+1;
        T1(k) = E(i,1)+E(i+1,1)+E(i+2,1)+E(i+3,1);
    end
    T=T1'
    for i = 1:(size(T,1)-1)
        if T(i,1) == 191
            T(i,1)= 1;
        end
        if T(i,1) == 190
            T(i,1)= 0.25;
        end
        if T(i,1) == 63
            T(i,1)= 0.5;
        end
        if T(i,1) == 127
            T(i,1)= 0.75;
        end
    end
    for i = 1:(size(T,1)/64)
        Tmlc(i,1:64)=T((1+64*(i-1)):(64+64*(i-1)));
    end
    J1 = 5.02;
    fclose('all')
else
    if isequal(char('txt'),D1)
    N = importdata(file11);
    Tmlc = N;
    else
    A = dicominfo(file11);
    handles.dicom = A;
    guidata(gcbo, handles);
    %The folder for saving the file
    file = A.PatientID; 
    Fi1 = mfilename('fullpath');
    Fi = [cd,'\', file '.xlsx'];
    i = 1;
    i_str=num2str(i);

    %Sinogram export
    S = A.BeamSequence.Item_1.NumberOfControlPoints;
    a = 0;
    for i=1:(S-1)
        B = ['A.BeamSequence.Item_1.ControlPointSequence.Item_' num2str(i) '.Private_300d_10a7'];
        C = eval(B);
        numOfRows = size(C, 1);
        if  numOfRows > a
            a = numOfRows;
        end
    end

    waitbar(.2,f,'Please, wait...');
    pause(1)

    v1 = A.BeamSequence.Item_1.Private_300d_1080;
    v2 = v1';
    v3 = char (v2);
    V = str2double(v3);
    
    t1 = A.BeamSequence.Item_1.Private_300d_1040;
    t2 = t1';
    t3 = char (t2);
    T = str2double(t3);
    
    p = A.BeamSequence.Item_1.NumberOfControlPoints;
    x0 = 0;
    x0 = V*(T*p/51+51*9)/10;

    
    D1=zeros(a,100);
    for i=1:(S-1)
        B = ['A.BeamSequence.Item_1.ControlPointSequence.Item_' num2str(i) '.Private_300d_10a7'];
        C = eval(B);
        numOfRows = size(C, 1);
        if  numOfRows == 0
            D1(1:numOfRows,i)= zeros (1:numOfRows,i);
        else
            D1(1:numOfRows,i)=C(1:numOfRows,1);
        end
    end

    D = D1';
    V = char (D);
    X = cellstr(V);

    %Field size value
    J = A.BeamSequence.Item_1.ControlPointSequence.Item_1.BeamLimitingDevicePositionSequence.Item_2.LeafJawPositions;

    for s = 1 : size(V, 1)
        V1(s,:) = strrep(V (s,:), '.', ',');
    end

    for p = 1 : size(V, 1)
        C = strsplit(V1(p,:),'\');
        W(p, 1:64) = C;
    end


    % Saving the sinogram in a file
    xlswrite(Fi,W);
    fclose('all'); 
    Tmlc = xlsread(Fi);
    delete(Fi); 
    numOfRows = size(Tmlc, 1); 
    numOfCols = size(Tmlc, 2); 
    Tmlc(isnan(Tmlc))=0;
    J = A.BeamSequence.Item_1.ControlPointSequence.Item_1.BeamLimitingDevicePositionSequence.Item_2.LeafJawPositions;
    J1 = round((J(2,1) - J(1,1))/10*100)/100;
    end
end
Tmlc';
waitbar(.4,f,'Please, wait...');
pause(1)

for i = 1:size(Tmlc, 2)
    for p = 1:size(Tmlc, 1)
        if Tmlc(p,i)>1
            Tmlc(p,i)=0;
        end
    end
end
Told = Tmlc;
size(Told)

l = get(handles.checkbox3, 'Value');
if l == 1
    str1=get(handles.popupmenu4,'String');
    str2=str1(get(handles.popupmenu4,'Value'));
    J1 = char(str2);
    J1=str2double(J1);
end
str3= ['J = ', num2str(J1)];
set(handles.text22, 'String',str3);
if J1 == 5.02
%Accounting for the opening time of the leaf itself and two adjacent ones
for p = 1:size(Tmlc, 1)
    for i = 2:(size(Tmlc, 2)-1)
        b(p,i) = 0.69526*Tmlc(p,i);
        k(p,i)=0.10269+0.06569*Tmlc(p,i);
        N(p,i) = (Tmlc(p,i-1)+Tmlc(p,i+1))/2;
        T(p,i) = k(p,i)*N(p,i)+b(p,i);
    end
end

%Accounting for the number of open adjacent leaves
for p = 1:size(Tmlc, 1)
    for i = 2:(size(Tmlc, 2)-1)
        h1=0;
        h2=0;
        Tr(p,i) = 0;
        for k = 1:(i-1)
            h1 = h1+1;
            Tr(p,i) = Tr(p,i)+Tmlc(p,i-k);
             if Tmlc(p,i-k)==0
                break
            end
        end
        for k = (i+1):(size(Tmlc, 2)-1)
             h1 = h1+1;
            Tr(p,i) = Tr(p,i)+Tmlc(p,k);
            if Tmlc(p,k)==0
                break
            end
        end  
        hh(p,i) = (h1-2)/2;
        if hh(p,i) == 0 || hh(p,i) == -0.5 || hh(p,i) == -1
            Rr(p,i) = 0;
        else
        Trs(p,i) = Tr(p,i)*100/h1;
        Aa(p,i) = -0.0021*Trs(p,i)+0.00102;
        kk = 2.2781075;
        x00(p,i) = 0.0014*Trs(p,i)+0.00013;
        Rr(p,i) = Aa(p,i)*exp(-hh(p,i)/kk)+x00(p,i);
        end
    end
end
for p = 1:size(Tmlc, 1)
    for i = 2:(size(Tmlc, 2)-1)
        if Tmlc(p,i)~=0 
            Tmlc(p,i) = T(p,i)+Rr(p,i);
        end
        if Tmlc(p,i)>1
            Tmlc(p,i) = 1;
        end
    end
end
end

if J1 == 2.5
%Accounting for the opening time of the leaf itself and two adjacent ones
for p = 1:size(Tmlc, 1)
    for i = 2:(size(Tmlc, 2)-1)
        b(p,i) = 0.76177*Tmlc(p,i);
        k(p,i)=0.09+0.06566*Tmlc(p,i);
        N(p,i) = (Tmlc(p,i-1)+Tmlc(p,i+1))/2;
        T(p,i) = k(p,i)*N(p,i)+b(p,i);
    end
end

%Accounting for the number of open adjacent leaves
for p = 1:size(Tmlc, 1)
    for i = 2:(size(Tmlc, 2)-1)
        h1=0;
        h2=0;
        Tr(p,i) = 0;
        for k = 1:(i-1)
            h1 = h1+1;
            Tr(p,i) = Tr(p,i)+Tmlc(p,i-k);
            if Tmlc(p,i-k)==0
                break
            end
        end
        for k = (i+1):(size(Tmlc, 2)-1)
             h1 = h1+1;
            Tr(p,i) = Tr(p,i)+Tmlc(p,k);
            if Tmlc(p,k)==0
                break
            end
        end  
        hh(p,i) = (h1-2)/2;
        if hh(p,i) == 0 || hh(p,i) == -0.5 || hh(p,i) == -1
            Rr(p,i) = 0;
        else
        Trs(p,i) = Tr(p,i)*100/h1;
        Aa(p,i) = -0.0013*Trs(p,i)-0.00023;
        kk(p,i) = 0.0009*Trs(p,i)+0.00021;
        x00 = 2.14;
        Rr(p,i) = Aa(p,i)*exp(-hh(p,i)/x00)+kk(p,i);
        end
    end
end
for p = 1:size(Tmlc, 1)
    for i = 2:(size(Tmlc, 2)-1)
        if Tmlc(p,i)~=0 
            Tmlc(p,i) = T(p,i)+Rr(p,i);
        end
        if Tmlc(p,i)>1
            Tmlc(p,i) = 1;
        end
    end
end
end

if J1 == 1.05
%Accounting for the opening time of the leaf itself and two adjacent ones
for p = 1:size(Tmlc, 1)
    for i = 2:(size(Tmlc, 2)-1)
        b(p,i) = 0.80311*Tmlc(p,i);
        k(p,i)=0.07203+0.06823*Tmlc(p,i);
        N(p,i) = (Tmlc(p,i-1)+Tmlc(p,i+1))/2;
        T(p,i) = k(p,i)*N(p,i)+b(p,i);
    end
end

%Accounting for the number of open adjacent leaves
for p = 1:size(Tmlc, 1)
    for i = 2:(size(Tmlc, 2)-1)
        h1=0;
        h2=0;
        Tr(p,i) = 0;
        for k = 1:(i-1)
            h1 = h1+1;
            Tr(p,i) = Tr(p,i)+Tmlc(p,i-k);
            if Tmlc(p,i-k)==0
                break
            end
        end
        for k = (i+1):(size(Tmlc, 2)-1)
             h1 = h1+1;
            Tr(p,i) = Tr(p,i)+Tmlc(p,k);
            if Tmlc(p,k)==0
                break
            end
        end  
        hh(p,i) = (h1-2)/2;
        if hh(p,i) == 0 || hh(p,i) == -0.5 || hh(p,i) == -1
            Rr(p,i) = 0;
        else
        Trs(p,i) = Tr(p,i)*100/h1;
        Aa(p,i) = -0.00085*Trs(p,i)+0.0009;
        kk(p,i) = 0.00056*Trs(p,i)+0.0019;
        x00 = 2.09;
        Rr(p,i) = Aa(p,i)*exp(-hh(p,i)/x00)+kk(p,i);
        end
    end
end

for p = 1:size(Tmlc, 1)
    for i = 2:(size(Tmlc, 2)-1)
        if Tmlc(p,i)~=0 
            Tmlc(p,i) = T(p,i)+Rr(p,i);
        end
        if Tmlc(p,i)>1
            Tmlc(p,i) = 1;
        end
    end
end
end
Told;
dif = zeros(size(Tmlc, 1),size(Tmlc, 2));
for p = 1:size(Tmlc, 1)
    for i = 1:(size(Tmlc, 2))
        if Tmlc(p,i)==0
            dif(p,i)=0;
        else
                dif(p,i) = Told(p,i)/Tmlc(p,i);
        end
    end
end
handles.dif = dif;
guidata(gcbo, handles);

size(Tmlc);
size(Tmlc);
b = 0.35;    %penumbra width
L = 0.625;  %leaf length in the isocenter
F = zeros(size(Tmlc, 1),size(Tmlc, 2));
i = 0;
d2 = 0.05;  %Interpolation step

%Interpolate by X
for p = 1 : size(Tmlc, 1)
    for n = 1 : size(Tmlc, 2)
        for x2 = 0:d2:40
            x3 = x2 - 20;
            i = i+1;
            if x3 < -20+n*L+b/2 & x3 >= -20+n*L-b/2
                F(p, i) = Tmlc(p, n);
            end
            if x3 <= -20+(n+1)*L-b/2 & x3 >= -20+n*L+b/2
               if n ~= 64
                    F(p, i) = (Tmlc(p, n+1)*(x3+20-n*L-b/2) - Tmlc(p, n)*(x3+20-n*L-L+b/2))/(L-b);
               else
                    F(p, i) = 0;
               end
            end
        end
        i = 0;
    end
end

waitbar(.6,f,'Please, wait...');
pause(1)

%Field size accounting
path = cd;
if J1 == 2.5
    Mref = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw1_ExitDetector.csv']);
    Dref2 = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw1_vars.csv']);
end
    
if J1 == 5.02
    Mref = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw2.1_ExitDetector.csv']);
    Dref2 = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw2.1_vars.csv']);
end
    
if J1 == 1.05
    Mref = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw0.35_ExitDetector.csv']);
    Dref2 = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw0.35_vars.csv']);
end
    
handles.data2 = Mref;
guidata(gcbo, handles);
handles.data3 = Dref2;
guidata(gcbo, handles);
Mrefraw = Mref.data (:,30:559); 
B1 = 5168.4;
for i = 1:size(Mrefraw, 1)
    for k = 1:size(Mrefraw, 2)
        Mr2(i,k) = (Mrefraw(i,k)-B1)/Mrefraw(2900,k);
    end
end

%Interpolate by 800 points
d1 = 0.66;    %interpolation step
for p = 1:size(Mr2, 1)
    for a = 1:size(Mr2, 2)
        i = 0;
        if a ~= size(Mr2, 2)
            k = Mr2 (p,a+1)-Mr2(p,a);
            b = Mr2(p,a) - k*a;
            for x1 = 0.66:d1:576
                i = i + 1;
                if x1 >= a & x1 <= a+1
                    Mr1(p, i)=k.*x1+b;
                end
            end
        else
        end
     end
end

if J1 == 2.5
    G = 0.232;
end
if J1 == 5.02
    G = 0.373;
end
if J1 == 1.05
    G = 0.186;
end
    
x3 = -1 :0.2:1;
Klor = G./(pi*x3.^2+pi*G.^2/4); 
MaxKlor= max(Klor);
Knorm = Klor/MaxKlor;
Doze = F;

waitbar(.8, f,'Please, wait...');
pause(1)

F1 = F;

% Accounting for dispersion at the ends of the field
F11 = F1;
% On the left
F21=zeros(size(F1,1),size(F1,2));
for p = 1:size(F1,1)
    for i = 10:(size(F1,2))
        if i>=130
            if F1(p,i)>F1(p,i-1) && F1(p,i-1)==0  && F1(p,i-5)==0
                for n2 = 1:128
                    if F1(p,i-n2)==0 
                        M1 = Mr1(620,353);
                        M2 = F1(p,i+7);
                        if F1(p,i-n2-1)==0
                            F21(p,i-n2+1+7)= F21(p,i-n2+1+7) + Mr1(620,353-n2+1)*M2/M1;
                            F11(p,i-n2+1+7)= F21(p,i-n2+1+7);
                        else
                            for n4=1:14
                                F21(p,i-n2-n4+1+7+1)= F21(p,i-n2-n4+1+7+1) + Mr1(620,353-n2-n4+1+1)*M2/M1;
                                F11(p,i-n2-n4+1+7+1)= F21(p,i-n2-n4+1+7+1);
                            end
                        end
                     end
                end
            end
        else
           if F1(p,i)>F1(p,i-1) && F1(p,i-1)==0  && F1(p,i-5)==0
                for n2 = 1:(i-2)
                    if F1(p,i-n2)==0 
                        M1 = Mr1(620,351);
                        M2 = F1(p,i+7);
                        if F1(p,i-n2-1)==0
                            F21(p,i-n2+1+7)= F21(p,i-n2+1+7) + Mr1(620,351-n2+1)*M2/M1;
                            F11(p,i-n2+1+7)= F21(p,i-n2+1+7);
                        else
                            for n4=1:14
                                F21(p,i-n2-n4+1+7+1)= F21(p,i-n2-n4+1+7+1) + Mr1(620,351-n2-n4+1+1)*M2/M1;
                                F11(p,i-n2-n4+1+7+1)= F21(p,i-n2-n4+1+7+1);
                            end
                        end
                     end
                end
           end 
        end
    end
end

% On the right
F22=zeros(size(F1,1),size(F1,2));
for p = 1:size(F1,1)
    for i = 1:(size(F1,2)-10)
        if i<=(size(F1,2)-130)
            if F1(p,i)>F1(p,i+1) && F1(p,i+1)==0 && F1(p,i+5)==0
                for n2 = 1:129
                    if F1(p,i+n2)==0
                        M1 = Mr1(620,453);
                        M2 = F1(p,i-7);
                        F22(p,i+n2-1-6)=Mr1(620,453+n2-1)*M2/M1;
                        if F21(p,i+n2-1-6)==0
                            F11(p,i+n2-1-6) = F22(p,i+n2-1-6);
                        else
                            F11(p,i+n2-1-6) = F22(p,i+n2-1-6)+F21(p,i+n2-1-6);
                        end
                    end
                end
            end
        else
         if F1(p,i)>F1(p,i+1) && F1(p,i+1)==0 && F1(p,i+5)==0
                for n2 = 1:(size(F1,2)-i)
                    if F1(p,i+n2)==0
                        M1 = Mr1(620,453);
                        M2 = F1(p,i-7);
                        F22(p,i+n2-1-6)=Mr1(620,453+n2-1)*M2/M1;
                        if F21(p,i+n2-1-6)==0
                            F11(p,i+n2-1-6) = F22(p,i+n2-1-6);
                        else
                            F11(p,i+n2-1-6) = F22(p,i+n2-1-6)+F21(p,i+n2-1-6);
                        end
                    end
                end
         end 
        end
    end
end


 
for p = 1:size(F11,1)
    for i = 1:(size(F11,2))
        if F11(p,i)<0
            F11(p,i)=0;
        end
        if F11(p,i)>1
            F11(p,i)=1;
        end
    end
end

for i = 1:64
    n11 = i*8.32505+19.62351;
    n22 = round(n11);
    n3=n22/576*800;
    n3=round(n3);
    for p = 1:size(F11,1)
        if Tmlc(p,i)==0
            dif2(p,i)=0;
        else
            dif2(p,i) = F11(p,n3)/Tmlc(p,i);
        end
    end
end
handles.dif2 = dif2;
guidata(gcbo, handles); 
waitbar(1,f,'Please, wait...');
pause(1)
%Create a graf for sinogram
X3= 0:d2:40;
Y2 = 1:size(F11, 1);

if l == 0
handles.x0 = x0;
guidata(gcbo, handles);
end

axes(handles.axes1);
size(Told)
ynew = 1:64;
xnew = 1:size(Told, 1);

handles.Told = Told;
guidata(gcbo, handles);
contourf(Y2,X3,F11', 'LineStyle', 'none');
caxis([0 1]);
colorbar(handles.axes1);
datacursormode on

grid on; 
colormap jet;
xlabel('Projection number');
ylabel ('X, cm');
title('Data from Planning System');
end 


%Saving F1
handles.data1 = F11;
guidata(gcbo, handles); 

close(f);


% --- Executes on button press in obzor2.
function obzor2_Callback(hObject, eventdata, handles)
[FileName,PathName] = uigetfile('*.*' ,'Pick a file','MultiSelect','on');
x = FileName;
x1 = char(FileName(1,1));
x2 = char(FileName(1,2));
file12 = strcat(PathName,x1);
file22 = strcat(PathName,x2);
set(handles.text17, 'String', x1); 

f = waitbar(0,'Пожалуйста, подождите...');
pause(.5)
                            
M = importdata(file12);   %data from detectors
D2 = importdata(file22);

% Reference
v = get(handles.checkbox4,'Value');
if v == 1 
l = get(handles.checkbox3, 'Value');

    str1=get(handles.popupmenu4,'String');
    str2=str1(get(handles.popupmenu4,'Value'));
    J1 = char(str2);
    J1=str2double(J1);

str3= ['J = ', num2str(J1)];
set(handles.text22, 'String',str3);    
path = cd;
if J1 == 2.5
    Mref = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw1_ExitDetector.csv']);
    Dref2 = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw1_vars.csv']);
end
    
if J1 == 5.02
    Mref = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw2.1_ExitDetector.csv']);
    Dref2 = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw2.1_vars.csv']);
end
    
if J1 == 1.05
    Mref = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw0.35_ExitDetector.csv']);
    Dref2 = importdata([path, '\Data from Detectors\Reference\detData_sin_jaw0.35_vars.csv']);
end
else
  Mref = handles.data2; 
  Dref2 = handles.data3; 
end
waitbar(.25,f,'Please, wait...');
pause(1)
                            
%% Data from detectors
F1 = handles.data1; 

nom = size(M.data,1)-size(F1,1)+1
Mraw = M.data (1:(size(M.data,1)),30:559);
Mrefraw = Mref.data (:,30:559);

% Normalization by reference value
Dref = Dref2.data(:,3);
D = D2.data(:,3);
SR = mean(D);   %The average dose value of the monitor chamber
SRref = mean(Dref); %The average reference dose value of the monitor chamber
B = 5168.4; %Background value of the signal

for i = 1:size(Mraw, 1)
    for k = 1:size(Mraw, 2)
        Mnorm(i,k) = (Mraw(i,k)-Mrefraw(2800,k))*SRref/SR/Mrefraw(2900,k);
    end
end

waitbar(.50,f,'Please, wait...');
pause(1)
                            
% Projecting onto the X-axis
a1 = 0.064;  %Length of one detector cell
for a = 1:size(Mraw, 2)
    X(a) = 132*sin(a*a1/99)/(33/99+cos(a*a1/99));
end

% Interpolate by 800 points
d1 = 0.66;    %interpolation step
for p = 1:size(Mnorm, 1)
    for a = 1:size(Mnorm, 2)
        i = 0;
        if a ~= size(Mnorm, 2)
        k = Mnorm (p,a+1)-Mnorm(p,a);
        b = Mnorm(p,a) - k*a;
            for x1 = 0.66:d1:576
            i = i + 1;
                if x1 >= a & x1 <= a+1
                Y(p, i)=k.*x1+b;
                end
            end
        else
        end
    end
end

%Checking for matching projections
for a = 1:size(Y, 2)
    l(a) = Y(size(Y, 1),a)-0.06;
end
l1 = max(l(:))
size(Y,1)
size(F1, 1)
size(Y,2)
if size(F1, 1)>size(Y,1)
n=0;
        for i = 1:50
            if Y(i,286)<0.003
                Y(i,286)
            n = n+1 
            else
                break
            end
        end
        
        for i = n:(n+4)
        t(i) = abs(Y(size(Y,1)-i,286)-F1((size(F1,1)-size(Y,1)),286));
        end
        [t,n1]=min(t)
        n2=n+n1
        Y22(1:(size(Y, 1)-n2+1),:)= Y((n2):((size(Y,1))),:);
        Y22((size(Y, 1)-n2+2):(size(F1,1)),:) = 0;   
else
    if l1>0
        Y22(1:size(F1, 1),:)=Y((nom):(size(Y,1)),:);
    else
        if nom > 1
        Y22(1:size(F1, 1),:)=Y((nom-1):(size(Y,1)-1),:);
        else
        Y22(1:size(F1, 1),:)= Y(1:(size(Y,1)),:);   
        end
    end
end
%Renaming matrices
for p = 1:size(F1,1)
    for i = 1:size(F1,2)
    	A12(p,i) = Y22(p,i);
    	A2(p,i) = F1(p,i);
    end
end
            
waitbar(.75,f,'Please, wait...');
pause(1)

for p = 1:size(A12, 1)
    for i = 1:size(A12, 2) 
        A1 (p, i) = A12 (p, size(A12, 2) + 1 - i);
    end
end

waitbar(1,f,'Пожалуйста, подождите...');
pause(1)
close(f); 
%%Convert to sinogram
v = get(handles.checkbox4,'Value');
v2 = get(handles.checkbox3,'Value');
if v == 0
    if v2 == 0
A = handles.dicom;
dif = handles.dif;
Told = handles.Told;

size(A1);
size(dif);
A1(217,209);
        d2 = 0.05;
for p = 1:size(A1,1)
for i = 1:64
    n11 = i*8.32505+19.62351;
    n22 = round(n11);

    n3=i*800/64;
    t1=round(n3);
    n12 = (i+1)*8.32505+19.62351;
    n23 = round(n12);

    n4=(i+1)*800/64;
    t2=round(n4);
    k=0;
    f=0;
    f1=0;
    t3 = t2-t1;
    tmin = round((t1-t3));
    tmax = round((t1+t3));
    n=tmax-tmin;
    Tr=ones(n);
    for t = tmin:tmax
        k=k+1;
        if t>=1 & t<=(size(A1,2)-1)
        Tr(k) = abs(A1(p,t) - F1(p,t1));
        end
    end
    [y,f]= min(Tr(:));
    f1 = round((t1-t3))+f;
    p;
    i;
    size(A1);
    size(dif);
    det(p,i) = A1(p,f1)*dif(p,i);
    if i == 16 & p == 217
        F1(217,t1);
        f1;
        A1(217,f1);
        dif(217,17);
        det(217,17);
    end
end
end
   
%Creating a graf data from exit detectors
x0 = handles.x0;
x0;
Y2 = 0:x0/(size(A1, 1)-1):x0;
d3 = 0.05;
X3= 0:d3:40;
Y1 = 1:size(A1, 1);   
axes(handles.axes2);
size(A1)
size(F1)
size(A12)
Y22 = 1:size(A1, 1);
contourf(Y22,X3,A1', 'LineStyle', 'none');
size(det);
ynew = 1:64;
xnew = 1:size(det,1);
det';
handles.det = det;
guidata(gcbo, handles);
caxis([0 1]);
colorbar(handles.axes2);
colormap jet;
grid on; 
xlabel('Projection number');
ylabel ('X, cm');
title('Data from exit detectors');
    end
end

%Creating a graf data from exit detectors
Y2 = 1:size(A1, 1);
d3 = 0.05;
X3= 0:d3:40;
Y1 = 1:size(A1, 2);   
axes(handles.axes2);
contourf(Y2,X3,A1', 'LineStyle', 'none');
caxis([0 1]);
colorbar(handles.axes2);
colormap jet;
grid on; 
xlabel('Projection number');
ylabel ('X, cm');
title('Data from exit detectors');



%Saving A1, A2
handles.data1 = A1;
guidata(gcbo, handles);
handles.data2 = A2;
guidata(gcbo, handles);

%1D gamma-index for leaves
v = get(handles.checkbox4,'Value');
if v == 1 
    DTA = 3;
    dosed = 0.03;
else
        DTA = 5;
        dosed = 0.05;
end

for i = 1:64
    n11 = i*8.32505+19.62351;
    n22 = round(n11);
    n3=n22/576*800;
 n3=round(n3);
    for p1 = 1:size(A1,1)
        for p2 = 1:size(A1,1)
            r2 = (p1-p2)^2;
            d2 = (A1(p1,n3)-A2(p2,n3))^2;
            Ga2(p2,i)= r2/(DTA^2)+d2/dosed^2; 

        end
        G1(p1,i)=min(Ga2(:,i));
        G(p1,i) = sqrt(G1(p1,i));
    end
    Gm(i) = mean(G(:,i));
    Gm(i) = roundn(Gm(i),-2);
    St(i) = std(G(:,i));
    St(i) = roundn(St(i),-2);
    n1(i) = 0;
    n2(i) = 0;
    for p = 1:size(A1,1)
        if A1(p,n3) > 0.1
        n2(i) = n2(i) + 1;
        if G(p,i)<1
            n1(i) = n1(i) + 1;
        end
        end
    end
    gamma1(i) = n1(i)/n2(i)*100;
    gamma(i) = roundn(gamma1(i),-1);
    if gamma(i) == 0
    gamma(i) = 100;
end
end
    Gm(i) = mean(G(:,i));
    Gm(i) = roundn(Gm(i),-2);
    St(i) = std(G(:,i));
    St(i) = roundn(St(i),-2);
handles.Gm = Gm;
guidata(gcbo, handles);
handles.St = St;
guidata(gcbo, handles);
handles.gamma = gamma;
guidata(gcbo, handles);


axes(handles.axes4);
x4 = 1:64;
size(gamma);
plot(x4, gamma, 'Marker', '.', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'black');
ymin = min(gamma);
ymin1 = ymin - 10;
grid on;
axis([1 64 0 100]);
mylegend  = ['1D gamma-index for leaves'];
xlabel('Leaf number');
ylabel('Gamma-index, %'); 
title(mylegend);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_Callback(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A11 = handles.data1;  
A22 = handles.data2;

%Creating a slider for different projections
set(hObject, 'Max', size(A11,1), 'Min', 1);
a = get(hObject,'Value');
p = round(a);
myString = ['p = ', num2str(p)];
handles.data3 = p;
guidata(gcbo, handles);

v = get(handles.checkbox4,'Value');

%Construction of two datasets on the selected projection
axes(handles.axes6);
size(A11,1);
size(A11,2);
Y6 = A11(p,1:801);
Y7 = A22(p,1:801);

X6 = 1:801;
v2 = get(handles.checkbox3,'Value');
if v2 == 0
    if v == 0 
Told = handles.Told;
det = handles.det;
det';
size(Told);
size(det);
Xnew = 1:64;
Ynew1 = Told(p,1:64);
Ynew2 = det(p,1:64);
plot(Xnew,Ynew1*100,Xnew,Ynew2*100, 'LineWidth', 0.8); 
axis([0 64 0 100]);
else
axes(handles.axes6);
size(A11,1);
size(A11,2);
Y6 = A11(p,1:801);
Y7 = A22(p,1:801);
% X6 = 0:0.05:40;
X6 = 1:801;
plot(X6,Y7*100,X6,Y6*100, 'LineWidth', 0.8); 
axis([0 800 0 100]);
    end
else
    axes(handles.axes6);
size(A11,1);
size(A11,2);
Y6 = A11(p,1:801);
Y7 = A22(p,1:801);
% X6 = 0:0.05:40;
X6 = 1:801;
  plot(X6,Y7*100,X6,Y6*100, 'LineWidth', 0.8); 
axis([0 800 0 100]);  
end

grid on
v = get(handles.checkbox4,'Value');
if v == 1
xlabel('X, cm');
ylabel ('Leaf open time, %');
legend('Data from ED 1','Data from ED 2')  ; 
else
xlabel('X, cm');
ylabel('Leaf open time, %'); 
legend('Data from PS','Data from ED')  ;
end
set(handles.text, 'String', myString);       
 
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties.
function slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Отображение числа проекции
p1 = str2double(get(hObject,'String'));
handles.datap1 = p1;
guidata(gcbo, handles);

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x1 = str2double(get(hObject,'String'));
handles.datax1 = x1;
guidata(gcbo, handles);
% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit6_Callback(hObject, eventdata, handles)
p2 = str2double(get(hObject,'String'));
handles.datap2 = p2;
guidata(gcbo, handles);
function edit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit7_Callback(hObject, eventdata, handles)
x2 = str2double(get(hObject,'String'));
handles.datax2 = x2;
guidata(gcbo, handles);

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% 2D Gamma-Analysis
function pushbutton4_Callback(hObject, eventdata, handles)

%Data on the range of calculation of the gamma index
p1 = handles.datap1;
p2 = handles.datap2;
x1 = handles.datax1;
x2 = handles.datax2;
A11 = handles.data1;  
A22 = handles.data2;

%Calculation of Gamma
DTA1 = get(handles.edit8, 'String');
DTA = str2num(DTA1);
dosed1 = get(handles.edit9, 'String')
dosed = str2num(dosed1)
d1 = 0.055;
x11 = round(x1/d1);
x22 = round(x2/d1);

A1 = A11(p1:p2,x11:x22);
A2 = A22(p1:p2,x11:x22);
size1 = size(A2);
size2 = size(A1);

%2D gamma
if size1 == size2  
    for i = 1:size(A1,1)
        for j = 1:size(A1,2)
            for k = 1:size(A1,1)
                for l = 1:size(A1,2)
                    r2 = (i-k)^2+(j-l)^2;
                    d2 = (A1(i,j)-A2(k,l))^2;
                    Ga(k,l) = r2/(DTA.^2)+d2/dosed.^2;
                end
            end
            G(i,j)=min(min(Ga));
        end
    end
    G = sqrt(G);   
else
    fprintf = ('Size of the datasets do not match' /n');
end

for i = 1:size(A1,1)
        for j = 1:size(A1,2)
            if (A1(i,j)-A2(i,j))<0
                G(i,j)=-G(i,j);
            end
        end
end
% Colormap for difference
axes(handles.axes3);
Z = G';
x21 = x11*0.055;
x22 = x22*0.055;
Y6 = p1:p2;
X6 = x21:0.055:x22;
contourf(Y6,X6,Z, 'LineStyle', 'none', 'Levelstep',0.1);
caxis([-1.5 1.5]);
colormap jet;
colorbar(handles.axes3);
grid on; 
xlabel('Projection number') 
ylabel ('X, cm')
title('Difference, %')

% Histogram for gamma
axes(handles.axes5);
nbins = 100;
Zmax = max(max(Z));
xbins = [0 Zmax/2];
histogram(Z(:), nbins,'Normalization','probability');
xlim([0 Zmax/2]);
grid on
xlabel('Gamma') 
ylabel ('Probability, %')
title('Histogram of the gamma')

% Gamma calculation
n1 = 0;
n2 = 0;
n3 = 0;
for p = 1:size(A1,2)
    for i = 1:size(A1,1)
        if A1(i,p) > 0.05
            n2 = n2 + 1;
            if abs(Z(p,i))<1
                n1 = n1 + 1;
            end
        end
    end
end

for p = 1:(size(A11,1))
    for i = 1:(size(A11,2))
        A33(p,i) = abs(A11(p,i)-A22(p,i));
    end
end
min(min(A22));
max(max(A11));
D = mean(mean(A33));
myString5 = ['Mean error = ', num2str(D)];
set(handles.text21, 'String', myString5);
gamma1 = n1/n2*100;
gamma = roundn(gamma1,-1);
myString2 = ['Number of points that passed the criteria = ', num2str(gamma), '%'];  
myString3 = ['Number of all points = ', num2str(n2)];    
set(handles.text2, 'String', myString2);           
set(handles.text20, 'String', myString3);           

function text16_ButtonDownFcn(hObject, eventdata, handles)

% Slider for analyzing the movement of leaf
function slider4_Callback(hObject, eventdata, handles)
A11 = handles.data1;  
A22 = handles.data2;

% Slider for analyzing different projection
set(hObject, 'Max', 64, 'Min', 1);
a = get(hObject,'Value');
n = round(a)
myString = ['n = ', num2str(n)];
n1 = n*8.32505+19.62351;
n2 = round(n1)

% Construction of two distributions on the selected projection
n3=n2/577*800;
n3=round(n3);
size(A22,2);
axes(handles.axes12);
Y7 = A22(1:size(A22,1),n3);
X6 = 1:size(A22,1);
Y8 = A11(1:size(A22,1),n3);
X8 = 1:size(A22,1);

v = get(handles.checkbox4,'Value');
v2 = get(handles.checkbox3,'Value');
if v2 == 0
    if v == 0 
        Told = handles.Told;
        det = handles.det;
        det';
        size(Told);
        size(det);
        Xnew = 1:64;
        n3;
        Ynew1 = Told(1:size(A22,1),n);
        Ynew2 = det(1:size(A22,1),n);
        plot(X6,Ynew1*100,X8,Ynew2*100);
        axis([1 200 0 110]);
    else
        plot(X6,Y7*100,X8,Y8*100); 
        axis([1 200 0 110]);
    end
else
    plot(X6,Y7*100,X8,Y8*100); 
    axis([1 200 0 110]);
end

grid on
v = get(handles.checkbox4,'Value');
if v == 1
legend('Data from ED 1','Data from ED 2')   
else
legend('Data from PS','Data from ED')  
end 
xlabel('Projection number') 
ylabel ('Leaf open time, %')
set(handles.text18, 'String', myString);  
Gm = handles.Gm;  
St = handles.St
gamma = handles.gamma;  
myString2 = ['Gamma = ', num2str(Gm(n)), ' ± ', num2str(St(n)), ', gamma = ', num2str(gamma(n)), '%'];
set(handles.text19, 'String', myString2);

% Slider for analyzing the movement of each petal
function slider4_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% Error for distance
function edit8_Callback(hObject, eventdata, handles)

function edit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Error for the dose
function edit9_Callback(hObject, eventdata, handles)

function edit9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Open file *bin
function checkbox3_Callback(hObject, eventdata, handles)
v = get(hObject,'Value');
if v == 1
    set(handles.popupmenu4, 'Enable', 'on');
else
    set(handles.popupmenu4, 'Enable', 'off');
end

% Choosing the field size
function popupmenu4_Callback(hObject, eventdata, handles)
function popupmenu4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% 1D gamma
function pushbutton5_Callback(hObject, eventdata, handles)
% Data on the range of calculation of the gamma index
v = get(handles.checkbox4,'Value');
v2 = get(handles.checkbox3,'Value');
if v2 == 0
if v == 0 
x0 = handles.x0; 
A11 = handles.data1;  
A22 = handles.data2;
size(A11)
size(A22)
p11 = handles.datap1;
p12 = handles.datap1*size(A11,1)/x0;
p1 = ceil(p12);
p21 = handles.datap2;
p22 = handles.datap2*size(A11,1)/x0;
p2 = round(p22);
x1 = handles.datax1;
x2 = handles.datax2;
else
A11 = handles.data1;  
A22 = handles.data2;
size(A11)
size(A22)
p11 = handles.datap1;
p1 = ceil(p11);
p21 = handles.datap2;
p2 = round(p21);
x1 = handles.datax1;
x2 = handles.datax2;

end
else
    A11 = handles.data1;  
A22 = handles.data2;
size(A11)
size(A22)
p11 = handles.datap1;
p1 = ceil(p11);
p21 = handles.datap2;
p2 = round(p21);
x1 = handles.datax1;
x2 = handles.datax2;
end
% Gamma calculation
DTA1 = get(handles.edit8, 'String');
DTA = str2num(DTA1);
dosed1 = get(handles.edit9, 'String');
dosed = str2num(dosed1);
d1 = 0.055;
x11 = round(x1/d1);
x22 = round(x2/d1);
A1 = A11(p11:p21,x11:x22);
A2 = A22(p11:p21,x11:x22);
size(A1)
% 1D gamma
size1 = size(A2);
size2 = size(A1);
if size1 == size2  
    for p = 1:size(A1,1)
        for j = 1:size(A1,2)
            for k = 1:size(A1,2)
                    r2 = ((k-j)*0.05*10)^2;
                    d2 = (A1(p,j)-A2(p,k))^2;
                    Ga(p,k) = r2/(DTA.^2)+d2/dosed.^2;
            end
            G(p,j)=min(Ga(p,:));
            G2(p,j) = sqrt(G(p,j));
            if (A1(p,j)-A2(p,j))<0
                G2(p,j) = - G2(p,j);
            end
        end
    end
else
    fprintf = ('Размеры матриц не совпадают/n');
end

% Colormap for difference
axes(handles.axes3);
Z = G2';
Zmax = max(max(Z));
size(Z);
for p = 1 : (size(Z,2)/2)
    for i = 1 : fix(size(Z,1)*0.11)
    v=fix(i/0.11);
    w=p*2;
    Z1(i,p)=Z(v,w);
end
end
x21 = x11*0.055;
x22 = x22*0.055;
v = get(handles.checkbox4,'Value');
if v2 == 0
    if v == 0 
    Y6 = p11:x0/(size(A11,1)-1):p21;
    Y66 = p11:p21;
    size(Z);
    Y7 = p1:2:(p2-2);
    X6 = x21:0.055:x22;
    X7 = x21:0.5:(x22);
else
    Y66 = p11:p21;
    X6 = x21:0.055:x22;
    end
else
  Y66 = p11:p21;
    X6 = x21:0.055:x22;
    end
contourf(Y66,X6,Z, 'LineStyle', 'none', 'Levelstep',0.1);
caxis([-2 2]);
colormap jet;
colorbar(handles.axes3);
grid on; 
xlabel('Projection number'); 
ylabel ('X, cm');
title('Difference, %');

% Gamma calculation
n1 = 0;
n2 = 0;
n3 = 0;
for p = 1:size(A1,2)
    for i = 1:size(A1,1)
        if A1(i,p) > 0.05
            n2 = n2 + 1;
            if abs(Z(p,i))<1
                n3 = n3 + 1;
            end
        end
    end
end

for p = 1:(size(A1,1))
    for i = 1:(size(A1,2))
        A33(p,i) = (A1(p,i)-A2(p,i));
    end
end
min(min(A22));
max(max(A11));
D = mean(mean(A33));
myString5 = ['Mean error = ', num2str(D)];
set(handles.text21, 'String', myString5);

% Histogram for gamma
axes(handles.axes5);
nbins = 100;
Zmax = max(max(Z));
xbins = [0 Zmax/2];
histogram(Z(:), nbins,'Normalization','probability');
xlim([0 Zmax/2]);
grid on
xlabel('Gamma-index') 
ylabel ('Probability, %')
title('Histogram of the gamma')

gamma2 = n3/n2*100;
gamma3 = roundn(gamma2,-1);
myString2 = ['Number of points that passed the criteria = ', num2str(gamma3), '%'];  
myString3 = ['Number of all points = ', num2str(n2)];    
set(handles.text2, 'String', myString2);           
set(handles.text20, 'String', myString3);  


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v=get(hObject,'Value')
if v==1
WW1 = 'Load data from detectors 1';
WW2 = 'Load data from detectors 2';
  set(handles.text16, 'String', WW1); 
  set(handles.text17, 'String', WW2); 
  set(handles.popupmenu4, 'Enable', 'on');
else
WW1 = 'Load data from planning system';
WW2 = 'Load data from detectors';
  set(handles.text16, 'String', WW1); 
  set(handles.text17, 'String', WW2);
  set(handles.popupmenu4, 'Enable', 'off');
end
% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
A = handles.dicom;
dif = handles.dif2;
A1 = handles.data1;
name = handles.name;
size(A1);
size(dif);
det = handles.det;

for p = 1:(size(det,1))
    for i = 1:64
        N2{p,i} =num2str(det(p,i));
        N3{i} =num2str(det(p,i));
    end
    T2=strjoin(N3,'\\');
    T3{p}=T2;
end
X2 = T3';
V2 = char(X2);
D2 = uint8(V2);

for i=1:(size(det,1))
        B = ['A.BeamSequence.Item_1.ControlPointSequence.Item_' num2str(i) '.Private_300d_10a7']; 
        C = eval(B);
        numOfRows = size(C, 1);

            Cnew = D2(i,:)';
            eval (['A.BeamSequence.Item_1.ControlPointSequence.Item_' num2str(i) '.Private_300d_10a7 = Cnew']);
end
name2 = ['new' name];
dicomwrite ([], name2, A, 'CreateMode', 'copy');


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
