% MTRN4010, Project 1
% Author: Runkai Zhao, z5146927
% Program: This is a main function, including to implement PartB E and F
% Citation: Some parts and namings of this program is refered from the sample code
% provided by lecturer, like extracting data from raw measurements.

function MyProgram(DataFileName)

clc(); close all; clf;
% In case the caller does not specify the input argument, we propose a default one.
if ~exist('DataFileName','var'), DataFileNameLaser ='Laser__2C.mat'; end
if ~exist('DataFileName','var'), DataFileNameGyroscope ='IMU_dataC'; end
if ~exist('DataFileName','var'), DataFileNameSpeeder ='Speed_dataC'; end
load(DataFileNameLaser); 
load(DataFileNameGyroscope);
load(DataFileNameSpeeder); % return Vel struct

global theta;
global vehicleLocation;
global velocity;
global angularRate;
global flag;
flag.flagPause = 0;
global initialX;
global initialY;


% get the yaw angle and the position of local coordinate from the speed
% encoder and gyroscope
[theta, vehicleLocation] = EstAttitude();
% select the some angualr rate from IMU
[angularRate, indexInterest] = ExtractAngualrRate();

% update the yaw angle and the local position
theta = theta(indexInterest);
vehicleLocation.X = vehicleLocation.X(indexInterest);
vehicleLocation.Y = vehicleLocation.Y(indexInterest);

% update the velocity
velocity = Vel.speeds(indexInterest);

%initialization of EKF
initilizationEKF();

%Create graphical object for refreshing data
figure(1);clf();
layoutHandle = tiledlayout(1,2);
nexttile;
handle1 = plot(0,0,'b.');      % to be used for showing the laser points
axis([-10,10,0,20]);                         % focuses the plot on this region (of interest, close to the robot)
xlabel('X(m)');
ylabel('Y(m)');
title('Global Coordinate Frame');
hold on;
handle2 = plot(0,0,'*');           % create an empty title..
hold on;
handle3 = plot(0,0,'b*');
zoom on;grid on;

%figure(2);
nexttile;
handle4 = plot(0,0, 'b.');

axis([-5,5,0,10]); 
xlabel('X(m)');
ylabel('Y(m)');
s = sprintf('Local Coordinate Frame\n(red circle means to exceed a tolerance of 40cm)');
title(s);
hold on; zoom on; grid on;
handle5 = plot(0,0,'*');  
handle6 = plot(0,0);
handle7 = plot(0,0, 'ro');
set(gcf, 'position', [300 300 1000 400]);
legend('Currently', 'Initially');
s = sprintf('dead reckoning');
title(layoutHandle, s);

% create handles for text targets.
for i = 1:5
    textHandle(i) = text(0, 0, '', 'FontSize', 8, 'Color', 'green');
end 
zoom on; grid on;

uicontrol('Style','pushbutton','String','Pause/Cont.','Position',[10,1,80,20],'Callback',{@MyCallBackA,1});
uicontrol('Style','pushbutton','String','END Now','Position',[90,1,80,20],'Callback',{@MyCallBackA,2});


N = dataL.N;        
global Xe_History;
Xe_History= zeros(2,N);

for i=1:3:N
% i = 2101;
    if (flag.flagPause), pause(0.2); continue; end
    
    scan_i = dataL.Scans(:,i);
    disp(i);
    t =  0.08*3;
    
    [xGlobal, yGlobal, NumCurOOI, MinDistncIndice] = ProcessScan(scan_i, handle1,handle2,handle3, handle4, handle5, handle6, handle7,  i, t);
    
    % cleaner text objects
    for i = 1:5
        set(textHandle(i), 'string', '');
    end 
    
    % label current OOI
    if NumCurOOI ==  0 || NumCurOOI > 5
        disp("Something weird happened: # of OOI is 0 or greater than 5.");
    else
        for counterCurOOI = 1:NumCurOOI  
            set(textHandle(counterCurOOI), 'String', num2str(MinDistncIndice(counterCurOOI)), 'Position', [xGlobal(counterCurOOI)-0.5; yGlobal(counterCurOOI)-0.5; 0], 'FontSize', 8, 'Color', 'green');
        end
    end
    
    pause(0.01) ;                   % wait for ~10ms
end

figure(3);hold on;
plot(vehicleLocation.X, vehicleLocation.Y, 'b');
plot(Xe_History(1,:), Xe_History(2,:), 'r');

disp('Done. Bye.');

return;
end


% Process the data from the separated files
function [xGlobal, yGlobal, NumCurOOI, MinDistncIndice] = ProcessScan(scan, handle1,handle2,handle3, handle4, handle5, handle6, handle7, i, t)
global initialX;
global initialY;
global velocity;
global angularRate;

%global t;
% Extract range and intensity information, from raw measurements.
% Each "pixel" is represented by its range and intensity of reflection.
% It is a 16 bits number whose bits 0-12 define the distance (i.e. the range)
% in cm (a number 0<=r<2^13), and bits 13-15 indicate the intensity 
%( a number 0<=i<8 ).

% We extract range and intensity, here.
%useful masks, for dealing with bits.
mask1FFF = uint16(2^13-1);
maskE000 = bitshift(uint16(7),13)  ;

intensities = bitand(scan,maskE000);

ranges    = single(bitand(scan,mask1FFF))*0.01; 
% Ranges expressed in meters, and unsing floating point format (e.g. single).

% 2D points, expressed in Cartesian. From the sensor's perpective.
angles = [0:360]'*0.5* pi/180 ;         % associated angle, for each individual range in a scan
X = cos(angles).*ranges;
Y = sin(angles).*ranges;    

index_NotInst = find(intensities==0);
set(handle1, 'xdata',X(index_NotInst),'ydata',Y(index_NotInst), 'Color','blue','Marker','.');
ii = find(intensities~=0);
set(handle2, 'xdata',X(ii),'ydata',Y(ii), 'Color','red','Marker','.');

OOIs = ExtractOOIs(ranges,intensities,handle3) ;

if i == 1
    [initialX, initialY] = TansformCoordinate(OOIs, 1);
    set(handle5, 'xdata',initialX,'ydata',initialY, 'Color','blue','Marker','*'); 
    for q = 1:5
        text(initialX(q)+0.2, initialY(q)+0.2, num2str(q), 'FontSize', 8, 'Color', 'blue');
    end
end

if OOIs.Alarming == 1
    disp("There is no OOI or more than 5 #.");
    OOIs.Alarming = 0;
    NumCurOOI = 0;
    xGlobal = 0;
    yGlobal = 0;
    MinDistncIndice = 0;
else
    velocity_i = velocity(i);
    angularRate_i = angularRate(i);
    EKFmethod(velocity_i, angularRate_i, OOIs, i, t);
    [xGlobal, yGlobal,  xGlobalEKF, yGlobalEKF] = TansformCoordinate(OOIs, i); 
    set(handle4, 'xdata',xGlobal,'ydata',yGlobal, 'Color','green','Marker','+', 'MarkerSize', 12); % display dead reckoning results
    set(handle7, 'xdata',xGlobalEKF,'ydata',yGlobalEKF, 'Color','red','Marker','+', 'MarkerSize', 12);
    
    % find the distance between the initial and currents poles
    [~, NumCurOOI] = size(xGlobal);

    distance = zeros(NumCurOOI, 5);
    for a = 1:NumCurOOI
        for j = 1:5
            distance(a, j) = sqrt((initialX(j)-xGlobal(a))^2+(initialY(j)-yGlobal(a))^2);
        end
    end

    MinDistncIndice = zeros(1, NumCurOOI);
    MinDistance = zeros(1, NumCurOOI);
    for a = 1:NumCurOOI
        MinDistncIndice(a) = find(distance(a,:) == min(distance(a,:)));
        MinDistance(a) = distance(a, MinDistncIndice(a));
    end
    
    % for label OOI
    [~, Greater40exist] = find(MinDistance>0.4);
    if Greater40exist ~= 0
        disp('wrong');
%         text(xGlobal(Greater40exist)+0.2, yGlobal(Greater40exist)-0.2, 'string', 'Unexpected', 'FontSize', 12, 'Color', 'blue', 'FontWeight', 'bold');
        set(handle6, 'xdata', xGlobal(Greater40exist),'ydata', yGlobal(Greater40exist), 'Color','red','Marker','o', 'MarkerSize', 15);
        pause(0.1);
    else
       set(handle6, 'xdata', 0,'ydata',0, 'Color','white', 'MarkerSize', 15);
    end    
end



return;
end
   
% This function is to implement PartE
% Transform the local-fram coordinates to global-fram
function [xGlobal, yGlobal,  xGlobalEKF, yGlobalEKF] = TansformCoordinate(r, i) 
    global theta;
    global vehicleLocation;

    %first, r cotains the information about OOI
    x_b = r.Centers(1,:);
    y_b = r.Centers(2,:)+0.46; % in meters
    
    %transform to global convention
    angle = theta(i)*(pi/180) - pi/2;
    x_r = cos(angle).*(x_b)-sin(angle).*(y_b);
    y_r = sin(angle).*(x_b)+cos(angle).*(y_b);
    
    %plus point modeled by kinematic model
    xGlobal = x_r + vehicleLocation.X(i);
    yGlobal = y_r + vehicleLocation.Y(i);
    
    global ContentEKF;
    xGlobalEKF = x_r + ContentEKF.vehicleStateEKF(1);
    yGlobalEKF = y_r + ContentEKF.vehicleStateEKF(2);
    
   global Xe_History;
   Xe_History(1,i) = ContentEKF.vehicleStateEKF(1);
   Xe_History(2,i) = ContentEKF.vehicleStateEKF(2);
return;
end

function initilizationEKF()
global ContentEKF;
ContentEKF.vehicleStateEKF = [ 0; 0;pi/2 ];
ContentEKF.covariance = zeros(3,3);
return;
end

% To manage figure GUI
function MyCallBackA(~,~,x)   
    global flag;
        
    if (x==1)
       flag.flagPause = ~flag.flagPause; %Switch ON->OFF->ON -> and so on.
       return;
    end
    if (x==2)
        
        disp('you pressed "END NOW"');
        uiwait(msgbox('Ooops, you still need to implement this command!','?','modal'));
        
        % students complete this.
        return;
    end
    return;    
end



