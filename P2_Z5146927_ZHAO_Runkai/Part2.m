% Author: Runkai Zhao, z5146927
% Program: Solution for AAS, T1.2020, Project2.Part2
% Function for the mian manipulation to data fusion and implement drawings
% Citation: Referred to the sample code provided by the lecturer of MTRN4010.

function Part2()
    % loaf full dataset
    Data =  load('All01.mat'); Data=Data.All;
    L=numel(Data.times);
    
    % create some useful buffer.
    Buf=[];  
    Buf.X = zeros(3,L);
    Buf.X_dr = zeros(3,L);
    
    X = [0;0;pi/2];         % initialize estimated pose.
    X_dr = [0;0;pi/2];
    P = zeros(3,3); % covariance of state
    Q = diag( [ (0.01)^2 ,(0.01)^2 , (1*pi/180)^2]) ; % process model approximation covariance
        
    speed=0; Wz=0;  % initialize current speed and angular rate inputs.


    % get bias of gyro; off-line approach.
    biasWz = mean(Data.Wz(1:2500));  % platform was static, during those samples. (you may see see the data, if you want.)

    % Native timestamps are in uint32, 1 count= 100microseconds.
    tLastPrediction = 0.0001*double(Data.times(1));     %scale to seconds. 
    
    figure(1);clf;
    xlabel('X(m)');ylabel('Y(m)');axis([-5,5,-1,10]);hold on; zoom on;
    handle1 = plot(0,0, '*');
    handle2 = plot(0,0, 'g+');
    handle3 = plot(0,0, 'r+');
    handle4 = plot(0,0, 'go');
    handle5 = plot(0,0, 'ro');
    handle6 = plot(0,0, 'g.');
    handle7 = plot(0,0, 'r.');
    handle8 = plot(0,0, 'b.');
    legend('CoorRef', 'EKFobjects', 'DRobjects', 'EKFlocation', 'DRlocation', 'EKFtrajectory', 'DRtrajectory', 'OtherOOIs', 'location', 'northwest');
    title("EKF method vs DR");
    % create some objects for texting
    for i = 1:5
       textHandle(i) = text(0, 0, '', 'FontSize', 8);
    end 
   
    global initial;
    initial.X = zeros(1,5);
    initial.Y = zeros(1,5);
    scan = Data.Lidar(:,1); 
    OOIs = ProcessScan(scan);
    initial.X = OOIs.Centers(1,:);
    initial.Y = OOIs.Centers(2,:)+0.46;
    DrawInitialMarkers(handle1, initial.X, initial.Y);
    
     for i = 1:5
         s = sprintf('#%d', i); 
        set(textHandle(i), 'String', s, 'Position', [initial.X(i)-0.2; initial.Y(i)-0.3; 0], 'FontSize', 8, 'Color', 'blue');
     end    
    
% Main loop; navigate the data, monotonically. 
%-------------------------------------------
for i=1:5:L
    disp(i);
    what = Data.what(i);        %type of new data sample, 
    t = 0.0001*double(Data.times(i));     %scale to seconds. 
    % Native timestamps are in uint32, 1 count= 100microseconds.
    u = Data.uu(i);                 % index in array of that of data type.

    dt = t-tLastPrediction; tLastPrediction=t;          % horizon since last prediction.
    
    X=DoPrediction(X,dt, speed, Wz);                  % do prediction.
    set(handle4, 'xdata',X(1),'ydata',X(2), 'Color','green','Marker','o', 'MarkerSize', 12);
    Buf.X(:,i)=X;  
    X_dr = DoPrediction(X_dr,dt, speed, Wz);    % just use the process model to esimate the next state without considering noises
    set(handle5, 'xdata',X_dr(1),'ydata',X_dr(2), 'Color','red','Marker','o', 'MarkerSize', 10);
    Buf.X_dr(:,i) = X_dr;
    
    set(handle6, 'xdata',Buf.X(1,:),'ydata',Buf.X(2,:), 'Color','green','Marker','.');
    set(handle7, 'xdata',Buf.X_dr(1,:),'ydata',Buf.X_dr(2,:), 'Color','red','Marker','.');
    
    %------------------------
    % depending on the data type, I will take different actions.
    if (what==3)  %type=3: dead reckoning data (speed and agular rate)
        % new speed and angular rate measurement.
        % remember speed and Wz, to be used in subsequent prediction.
        speed = Data.speeds(u);
        Wz=Data.Wz(u)-biasWz;
        fprintf('New inputs: V=[%.2f]m/s;Wz=[%.2f]deg/sec],t=[%.3f]sec\n',speed,Wz*180/pi,t);
        continue ;
    end
    
    %------------------------
    if what==2  % LiDAR        % type=2: LiDAR scan
        scan = Data.Lidar(:,u); 
        OOIs = ProcessScan(scan);

        if OOIs.N == 0 || OOIs.N > 6 || X(1) > 5 || X(2) > 5
            disp("no valid obervation detected");
        else
            [X, P] = ProcessLidar(X, P, OOIs, speed, dt);
            DrawNotInstOOI(handle8, scan, X);
            DrawEKFMarkers(handle2, OOIs, X);
            DrawDRMarkers(handle3, OOIs, X_dr);
        end
    end
    %------------------------
    pause(0.000000000000000000000001) % pause longer if you want to zoom on;
end              % end main loop.
%-------------------------------------------

end   % end Main() 

% --------------------------------------------------------------
function X=DoPrediction(X0,dt, v, W)  % Implements: one step approximation X(t+dt) = X(t)+dt* dX/dt(t,X,w,v)
    X=X0;
    X = X + dt*[ v*cos(X(3)) ;  v*sin(X(3)) ; W ];
end

function DrawNotInstOOI(handle, scan, X)
    mask1FFF = uint16(2^13-1);
    maskE000 = bitshift(uint16(7),13);
    intensities = bitand(scan,maskE000);
    ranges    = single(bitand(scan,mask1FFF))*0.01; 
    angles = [0:360]'*0.5* pi/180 ;         % associated angle, for each individual range in a scan
    
    xNotInst = cos(angles).*ranges;
    yNotInst = sin(angles).*ranges;    
    indexNotInst = find(intensities==0);
    xNotInst = xNotInst(indexNotInst);
    yNotInst = yNotInst(indexNotInst);
    
    x_b = xNotInst;
    y_b = yNotInst+0.46; % in meters

    %transform to global convention
    angle = X(3) - pi/2;
    x_r = cos(angle).*(x_b)-sin(angle).*(y_b);
    y_r = sin(angle).*(x_b)+cos(angle).*(y_b);
    
    %plus point modeled by kinematic model
    x_NotInst = x_r + X(1);
    y_NotInst = y_r + X(2);
    
     set(handle, 'xdata',x_NotInst,'ydata',y_NotInst, 'Color','blue','Marker','.');
    
end  

function DrawInitialMarkers(handle1, X, Y)
    set(handle1, 'xdata',X,'ydata',Y, 'Color','blue','Marker','*');
end

function DrawEKFMarkers(handle, OOIs, X)
    x_b = OOIs.Centers(1,:);
    y_b = OOIs.Centers(2,:)+0.46; % in meters

    %transform to global convention
    angle = X(3) - pi/2;
    x_r = cos(angle).*(x_b)-sin(angle).*(y_b);
    y_r = sin(angle).*(x_b)+cos(angle).*(y_b);
    
    %plus point modeled by kinematic model
    OOI_Xglobal = x_r + X(1);
    OOI_Yglobal = y_r + X(2);
    
     set(handle, 'xdata',OOI_Xglobal,'ydata',OOI_Yglobal, 'Color','green','Marker','+');
end

function DrawDRMarkers(handle, OOIs, X)
    x_b = OOIs.Centers(1,:);
    y_b = OOIs.Centers(2,:)+0.46; % in meters

    %transform to global convention
    angle = X(3) - pi/2;
    x_r = cos(angle).*(x_b)-sin(angle).*(y_b);
    y_r = sin(angle).*(x_b)+cos(angle).*(y_b);
    
    %plus point modeled by kinematic model
    OOI_Xglobal = x_r + X(1);
    OOI_Yglobal = y_r + X(2);
    
     set(handle, 'xdata',OOI_Xglobal,'ydata',OOI_Yglobal, 'Color','red','Marker','+');
end
