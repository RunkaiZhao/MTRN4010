% Author: Runkai Zhao, z5146927
%Program: Solution for AAS, T1.2020, Project2.Part1
%Citation: Referred to the sample code provided by the lecturer of MTRN4010.

function Part1()
% Standard deviation of the error in the angular rate sensor. 
stdDevGyro = 2*pi/180 ;    

% Standard deviation of the error in the speed's measurementsiy
stdDevSpeed = 0.15 ;

% ... errors in the range measurements (25cm, standard dev.)
sdev_rangeMeasurement = 0.25 ;   

Dt=0.05 ;                       % "sample time", 50ms
Li = 5000 ;                     % "experiment" duration, in iterations 
DtObservations=0.250 ;          % laser "sample time" (for observations), 4Hz, approximately

n_usedLanmarks = 5 ; 
global NavigationMap;
NavigationMap = CreateSomeMap(n_usedLanmarks) ;  %creates a artificial map!
    
Xe = [ 0; 0;pi/2 ] ; % expected value
Xdr = [ 0; 0;pi/2 ] ;
P = zeros(3,3) ;   % initial covariance

Xreal_History= zeros(3,Li) ;
Xe_History= zeros(3,Li) ;
XeDR_History= zeros(3,Li) ;

Q = diag( [ (0.01)^2 ,(0.01)^2 , (1*pi/180)^2]) ; %uncertainty of process model

time=0 ;
% initialize the simulator of process (the "real system").
InitSimulation(stdDevSpeed,stdDevGyro,sdev_rangeMeasurement,DtObservations); % create a struct containing info

clc();
disp('Running full simulation, OFF line. Please wait...');

for i = 1:Li
    time = time + Dt;
    SimuPlatform(time,Dt);      % because NOW I do not have a real system, I simulate it.
                                                 % update speed, angular rate, state and current time
    [Noisy_speed,Noisy_GyroZ]=GetProcessModelInputs(); % add noise into the input
    
    Xdr = RunProcessModel(Xdr,Noisy_speed,Noisy_GyroZ,Dt);
    
    % start of EKF
    %Prediction step:
    % covariance P(k|k)=>P(k+1|k)
    J = [ [1,0,-Dt*Noisy_speed*sin(Xe(3))  ]  ; [0,1,Dt*Noisy_speed*cos(Xe(3))] ;    [ 0,0,1 ] ] ;  % 3 by 3 matrix
    Fu = [Dt*cos(Xe(3)), 0; Dt*sin(Xe(3)), 0; 0, Dt]; % 3 by 2 matrix
    Pu = diag([(stdDevSpeed)^2, stdDevGyro^2]); % 2 by 2 matrix, covariance of noise in input
    % then I calculate the new coveraince based on currect stuff(k|k), after the prediction P(K+1|K) = J*P(K|K)*J'+Fu*Pu*Fu'+Q ;
    P = J*P*J'+Fu*Pu*Fu'+Q;   
    %P = J*P*J'+Fu*Pu*Fu'
    %expected value: X(k|k)=>X(k+1|k), inset the currect state into process
    %model equation.
    Xe = RunProcessModel(Xe,Noisy_speed,Noisy_GyroZ,Dt) ; % beside, we assume that the noise is zero mean.
    % so, here/now, the variable "Xe" contains "X^(k+1|k)"   ("X^" means "X hat", i.e. predicted expected value of X(k+1) )
    % and the variable "P" is "P(k+1|k)".     
    % The predition step, for this iteration, is done.
    
    % .. Get range measuremens, if those are available.
    [nDetectedLandmarks,MeasuredRanges, MeasuredBearings, IDs]=GetObservationMeasurements();
    
    % Update steps
    for u = 1:nDetectedLandmarks
        ID = IDs(u);            % landmark ID?    (in "real life", this is provided by the "data association")
            
        % some auxiliary variables.
        eDX = (NavigationMap.landmarks(ID,1)-Xe(1));      % (Xa-X)
        eDY = (NavigationMap.landmarks(ID,2)-Xe(2));      % (Ya-Y)
        eRange = sqrt( eDX*eDX + eDY*eDY ) ; %   estimated range
        eBearing = atan2(eDY, eDX)-Xe(3)+pi/2; % estimated bearing(angle) in radian
        % H matrix, evaluated at prior state X^(k+1|k)
        H = [-eDX/eRange , -eDY/eRange , 0; eDY/eRange^2, -eDX/eRange^2, -1]; %2 by 3 matrix
        % Evaluate residual (innovation)
        Z = [MeasuredRanges(u); MeasuredBearings(u)] - [eRange; eBearing];
        % covariance of the noise in measurements
        R = [0.25^2, 0; 0, (1.5*pi/180)^2];
        % intermediate equations for EKF same as KF
        S = R + H*P*H' ;
        K = P*H'*inv(S);
        % finally, we do it...We obtain  X(k+1|k+1) and P(k+1|k+1)
         Xe = Xe+K*Z ;       % update the  expected value
         P = P-K*H*P ;       % update the Covariance % i.e. "P = P-P*H'*S^-1*H*P"  )
        % update done, sigle EKF done
        
    end
    Xe_History(:,i)    = Xe;
    Xreal_History(:,i) = GetCurrentSimulatedState() ;
    XeDR_History(:,i)  = Xdr ;

end %end of loop
        
fprintf('Done. Showing results, now..\n');
% now, we can see the resutls.
% PLOT some results. 
SomePlots(Xreal_History,Xe_History, XeDR_History, NavigationMap) ;
    
end

function map = CreateSomeMap(n_used)
    n_used = max(0,min(4,n_used));      % accepts no less than 1, no more than 4. 
    
    landmarks = [  [ -40,0 ];[ -0,-20 ];[ 10,10 ] ;[ 30,10 ]  ] ;   
    % you may modify this list, in case you want to add more landmarks to the navigation map.
    
    map.landmarks = landmarks(1:n_used,:) ;
    map.nLandmarks = n_used ;
return ;
end

function InitSimulation(stdDevSpeed,stdDevGyro,sdev_rangeMeasurement,DtObservations)
    global ContextSimulation;
    ContextSimulation.Xreal = [ 0; 0;pi/2 ] ;     % [x;y;phi]
    ContextSimulation.stdDevSpeed = stdDevSpeed;
    ContextSimulation.stdDevGyro = stdDevGyro;
    ContextSimulation.Xreal = [0;0;pi/2];
    ContextSimulation.speed=0;
    ContextSimulation.GyroZ=0;
    ContextSimulation.sdev_rangeMeasurement=sdev_rangeMeasurement;
    ContextSimulation.DtObservations=DtObservations;
    ContextSimulation.timeForNextObservation= 0;
    ContextSimulation.CurrSimulatedTime=0;
return;
end


function  SimuPlatform(time,Dt) % update speed, angular rate, state and current time
    global ContextSimulation;    
    
    % simulate some crazy driver for the car..
    [ContextSimulation.speed,ContextSimulation.GyroZ] = SimuControl(ContextSimulation.Xreal,time) ;      % read kinematic model inputs, ideal ones
    % .........................................
    % simulate one step of the "real system":  Xreal(time)
    ContextSimulation.Xreal = RunProcessModel(ContextSimulation.Xreal,ContextSimulation.speed,ContextSimulation.GyroZ,Dt) ;
    ContextSimulation.CurrSimulatedTime = ContextSimulation.CurrSimulatedTime+Dt;
return;
end

function [speed,GyroZ] = SimuControl(X,t)
    speed = 2 ;                                         % cruise speed, 2m/s  ( v ~ 7km/h), the speed is constant, may need to change later
    GyroZ = 3*pi/180 + sin(0.1*2*pi*t/50)*.02 ;         % some crazy driver moving the steering wheel...
return ;
end

function Xnext=RunProcessModel(X,speed,GyroZ,dt) 
    Xnext = X + dt*[ speed*cos(X(3)) ;  speed*sin(X(3)) ; GyroZ ] ;
return ;
end

function [Noisy_speed,Noisy_GyroZ]=GetProcessModelInputs()
    % .....................................................
    % add noise to simulate real conditions
    % WHY? to be realistic. When we measure things the measurements are polluted with noise, So I simulated that situation by adding random
    % noise to the perfect measurements (the ones I get from the simulated "real" platform.
    global ContextSimulation;
    Noisy_speed =ContextSimulation.speed+ContextSimulation.stdDevSpeed*randn(1) ;
    Noisy_GyroZ =ContextSimulation.GyroZ+ContextSimulation.stdDevGyro*randn(1);
return;
end

function [nDetectedLandmarks,MeasuredRanges, MeasuredBearings, IDs]=GetObservationMeasurements(map)
    global ContextSimulation NavigationMap;       
   
    if ContextSimulation.CurrSimulatedTime<ContextSimulation.timeForNextObservation,
        % no measurements of outputs at this time.
        nDetectedLandmarks=0;
        MeasuredRanges=[];
        MeasuredBearings = [];
        IDs=[];
        return ; 
    end;
    
    % next observation time
    ContextSimulation.timeForNextObservation = ContextSimulation.timeForNextObservation+ContextSimulation.DtObservations;
    
    % get simulated range measurements (actual system), in this case the
    % simulated platform.
        [RealRanges, RealBearing, IDs] = GetMeasurementsFomNearbyLandmarks(ContextSimulation.Xreal,NavigationMap) ;
                % ...................................................
        nDetectedLandmarks = length(RealRanges) ;
        
        if (nDetectedLandmarks<1)       % no detected landmarks...
            MasuredRanges=[]; MeasuredBearings = [];  IDs=[];    return ; 
        end
        
        % Here I corrupt the simulated measurements by adding some random noise (to simulate the noise that would be present in the reality)
        % this is the noise:
        noiseInRange= ContextSimulation.sdev_rangeMeasurement*randn(size(RealRanges));
        noiseInBearing= (1.5*pi/180)*randn(size(RealBearing));
        % here I add it to the perfect ranges' measurements
        MeasuredRanges = RealRanges +  noiseInRange ;
        MeasuredBearings = RealBearing +  noiseInBearing ;
        % so MasuredRanges are the measurements polluted with
        % noise. I get the "perfect measurements" from the simulated
        % platform.
        
        % in real life they arrive already polluted!
        
 return;
end

function [ranges, bearing, IDs] = GetMeasurementsFomNearbyLandmarks(X,map) % adding bearing
    
    if map.nLandmarks>0
        dx= map.landmarks(:,1) - X(1);
        dy= map.landmarks(:,2) - X(2) ;
        ranges = sqrt((dx.*dx + dy.*dy)) ;
        bearing = atan2(dy, dx) - X(3) + pi/2;
        IDs = [1:map.nLandmarks];
    else
        IDs=[];ranges=[]; bearing = [];
    end
    % I simulate I measure/detect all the landmarks, however there can be
    % cases where I see just the nearby ones.
    
return ;
end

function X=GetCurrentSimulatedState()
        global ContextSimulation;
        X=ContextSimulation.Xreal;  
 return;   end



function SomePlots(Xreal_History,Xe_History, Xdr_History, map)



figure(2) ; clf ; hold on ;
plot(Xreal_History(1,:),Xreal_History(2,:),'b') ;
plot(Xe_History(1,:),Xe_History(2,:),'r') ;
plot(Xdr_History(1,:),Xdr_History(2,:),'m') ;

if (map.nLandmarks>0)
    plot(map.landmarks(:,1),map.landmarks(:,2),'*r') ;
    legend({'Real path','EKF Estimated path','DR Estimated path','Landmarks'});
else
    legend({'Real path','EKF Estimated path','DR Estimated path'});
end





%show vectors, for visualizing heading, at a small number of points.
ii = [1:225:length(Xe_History(1,:))] ;
m=10;
 quiver(Xreal_History(1,ii),Xreal_History(2,ii),m*cos(Xreal_History(3,ii)),m*sin(Xreal_History(3,ii)),'b','AutoScale','off','Marker','o' ) ;
 quiver(Xe_History(1,ii),Xe_History(2,ii),m*cos(Xe_History(3,ii)),m*sin(Xe_History(3,ii)),'r','AutoScale','off','Marker','+') ;
 quiver(Xdr_History(1,ii),Xdr_History(2,ii),m*cos(Xdr_History(3,ii)),m*sin(Xdr_History(3,ii)),'m','AutoScale','off','Marker','o' ) ;

axis equal ;

title('Path') ;
zoom on ; grid on; box on;

% --------- plot errors between EKF estimates and the real values
figure(3) ; clf ; 
Xe=Xe_History;
subplot(311) ; plot(Xreal_History(1,:)-Xe(1,:)) ;ylabel('x-xe (m)') ;
title('Performance EKF') ;
subplot(312) ; plot(Xreal_History(2,:)-Xe(2,:)) ;ylabel('y-ye (m)') ;
subplot(313) ; plot(180/pi*(Xreal_History(3,:)-Xe(3,:))) ;ylabel('heading error (deg)') ;

figure(4) ; clf ; 
Xe=Xdr_History;
subplot(311) ; plot(Xreal_History(1,:)-Xe(1,:)) ;ylabel('x-xe (m)') ;
title('Performance Dead Reckoning (usually, not good)') ;
subplot(312) ; plot(Xreal_History(2,:)-Xe(2,:)) ;ylabel('y-ye (m)') ;
subplot(313) ; plot(180/pi*(Xreal_History(3,:)-Xe(3,:))) ;ylabel('heading error (deg)') ;
Xe=[];

return ;
end