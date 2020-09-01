% Author: Runkai Zhao, z5146927
%Function for implementing EKF estimation

function EKFmethod(velocity, angularRate, OOIs, i, Dt)
% global vehicleStateEKF
% global covariance;
global ContentEKF;
vehicleStateEKF = ContentEKF.vehicleStateEKF;
covariance = ContentEKF.covariance;

% the input state have been inserted as arguments
stdDev_vel = 0.4;
stdDev_AngRate = 1.5*pi/180;
stdDev_rangeMeasurement = 0.2;
stdDev_bearingMeasurement = 1.5*pi/180;
Q = diag( [ (0.01)^2 ,(0.01)^2 , (1*pi/180)^2]) ; %uncertainty of process model

% start of EKF
%Prediction step:
% covariance P(k|k)=>P(k+1|k)
J = [ [1,0,-Dt*velocity*sin(vehicleStateEKF(3))]  ; [0,1,Dt*velocity*cos(vehicleStateEKF(3))] ;    [ 0,0,1 ] ] ;  % 3 by 3 matrix
Fu = [Dt*cos(vehicleStateEKF(3)), 0; Dt*sin(vehicleStateEKF(3)), 0; 0, Dt]; % 3 by 2 matrix
Pu = diag([(stdDev_vel)^2, stdDev_AngRate^2]); % 2 by 2 matrix, covariance of noise in input
% then I calculate the new coveraince based on currect stuff(k|k), after the prediction P(K+1|K) = J*P(K|K)*J'+Fu*Pu*Fu'+Q ;
covariance = J*covariance*J'+Fu*Pu*Fu'+Q;   
%P = J*P*J'+Fu*Pu*Fu'
%expected value: X(k|k)=>X(k+1|k), inset the currect state into process
%model equation.
vehicleStateEKF = RunProcessModel(vehicleStateEKF, velocity, angularRate, Dt) ; % beside, we assume that the noise is zero mean.
% so, here/now, the variable "Xe" contains "X^(k+1|k)"   ("X^" means "X hat", i.e. predicted expected value of X(k+1) )
% and the variable "P" is "P(k+1|k)".     
% The predition step, for this iteration, is done.

% Update steps
for u = 1:OOIs.N
        %ID = IDs(u);            % landmark ID?    (in "real life", this is provided by the "data association")
            
        % some auxiliary variables.
        eDX = OOIs.Centers(1,u)-vehicleStateEKF(1);      % (Xa-X)
        eDY = OOIs.Centers(2,u)-vehicleStateEKF(2);      % (Ya-Y)
        eRange = sqrt( eDX*eDX + eDY*eDY ) ; %   estimated range
        eBearing = atan2(eDY, eDX)-vehicleStateEKF(3)+pi/2; % estimated bearing(angle) in radian
        % H matrix, evaluated at prior state X^(k+1|k)
        H = [-eDX/eRange , -eDY/eRange , 0; eDY/eRange^2, -eDX/eRange^2, -1]; %2 by 3 matrix
        % Evaluate residual (innovation)
        Z = [OOIs.adjustedRanges(u); OOIs.angles(u)] - [eRange; eBearing];
        % covariance of the noise in measurements
        R = [stdDev_rangeMeasurement^2, 0; 0, stdDev_bearingMeasurement^2];
        % intermediate equations for EKF same as KF
        S = R + H*covariance*H' ;
        iS = inv(S);
        K = covariance*H'*iS;
        % finally, we do it...We obtain  X(k+1|k+1) and P(k+1|k+1)
         vehicleStateEKF = vehicleStateEKF+K*Z ;       % update the  expected value
         covariance = covariance-K*H*covariance ;       % update the Covariance % i.e. "P = P-P*H'*S^-1*H*P"  )
        % update done, sigle EKF done
        
end

ContentEKF.vehicleStateEKF = vehicleStateEKF ;
ContentEKF.covariance = covariance;

% global Xe_History;
% Xe_History(:,i) = vehicleStateEKF;

return;
end


function Xnext=RunProcessModel(X,speed,GyroZ,dt) 
    Xnext = X + dt*[ speed*cos(X(3)) ;  speed*sin(X(3)) ; GyroZ ] ;
return ;
end
