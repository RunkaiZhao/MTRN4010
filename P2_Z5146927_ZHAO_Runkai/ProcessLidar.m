% Author: Runkai Zhao, z5146927
% Program: separated solution for AAS, T1.2020, Project2.Part2
% Function for implementing EKF method to acquire the UGV state and
% covariance

function [X, P] = ProcessLidar(X, P, OOIs, v, dt)
   global initial;

    stdDev_rangeMeasurement = 0.2;
    stdDev_bearingMeasurement = 1.5*pi/180;
    stdDev_vel = 0.4;
    stdDev_AngRate = 1.5*pi/180; %input standard deviation
    Q = diag( [ (0.02)^2 ,(0.02)^2 , (0.01*pi/180)^2]) ;
    
    [OOI_Xglobal, OOI_Yglobal] = transformCoorFram(X, OOIs);
    
    %find index
    distance = zeros(OOIs.N, 5);
    for a = 1:OOIs.N
        for j = 1:5
            distance(a, j) = sqrt((initial.X(j)-OOI_Xglobal(a))^2+(initial.Y(j)-OOI_Yglobal(a))^2);
        end
    end
    index = zeros(1, OOIs.N);
    %MinDistance = zeros(1, OOIs.N);
    disp(distance);
    disp(OOIs.N);
    for a = 1:OOIs.N
        index(a) = find(distance(a,:) == min(distance(a,:)));
        %MinDistance(a) = distance(a, MinDistncIndice(a));
    end    
    disp("index: ");
    disp(index)
    
    J = [ [1,0,-dt*v*sin(X(3)) ]  ; [0,1,dt*v*cos(X(3))] ;    [ 0,0,1 ] ] ;
    Pu = [stdDev_vel^2 0;0 stdDev_AngRate^2];
    Fu = [dt*cos(X(3)) 0; dt*sin(X(3)) 0;0 dt];
    P = J*P*J'+ Fu*Pu*Fu' + Q;
    
    %if OOIs.N>0 || OOIs.N<6
        for u=1:OOIs.N
            eDX = (initial.X(index(u)) - X(1)); 
            eDY = (initial.Y(index(u)) - X(2));
            eRange = sqrt( eDX*eDX + eDY*eDY ) ; %   estimated range 
            eBearing = atan2(eDY, eDX)-X(3)+pi/2; % estimated bearing(angle) in radian. within the local frame
            H = [-eDX/eRange , -eDY/eRange , 0; eDY/eRange^2, -eDX/eRange^2, -1]; %2 by 3 matrix
            % Evaluate residual (innovation)
            MeasuredRanges = sqrt((OOI_Xglobal(u)-X(1))^2 + (OOI_Yglobal(u)-X(2))^2);
            MeasuredBearings = atan2(OOI_Yglobal(u)-X(2), OOI_Xglobal(u)-X(1)) - X(3)+pi/2;  
            Z = [MeasuredRanges; MeasuredBearings] - [eRange; eBearing];
            R = [[4*stdDev_rangeMeasurement^2 0];[0 stdDev_bearingMeasurement^2]];
            S = R + H*P*H';
            iS = 1\S;
            %iS = inv(S);
            K = P*H'*iS;
            X = X+K*Z;
            P = P-K*H*P;
       end

return;
end

function [OOI_Xglobal, OOI_Yglobal] = transformCoorFram(X, OOIs)

    x_b = OOIs.Centers(1,:);
    y_b = OOIs.Centers(2,:)+0.46; % in meters
%     disp(x_b);
%     disp(y_b);
    %transform to global convention
    angle = X(3) - pi/2;
    x_r = cos(angle).*(x_b)-sin(angle).*(y_b);
    y_r = sin(angle).*(x_b)+cos(angle).*(y_b);
    
    %plus point modeled by kinematic model
    OOI_Xglobal = x_r + X(1);
    OOI_Yglobal = y_r + X(2);
    
end
        