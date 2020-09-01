% Author: Runkai Zhao, z5146927
%Function for implementing Part D
% to acquire heading angle and current location

function [theta, vehicleLocation] = EstAttitude()

if ~exist('DataFileName','var'), DataFileNameGyroscope ='IMU_dataC'; end
if ~exist('DataFileName','var'), DataFileNameSpeeder ='Speed_dataC'; end
load(DataFileNameGyroscope); % return IMU struct
load(DataFileNameSpeeder); % return Vel struct

for i = 1:IMU.N
    t(i) =  double(IMU.times(i)-IMU.times(1))/10000;
end

%figure(1);
%plot(t,IMU.DATAf(6,:));

timeStamp = find(t==20);
bias = mean(IMU.DATAf(6, 1:timeStamp));

theta(1) = pi/2;
thetaWithoutBias = pi/2;
for i = 1:IMU.N-1
    theta(i+1) = theta(i) + (IMU.DATAf(6,i)-bias)*0.005; %0.005 is time step.
    thetaWithoutBias(i+1) = thetaWithoutBias(i) + (IMU.DATAf(6,i))*0.005;
end

theta = theta.*(180/pi);

velocity = Vel.speeds;
x(1) = 0;
for i = 1:Vel.N-1
    x(i+1) = x(i)+velocity(i)*cos(theta(i)*(pi/180))*0.005;
end

y(1) = 0;
for i = 1:Vel.N-1
    y(i+1) = y(i)+velocity(i)*sin(theta(i)*(pi/180))*0.005;
end

vehicleLocation.X = x;
vehicleLocation.Y = y;
return;
end