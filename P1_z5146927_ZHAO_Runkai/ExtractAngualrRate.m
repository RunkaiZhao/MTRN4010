% Author: Runkai Zhao, z5146927
% Function to implement Part C
% acquire amgular rate from gyroscope

function [angularRate, indexInterest] = ExtractAngualrRate()
    if ~exist('DataFileName','var'), DataFileNameGyroscope ='IMU_dataC'; end
    if ~exist('DataFileName','var'), DataFileNameLaser ='Laser__2C.mat'; end

    load(DataFileNameGyroscope); % return IMU struct
    load(DataFileNameLaser); % return Vel struct

    for i = 1:IMU.N
        t_IMU(i) =  double(IMU.times(i)-IMU.times(1))/10000;
    end

    for i = 1:dataL.N
        t_Laser(i) =  double(dataL.times(i)-dataL.times(1))/10000;
    end

    j = 1;
    for i = 1:IMU.N
        if t_IMU(i) > t_Laser(j) || t_IMU(i) == t_Laser(j)
            indexInterest(j) = i;
            j=j+1;
        end
    
    if j == dataL.N+1
        break;
    end
    end

    for i = 1:j-1
        angularRate(i) = IMU.DATAf(indexInterest(i));
    end
return;
end