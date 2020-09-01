% Author: Runkai Zhao, z5146927
% Function to implement Part A
% get all important properties of OOI

function r = ExtractOOIs(ranges,intensities,handle3)
    r.N = 0;
    r.Centers = [];
    r.Sizes   = [];
    
    i = 1;
    counter = 1;
    stop = 0;
    index = 1;
    
    while 1
       
        %% to get the start pixel of OOI. When the intensity is not zero or exceed internsities column, jump out from the loop
        while intensities(i) == 0
            i=i+1;
            counter = counter+1;
            if(counter > 361) % cannot find any OOI
                    stop = 1;
                    break;
            else
                    stop = 0;
            end
        end
        
        if (stop==1)
            stop = 0;
            break
        end
        Start(index)=i;
        
        %% get the end pixel of OOI
        while intensities(i)~=0
            i = i+1;
            counter = counter+1;
            if counter > 361
                break;
            end
        end
        End(index)=i-1;
        
        %% move to next loop
        index=index+1;
        
        %% condition to terminate this loop
        if i - 1 == 361
            break
        end
    end
    
    %% check if there is no OOI
    checker = ones(1,361);
    if (checker*double(intensities) == 0)
        disp("There is no OOI.");
        r.Alarming = 1;
    else
        r.Alarming = 0;
       
    %% the number of OOI---------------------------------------------------
        [~, r.N] = size(Start);
        if (r.N > 5) 
            r.Alarming = 1;
        end
    %% the average of the ranges--------------------------------------------
        for num=1:r.N
            AvrgDistances(num) = mean(ranges(Start(num):End(num)));
        end
        
        r.Ranges =  AvrgDistances;
    %% the average of the angle---------------------------------------------
        angles = [0:360]'*0.5* pi/180 ; % in radian
        for num=1:r.N
        AvrgAngles(num) = mean(angles(Start(num):End(num)));
        end
        
        r.angles = AvrgAngles;
    %% find each center of OOI----------------------------------------------
        for num=1:r.N
        r.Centers(1,num) = cos(AvrgAngles(num)).*AvrgDistances(num);
        r.Centers(2,num) = sin(AvrgAngles(num)).*AvrgDistances(num);
        end
        
     %% find adjusted ranges   
        for num=1:r.N
        adjustedRanges(num) = sqrt(r.Centers(1, num)^2 + (r.Centers(2, num) + 0.46)^2);
        end
        
        r.adjustedRanges = adjustedRanges;
    %% plot out the centers
        set(handle3, 'xdata',r.Centers(1,:),'ydata',r.Centers(2,:), 'Color','green','Marker','*');
        %display(r.N);
    end
return;
end