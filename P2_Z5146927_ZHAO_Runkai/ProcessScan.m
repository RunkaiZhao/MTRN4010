% Author: Runkai Zhao, z5146927
% Program: separated solution for AAS, T1.2020, Project2.Part2
% Function for process the scan data captured from laser

function OOIs = ProcessScan(scan)

mask1FFF = uint16(2^13-1);
maskE000 = bitshift(uint16(7),13);
intensities = bitand(scan,maskE000);
ranges    = single(bitand(scan,mask1FFF))*0.01; 

angles = [0:360]'*0.5* pi/180 ;         % associated angle, for each individual range in a scan

OOIs = ExtractFeatures(ranges,intensities);

return;
end