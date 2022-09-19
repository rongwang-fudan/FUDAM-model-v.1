% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.6.18

function [ output_dam ] = Calibration_DAM( dcoef, dpower, xy_damage )

n=size(xy_damage,1)-59;
output_dam=zeros(n,3);
output_dam(1:n,1)=xy_damage(1:n,1); % T
for i=1:n
    output_dam(i,2)=dcoef * (output_dam(i,1))^dpower; % predicted D
end
output_dam(1:n,3)=xy_damage(1:n,3); % observed D

end

