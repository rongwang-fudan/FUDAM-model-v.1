% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.10.18

function [ output_att ] = variableS8( S8, calrsav, L )
%   scov: outputs
%   LF: lost function of the simulations

output_att = zeros(396,13*16);
for i=1:13
    output_att(1:396,i)=S8(1:396,7,i); % Q - output
    output_att(1:396,i+13)=S8(1:396,12,i); % E - energy    
    output_att(1:396,i+13*2)=S8(1:396,15,i); % Omega
    output_att(1:396,i+13*4)=S8(1:396,14,i); % Energy price, $/kWh
    output_att(1:396,i+13*5)=S8(1:396,31,i); % Atmospheric temperature, C
    output_att(1:396,i+13*6)=S8(1:396,1,i); % EUE, $/kWh
    output_att(1:396,i+13*7)=S8(1:396,2,i); % EPE, PWh/(t$)^0.3 / (billion cap)^0.7
    output_att(1:396,i+13*8)=S8(1:396,3,i); % ENE, (t$)^0.7/(billion cap)^0.7
    output_att(1:396,i+13*9)=S8(1:396,4,i); % Marginal cost of CO2 emission abatements, $/tCO2
    for j=1:396
        if j<33
            rsav=calrsav(3);
        elseif j<38
            rsav=calrsav(4);
        else
            rsav=calrsav(5);
        end
%         output_att(j,i+13*3) = (1-rsav)*S8(j,7,i)/L(j)*1000; % Per capita consumption
        output_att(j,i+13*3) = (1-rsav)*S8(j,7,i)/L(j)*1000; % Per capita consumption
        
        output_att(j,i+13*10) = rsav*S8(j,7,i)*S8(j,9,i); % Investment allocation to energy, trill $
        output_att(j,i+13*11) = rsav*S8(j,7,i)*S8(j,9,i)*S8(j,16,i); % Investment allocation to carbon-emission-free energy energy, trill $
    end
    output_att(1:396,i+13*12)=S8(1:396,21,i); % Carbon-emission-free energy, PWh
    output_att(1:396,i+13*13)=S8(1:396,20,i); % CO2 emissions, Gt CO2
    output_att(1:396,i+13*14)=S8(1:396,6,i); % climate damage as a percentage of GDP
    output_att(1:396,i+13*15)=S8(1:396,5,i); % abatement cost as a percentage of GDP
end
for i=1:13
    for j=1:16
        output_att(:,i+13*(j-1)) = movavg(output_att(:,i+13*(j-1)),'linear',5);
        output_att(:,i+13*(j-1))=smooth(output_att(:,i+13*(j-1)),'sgolay',2);
    end
end

end
