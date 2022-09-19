function v = randvar( x, s, z )
%   Author:  Rong Wang // contact rongwang@fudan.edu.cn //
%   Date:  2021.8.7
%   Produce a random value based on the distribution
%   x(1): 1 for normal; 2 for uniform; 3 for log-normal; 4 for relative difference
%   x(2): central
%   x(3): lower
%   x(4): upper
%   s: 1, generate a random value; 2, select the percentile j
%   z: percentile

% xx=randn(10000,1);
% xra=prctile(xx,[0:5:100],1); % du/ds by year
% xra(1)=-1.96;
% xra(end)=1.96;
xra=[-1.96000000000000;-1.64813271378742;-1.28613372071428;-1.03696234321609;-0.844492884450751;-0.676925320813962;-0.527894367333415;-0.402920558568819;-0.265930233109539;-0.149229844668551;-0.0205257114860915;0.113957817061441;0.242432349790607;0.376949571488937;0.508044017944601;0.658551851844843;0.833469925620126;1.03016998482561;1.28025683193193;1.63438823247971;1.96000000000000];

if s==1
    if x(1)==1 || x(1)==4
        j = randn; % normal distribution
        while abs(j)>=1.96
            j = randn; % outliers
        end
        if j>0
            v = 1 + x(4) / x(2) * j;
        else
            v = 1 + x(3) / x(2) * j;
        end
    elseif x(1)==2 
        j = (randi(2001)-1)/1000-1; % uniform distribution
        if j>0
            v = 1 + x(4) / x(2) * j;
        else
            v = 1 + x(3) / x(2) * j;
        end
    elseif x(1)==3 
        j = randn; % log-normal distribution
        while abs(j)>=1.96
            j = randn; % outliers
        end
        if j>0
            v = exp(x(4) * j);
        else
            v = exp(x(3) * j);
        end
    end
else
    if x(1)==1 || x(1)==4 
%         j = (z/50-1) * 1.96; % normal distribution
        j = xra(floor(z/5)+1,1); % normal distribution
        if j>0
            v = 1 + x(4) / x(2) * j;
        else
            v = 1 + x(3) / x(2) * j;
        end
    elseif x(1)==2 
        j = z/50-1; % uniform distribution
        if j>0
            v = 1 + x(4) / x(2) * j;
        else
            v = 1 + x(3) / x(2) * j;
        end
    elseif x(1)==3 
%         j = (z/50-1) * 1.96; % normal distribution
        j = xra(floor(z/5)+1,1); % normal distribution
        if j>0
            v = exp(x(4) * j);
        else
            v = exp(x(3) * j);
        end
    end
end

end

