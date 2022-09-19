tic
clear;

load('files\cndata134.dat','-mat');
cndata=cndata2; clear cndata2;
cou_iform=load('files\country_information.txt');  % 222* 5 + 6 for annual mean temp + 7-18 for monthly mean temp

for i=1:222
    idx=find(cndata(:,1)==i);
    if size(idx,1)>0
        cou_iform(i,5)=idx(1)-1;
    end
end

% T_china=cou_iform(43,6);
% T_brazil=cou_iform(28,6);
% for mm=6:6
%     x=cou_iform(:,mm);
%     for cn=1:222
%         if cou_iform(cn,3)<=3
%             x(cn)=max(20,x(cn));
%         elseif cou_iform(cn,3)==5 || cou_iform(cn,3)==6
%             x(cn)=max(10,min(28,x(cn)));
%         else
%             if x(cn)<-990
%                 x(cn)=13;
%             else
%                 % correction to temperature in Figure 2 
%                 % Burke 2015. Global non-linear effect of temperature on economic production[J]. Nature.
%                 % Brazil 22C and China 14C
%                 x(cn)=max(5,min(28,14+(22-14)/(T_brazil-T_china)*(x(cn)-T_china)));
%             end
%         end
%     end
%     % correction to temperature in Figure 1
%     % Burke 2018. Large potential reduction in economic damages under UN mitigation targets. Nature.
%     % 90% Temperature between 10.5 and 18
%     ts=prctile(x,[5 25 50 75 95]); % du/ds by year
%     for cn=1:222
%         x(cn)=max(10.5,min(18,10.5+(18-10.5)/(ts(5)-ts(1))*(x(cn)-ts(1))));
% %         x(cn)=max(10,min(18,10.5+(16-10.5)/(ts(5)-ts(1))*(x(cn)-ts(1))));
%     end
%     x(43)=12.8; % Figure 2 in Burke 2015. Nature.
% end

% 1. Total Eastern and Southern Africa; 2. Total Northern Africa; 3. Total Western and Central Africa; 4. Total East Asia;
% 5. Total South and South-east Asia; 6. Total Western and Central Asia; 7. Total Europe;
% 8. Total Caribbean; 9. Total Central America; 10. Total North America; 11. Total Oceania; 12. Total South America.
for cn=1:222
    % applying temperature for countries missing data
    if cou_iform(cn,3)<=3 || cou_iform(cn,3)==5 || cou_iform(cn,3)==6 || cou_iform(cn,3)==8 || cou_iform(cn,3)==9
        cou_iform(cn,6)=16;
    elseif cou_iform(cn,3)==4
        cou_iform(cn,6)=16;
    elseif cou_iform(cn,3)==7
        cou_iform(cn,6)=13;
    elseif cou_iform(cn,3)==10
        cou_iform(cn,6)=14;
    elseif cou_iform(cn,3)==11 || cou_iform(cn,3)==12
        cou_iform(cn,6)=15;
    end
    % applying temperature for countries with data for 2016-2019 (max 22 C and min 9 C)
    if cou_iform(cn,19)>0
        cou_iform(cn,6)=min(22,max(9,cou_iform(cn,19)));
    end
    % applying minimal temperature for countries in middle east
    if cou_iform(cn,3)==6
        cou_iform(cn,6)=max(18,cou_iform(cn,6));
    end
end
cou_iform(132,6)=12; % mongolia

cou_iform2=cou_iform(:,1:7);

cou_iform2(:,7)=cou_iform(:,20);

% Data by country
load('files\cndata134.dat','-mat'); cndata=cndata2; clear cndata2;
temp134=zeros(134,1);
for cn=1:134
    temp134(cn,1)=cou_iform(cndata(cn+1,1),6);
end

save('files\cou_iform.dat','cou_iform2');
