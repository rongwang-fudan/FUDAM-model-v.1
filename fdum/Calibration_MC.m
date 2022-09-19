% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.5
tic
clear;

nmc=10000;
num1=21; % 21 for 0, 5, 10, ..., 100 percentiles; 11 for 0, 10, 20, ..., 100 percentiles
num2=18; % number of parameters
num3=num1*num2;
obserror=[1	0.1	2	0.4	0.45	0.1	0.1	0.1	0.01	1	1	1	0.04	0.01	0.04];

mcbias=zeros(num3,15);
constrains=1;

dds=zeros(21,1);

for mc=1:num3
    i=mod(mc-1,num1)+1; % percentile 1 - 21
    j=floor((mc-1)/num1)+1; % parameter 1 - 18
    load(strcat('D:\monte carlo',num2str(constrains),'\output_obsm-mc',num2str(mc),'.dat'),'-mat');
    output_obsm(1:45,7)=smooth(output_obsm(1:45,7),5); % LandC obs
    output_obsm(1:45,17)=smooth(output_obsm(1:45,17),5); % Energy price obs
    A=output_obsm(1:28,25:30); output_obsm(:,25:30)=0;
    output_obsm(1:11,25)=A(18:28,1); % Energy model
    output_obsm(1:11,26)=A(18:28,2)./mean(A(5:16,2),1); % Energy model
    output_obsm(1:11,27)=A(18:28,3); % Energy model
    output_obsm(1:11,28)=A(18:28,4)./mean(A(5:16,4),1); % GDP model
    output_obsm(1:11,29)=A(18:28,5); % Energy price model
    output_obsm(1:11,30)=A(18:28,6)./mean(A(5:16,6),1); % Energy price model     
    %
    devmc=zeros(15,1);
    devmc(1)=sqrt(sum((output_obsm(:,1)-output_obsm(:,2)).*(output_obsm(:,1)-output_obsm(:,2)),1)); %  Damage    
    devmc(2)=max(abs(output_obsm(:,3)-output_obsm(:,4)),[],1); %  temperature
    devmc(3)=max(abs(output_obsm(:,5)-output_obsm(:,6)),[],1); %  Atmco2
    devmc(4)=sqrt(mean((output_obsm(10:43,7)-output_obsm(10:43,8)).*(output_obsm(10:43,7)-output_obsm(10:43,8)),1)); %  LandC
    devmc(5)=sqrt(mean((output_obsm(5:43,9)-output_obsm(5:43,10)).*(output_obsm(5:43,9)-output_obsm(5:43,10)),1)); %  OceanC
    devmc(6)=max(abs(output_obsm(:,12)./output_obsm(:,11)-1),[],1); %  Capital 1971-2015
    devmc(7)=max(abs(output_obsm(:,14)./output_obsm(:,13)-1),[],1); %  GDP 1971-2015
    devmc(8)=max(abs(output_obsm(:,16)./output_obsm(:,15)-1),[],1); %  Energy 1971-2015
    devmc(9)=sqrt(mean((output_obsm(1:43,17)-output_obsm(1:43,18)).*(output_obsm(1:43,17)-output_obsm(1:43,18)),1)); %  Energy price 1971-2015    
    devmc(10)=sqrt(mean((output_obsm(1:26,19)-output_obsm(1:26,20)).*(output_obsm(1:26,19)-output_obsm(1:26,20)),1)); %  EUE rate 1971-2015
    devmc(11)=sqrt(mean((output_obsm(1:26,21)-output_obsm(1:26,22)).*(output_obsm(1:26,21)-output_obsm(1:26,22)),1)); %  EPE rate 1971-2015
    devmc(12)=sqrt(mean((output_obsm(1:26,23)-output_obsm(1:26,24)).*(output_obsm(1:26,23)-output_obsm(1:26,24)),1)); %  ENE rate 1971-2015      
    devmc(13)=sqrt(mean((output_obsm(1:11,25)-output_obsm(1:11,26)).*(output_obsm(1:11,25)-output_obsm(1:11,26)),1)); %  E in COVID9
    devmc(14)=sqrt(mean((output_obsm(1:11,27)-output_obsm(1:11,28)).*(output_obsm(1:11,27)-output_obsm(1:11,28)),1)); %  Y in COVID9
    devmc(15)=sqrt(mean((output_obsm(2:5,29)-output_obsm(2:5,30)).*(output_obsm(2:5,29)-output_obsm(2:5,30)),1)); %  P in COVID9
    mcbias(mc,1:15)=devmc(1:15,1);    
    if j==15
        dds(i,1)=mcbias(mc,15);
    end
end
mcbias((14*21+1):(14*21+21),10:12)=mcbias((13*21+1):(13*21+21),10:12); % no reason to constrain elasticity by IEC

LF=zeros(15,18);
para_range = zeros(18,21); % range of parameters
for mc=1:num3
    i=mod(mc-1,num1)+1; % percentile 1 - 21
    j=floor((mc-1)/num1)+1; % parameter 1 - 18
    load(strcat('D:\monte carlo',num2str(constrains),'\parameters-mc',num2str(mc),'.dat'),'-mat');
    for obs=1:15
        LF(obs,j)=max(LF(obs,j),mcbias(mc,obs));
    end
    para_range(j,i) = parameters(j,1);
end

% bias to observations
para_range2=zeros(18,21);
para_rannge_obs=zeros(11,16);
for p=1:18
    para=zeros(21,1);
    paraobs=zeros(15,21);
    p2=zeros(15,1);
    pj=0;
    for mc=1:num3
        i=mod(mc-1,num1)+1; % percentile 1 - 21
        j=floor((mc-1)/num1)+1; % parameter 1 - 18
        if j~=p
            continue;
        end
        load(strcat('D:\monte carlo',num2str(constrains),'\parameters-mc',num2str(mc),'.dat'),'-mat');
        lev=0;
        for obs=1:15
            if mcbias(mc,obs)>obserror(obs)
                lev=lev+1;
            else
                p2(obs,1)=p2(obs,1)+1;
                paraobs(obs,p2(obs,1))=parameters(p,1);
            end
        end
        if lev==0
            pj=pj+1;
            para(pj,1)=parameters(p,1);
        end
    end
    if pj>0
        x=prctile(para(1:pj),[0:5:100],1);
        if x(11)>0
            para_range2(p,1:21)=x(1:21,1);
        else
            para_range2(p,1:21)=x(21:-1:1,1);
        end
    else
        x=errors;
    end
    obs2=0;
    for obs=1:15
        if obs==1 || obs==10 || obs==11 || obs==12
            continue;
        end
        obs2=obs2+1;
        if p2(obs,1)>0
            p3 = prctile(paraobs(obs,1:p2(obs,1)),[5 95],2);
            para_rannge_obs(obs2,p) = p3(2)-p3(1);
        else
            para_rannge_obs(obs2,p) = -999;
        end
    end
end
para_rannge_obs=para_rannge_obs(:,12:17);

uncerts=zeros(5,18);
for para=1:18
    uncerts(1,para)=para_range2(para,11)/para_range(para,11);
    uncerts(2,para)=4;
    uncerts(3,para)=1;
    uncerts(4,para)=abs(para_range(para,11)-para_range(para,2))/abs(para_range(para,11))/1.75;
    uncerts(4,para)=min(uncerts(4,para),abs(para_range2(para,11)-para_range2(para,2))/abs(para_range2(para,11))/1.75);
    uncerts(5,para)=abs(para_range(para,20)-para_range(para,11))/abs(para_range(para,11))/1.75;
    uncerts(5,para)=min(uncerts(5,para),abs(para_range2(para,20)-para_range2(para,11))/abs(para_range2(para,11))/1.75);
end
uncerts2=uncerts';
% error in the regression model
uncerts2(1,1:5)=[1 4 1 0.156596346 0.156596346];
uncerts2(15,1)=1;
uncerts2(16,1:5)=[1 4 1 0.266 0.266];
uncerts2(17,1:5)=[1 4 1 0.269 0.269];





