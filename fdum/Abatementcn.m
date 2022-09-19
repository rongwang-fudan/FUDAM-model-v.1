% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.31
% Establishing the energy system by zone

function [ npv, tip, S ] = Abatementcn( FFlux, L, dpower, dcoef, tempopt, LR, covidyear, iec_cn, calrsav, deffs, techdiff, switcher, econcn, cndata, abtcn, tempdiff, weitzman, dicetfp, tipping, emi_neg, optlrsav, elasmu, rour )

global Egreen realtime S0 cou_iform theta2 alpha fudanccm_exo EFco2 burkecoef tempzone

% cou_iform % 1: id for 222 countries; 2: 2 developing/ 1 developed; 3: region id; 4 OECD; 5 id for 112 countries; 6 pi temperature

T = size(realtime,1);

% inertia of the adjustment of labor / investment allocation
inertia=2;

% initialize the DICE model
cn_num=size(cndata,1);
if fudanccm_exo==3
    % DICE
    TFPd = zeros(cn_num,49); % Total factor productivity in the DICE model
    EMId = zeros(cn_num,49); % Emission Intensity (emission/Y)
    for t=1:49
        for cn=1:cn_num
            TFPd(cn,t) = econcn(t,13,cn)/(econcn(t,10,cn)^alpha)/(cndata(cn,4*49-48+t)/1000)^(1-alpha);
            EMId(cn,t) = econcn(t,20,cn)/econcn(t,13,cn);
        end
    end
end

S=zeros(T,36,cn_num+1);
S(1,1:34,1)=S0(1,1:34);
S(1,35,1)=L(1);
S(1,1:23,2:(cn_num+1))=econcn(1,1:23,1:cn_num);
S(1,35,2:(cn_num+1))=cndata(1:cn_num,4*49-48+1); % Labor in 1971

% Tipping Interaction according to Cai 2016 NCC
tip_para=[0.063 50 15 0; 0.188 1500 10 0.067; 0.104 500 5 0.2; 0.163 50 5 1; 0.053 50 10 0.2];
tip_fij=[0 -0.235 0.125 0.55 0.121; 1.62 0 0.378 0.108 0; 0.107 0.246 0 0 0; 0 0 0 0 0; -0.083 0 0.5 2.059 0];
tip_prob=ones(5,4);
tip=ones(T,5);
t=1;
covid=1;
while realtime(t,1)<2301
    % using the simulation results before abatements
    if realtime(t,1)<2018 % using the historical data
        econ1 = S0(t,1:23);
        clim1 = S0(t,24:34);
        t=t+1;
        S(t,1:34,1) = S0(t,1:34);
        S(t,35,1) = L(t);
        S(t,1:23,2:(cn_num+1)) = econcn(t,1:23,1:cn_num);
        S(t,4,2:(cn_num+1))  = S0(t,4); % MAC of zero-carbon energy 
        S(t,22,2:(cn_num+1))  = S0(t,22);
        S(t,23,2:(cn_num+1))  = S0(t,23);
        S(t,35,2:(cn_num+1)) = econcn(t,24,1:cn_num);        
        %Initializing TFP adn Emission_factor in DICE
        if fudanccm_exo==3
            S(t,1,1)=EMId(1,floor(realtime(t,1))-1970);
            S(t,3,1)=TFPd(1,floor(realtime(t,1))-1970);
            S(t,1,2:(cn_num+1))=EMId(1:cn_num,floor(realtime(t,1))-1970);
            S(t,3,2:(cn_num+1))=TFPd(1:cn_num,floor(realtime(t,1))-1970);
        end
        continue;
    end
    
    %Endogeneous saving rate
    if t<33
        rsav=calrsav(3);
    elseif t<38
        rsav=calrsav(4);
    elseif t<121
        rsav=calrsav(5);
    else
        rsav=optlrsav(t-120);
    end
    
    %Fraction of investment allocated to carbon-emission-free energy: S transition
    cn_yabt=abtcn(1,1);
    cn_fabt=abtcn(1,2);
    if floor(realtime(t,1))>=cn_yabt
        cn_alen=abtcn(1,3);
    else
        cn_alen=100;
    end
    cn_t0=cn_yabt-cn_alen/100*(cn_yabt-2025);
    cn_Nabt=max(0,realtime(t,1)-cn_t0);
    fracinv=(Egreen(45,8)-cn_fabt)*exp(-(cn_Nabt^2)/2/cn_alen/cn_alen)+cn_fabt; 
    
    %Energy cost share (Omega) in the past 20 years
    omega=0; tt=0;
    for t2=1:t
        if (realtime(t,1)-realtime(t2,1))<20
            omega=omega+S(t2,15,1)*realtime(t2,2);
            tt=tt+realtime(t2,2);
        end
    end
    omega = omega/tt;
    
    %Investment for the previous calender year
    investment=0; tt2=0;
    nextyear=floor(realtime(t+1,1)+realtime(t+1,2)/2);
    for i=1:t
        if nextyear==(floor(realtime(i,1)+realtime(i,2)/2)+1)
            investment = investment + rsav * S(i,7,1) * realtime(i,2);
            tt2 = tt2 + realtime(i,2);
        end
    end
    investment=investment/tt2;
    
    %Change of efficiencies in COVID-19
    if realtime(t,1)>covidyear && covid<=11 && realtime(t,1)<(covidyear+1)
        deff=[1+deffs(1,covid),1+deffs(2,covid),1+deffs(3,covid)];
        covid=covid+1;
    else
        deff=[1,1,1];
    end
    
    %Global warming
    tempw = clim1(8); % warming in year t relative to 1850-1990
    tempw0 = S0(49,31); % warming in 2019 relative to 1850-1990
    
    %Tipping function according to Cai 2016 NCC
    for tip_i=1:5
        %Without interaction
        tip_prob(tip_i,1)=tip_prob(tip_i,1)*exp(-tip_para(tip_i,1)/100*max(0,tempw-1)); % Probabiity of no tipping
        if tip_prob(tip_i,1)<0.9
            tip_prob(tip_i,1)=0; % Record the time that a tipping might occur when the probablity is >10%
            if tip_prob(tip_i,3)==1
                tip_prob(tip_i,3)=t; % Time to reach a tipping point
            end
        end
        %Interaction
        interact=0;
        for tip_j=1:5
            if tip_prob(tip_j,2)==0
                interact=interact+tip_fij(tip_j,tip_i);
            end
        end
        tip_prob(tip_i,2)=tip_prob(tip_i,2)*exp(-tip_para(tip_i,1)/100*max(0,tempw-1)*(1+interact)); % Probabiity of no tipping
        if tip_prob(tip_i,2)<0.9
            tip_prob(tip_i,2)=0; % Record the time that a tipping might occur when the probablity is >10%
            if tip_prob(tip_i,4)==1
                tip_prob(tip_i,4)=t; % Time to reach a tipping point
            end
        end
        %Probabiity of no tipping
        tip(t,tip_i)=tip(t-1,tip_i)*exp(-tip_para(tip_i,1)/100*max(0,tempw-1)*(1+interact)); % Probabiity of no tipping
    end
    
    %Climate change catastrophe according to Weitzman 2012
    if weitzman==1
        D = 1 - (1+dcoef * tempw0^dpower+(tempw0/6.081)^6.754)/(1+dcoef * tempw^dpower+(tempw/6.081)^6.754);
    else
        D = 1 - (1+dcoef* tempw0^dpower)/(1+dcoef * tempw^dpower);
    end
    
    %Tipping function according to Cai 2016 NCC
    damage_tip=1;
    if tipping>0
        for tip_i=1:5
            if tip_prob(tip_i,min(2,tipping))==0
                if tipping<=2
                    damage_tip_i = tip_para(tip_i,3)/100 * min( 1,(t-tip_prob(tip_i,min(2,tipping)+2))/tip_para(tip_i,2) ) ;
                else
                    damage_tip_i = tip_para(tip_i,3)/100 * min( 1,(t-tip_prob(tip_i,min(2,tipping)+2))/ 10 ) ;
                end
                damage_tip = damage_tip * (1-damage_tip_i * (1-tip(t,tip_i)));
            end
        end
    end
    D = 1-(1-D)*damage_tip;
    
    %Economy
    econ2 = econdyn(t, L(t), econ1, fracinv, iec_cn(1,:), omega, investment, LR, D, inertia, deff, switcher);
    
    %National loop
    for cn=1:cn_num
        cn_labor = S(t,35,cn+1);
        cn_econ0 = S(t,1:23,cn+1);
        cn_omega = S(t,15,cn+1);
        cn_invest = rsav*S(t,7,cn+1);
        cn_yabt=abtcn(cn,1);
        cn_fabt=abtcn(cn,2);
        cn_growth=abtcn(cn,4); % the prescribed initital rate of per capita GDP growth in 2025
        if floor(realtime(t,1))>=cn_yabt
            cn_alen=abtcn(cn,3);
        else
            cn_alen=100;
        end
        cn_t0=cn_yabt-cn_alen/100*(cn_yabt-2025);
        cn_Nabt=max(0,realtime(t,1)-cn_t0);
        cn_fracinv=(cndata(cn,5*49-48+45)/cndata(cn,2*49-48+45)-cn_fabt)*exp(-(cn_Nabt^2)/2/cn_alen/cn_alen)+cn_fabt;
        
        % Absolute temperature
        cntemp0=tempopt; % global temperature in 2016-2019
        cntempw0=tempw0; % global warming in 2019 relative to 1850-1990
        cntempw=tempw; % global warming in year t relative to 1850-1990
        if tempdiff==3 && cn>1
            cnzone=cou_iform(cndata(cn,1),7);
            cntempw0=tempw0*tempzone(cnzone,1)+tempzone(cnzone,2); % warming by country in 2016-2019
            cntempw=tempw*tempzone(cnzone,1)+tempzone(cnzone,2); % warming by country in year t
        end
        if tempdiff>1 && cn>1
            cntemp0=cou_iform(cndata(cn,1),6); % absolute temperature by country in 2016-2019
        end
        cntemp=cntemp0-cntempw0+cntempw; % absolute temperature by country in year t
        
        % Damage
        if weitzman==1
            D = 1 -(1+dcoef*(abs(cntemp0-tempopt))^dpower+(tempw0/6.081)^6.754) /(1+dcoef*(abs(cntemp-tempopt))^dpower+(tempw/6.081)^6.754);
%         elseif fudanccm_exo==3 && dicetfp==0
%             D = dcoef*tempw^2;
        else
            D = 1 -(1+dcoef*(abs(cntemp0-tempopt))^dpower) /(1+dcoef*(abs(cntemp-tempopt))^dpower);
        end
        D = 1 -(1-D)*damage_tip;
        
        % Calibration of iec to match the growth of economy in the future
        if t<=121 && cn>1
            % EUE
            iec_cn(cn,16)=log10(cn_omega);
            iec_cn(cn,17)=0.01; % like DICE
            iec_cn(cn,28)=cn_omega;
            % ENE
            enerate=((1-0.3)*cn_growth-cn_omega*iec_cn(cn,17))/(1-cn_omega);
            iec_cn(cn,11)=log10(1-cn_omega)-enerate/iec_cn(cn,14);
            iec_cn(cn,24)=iec_cn(cn,11);
            iec_cn(cn,25)=0;
            iec_cn(cn,26)=log10(1-cn_omega);
            iec_cn(cn,27)=enerate;
        end
        
        % setting up TFP in DICE
        if fudanccm_exo==3
            tfprate = 0.016*exp(-0.001*(realtime(t,1)-2019));
            if dicetfp>=1
                if t>121 && cn>1
                    tfprate = (1-0.3)* cn_growth *exp(-0.001*(realtime(t,1)-2019));
                    % Burke 2015 function for TFP
                    if iec_cn(cn,29)==1
                        temp_tfp = burkecoef(1)*(cntemp^2-cntemp0^2)-tempopt*burkecoef(1)*2*(cntemp-cntemp0);
                    else
                        temp_tfp = burkecoef(2)*(cntemp^2-cntemp0^2)-tempopt*burkecoef(1)*2*(cntemp-cntemp0);
                    end
                    tfprate = tfprate - max(-tfprate, temp_tfp); % avoiding bad value in Monte Carlo
                    if dicetfp==1
                        tfprate = max(-0.02, tfprate);
                    else
                        tfprate = max(0, tfprate);
                    end
                end
            end
            cn_econ0(3) = cn_econ0(3) * (1+tfprate)^realtime(t,2); % TFP in DICE
        end

        %Economy
        S(t+1,1:23,cn+1)=econdyn(t, cn_labor, cn_econ0, cn_fracinv, iec_cn(cn,:), cn_omega, cn_invest, LR, D, inertia, deff, switcher);
        
        %Considering learning in DICE for dicetfp=3 only
        if fudanccm_exo==3 && dicetfp<3
            if t<=49
                S(t,4,cn+1) = 550; % The initial price  550 USD (t CO2)-1, which decreases at a constant rate of 0.5% y-1
            end
            if dicetfp<3
                S(t+1,4,cn+1) = S(t,4,cn+1) * (1-0.005)^realtime(t,2);
                S(t+1,19,cn+1) = S(t+1,4,cn+1) * S(t+1,18,cn+1)^(theta2-1);
            end
        end
        
        %Population
        S(t+1,35,cn+1)=S(t,35,cn+1)*L(t+1)/L(t); % labor
        
        %Utility
        if (t+1)>=121 && (t+1)<=400
            S(t+1,36,cn+1) = (((1-rsav)* S(t+1,7,cn+1)/S(t+1,35,cn+1)*1000)^(1-elasmu)/(1-elasmu)) * S(t+1,35,cn+1) /1000  * (1-rour*0.001)^(t-120); % discounted utility        
        end
    end
    
    %Negative emissions
    if tempw>1 && emi_neg>0
        emico2_pre=0; % CO2 emission Gt CO2 without abatement
        for cn=2:cn_num
            if fudanccm_exo==3
                emico2_pre=emico2_pre+S(t+1,13,cn+1) * S(t+1,1,cn+1);
            else
                emico2_pre=emico2_pre+S(t+1,12,cn+1) * EFco2(min(45,t)) * (1-Egreen(min(45,t+1),8));
            end
        end
        fract_neg=emi_neg/emico2_pre; % fraction of negative emissions relative to the fossil fuel emissions without abatement
        for cn=2:cn_num
            if fudanccm_exo==3
                S(t+1,20,cn+1) = S(t+1,20,cn+1) - S(t+1,13,cn+1) * S(t+1,1,cn+1) * S(t+1,18,cn+1) * fract_neg;
            else
                S(t+1,20,cn+1) = S(t+1,20,cn+1) - S(t+1,12,cn+1) * EFco2(min(45,t)) * (1-Egreen(min(45,t+1),8)) * S(t+1,18,cn+1) * fract_neg;
            end
            S(t+1,18,cn+1) = S(t+1,18,cn+1)*(1+fract_neg); % Fraction of emission abatement
            S(t+1,19,cn+1) = S(t+1,4,cn+1)*S(t+1,18,cn+1)^(theta2-1); % carbon price
        end
    end
    
    %CO2 emission Gt CO2
    emico2=sum(S(t+1,20,3:(cn_num+1)),3);
    
    %Additional emissions after a tipping function according to Cai 2016 NCC
    if tipping>0
        for tip_i=2:5
            if tip_prob(tip_i,min(2,tipping))==0
                if (t-tip_prob(tip_i,min(2,tipping)+2))<=tip_para(tip_i,2) || tip_i==5
                    emico2 = emico2 + tip_para(tip_i,4) * 3.666; % Gt CO2
                end
            end
        end
    end
    
    % Global / regional technology diffusion
    if techdiff==1
        ZCE=zeros(6,2);
        for cn=1:cn_num
            if cn>1
                zz=cou_iform(cndata(cn,1),4);
                ZCE(zz,1) = ZCE(zz,1) + S(t,21,cn+1); % green energy PJ
                ZCE(zz,2) = ZCE(zz,2) + S(t+1,21,cn+1); % green energy PJ
                ZCE(5,1) = ZCE(5,1) + S(t,21,cn+1); % green energy PJ
                ZCE(5,2) = ZCE(5,2) + S(t+1,21,cn+1); % green energy PJ
                ZCE(6,1) = ZCE(6,1) + S(t,12,cn+1); % energy PJ
                ZCE(6,2) = ZCE(6,2) + S(t+1,12,cn+1); % energy PJ
            end
        end
        if fudanccm_exo==1
            for cn=1:cn_num
                LRcn = log2(1-min(LR*2,max(LR,LR+log(ZCE(5,1)/3600/50)/log(4)*0.2))); % Effect of decarbonization on the learning rate
                S(t+1,4,cn+1) = S(t,4,cn+1)*(ZCE(6,2)/ZCE(6,1))^(theta2-1); % Effect of energy scale on the marginal abatement cost
                S(t+1,4,cn+1) = S(t+1,4,cn+1)*max(1+0.01*realtime(t,2),ZCE(5,2)/ZCE(5,1))^LRcn; % Effect of learning on the marginal abatement cost
                S(t+1,19,cn+1) = S(t+1,4,cn+1) * S(t+1,18,cn+1)^(theta2-1); % carbon price
            end
        end
    end
    
    % Climate module
    clim2=climdyn(t, clim1, FFlux, emico2);
    t=t+1;
    econ1=econ2;
    clim1=clim2;
    S(t,1:23,1)=econ2(1,1:23);
    S(t,24:34,1)=clim2(1,1:11);
    S(t,35,1)=L(t); % labor
end

for cn=1:(cn_num+1)
    S(:,1,cn) = S(:,1,cn) * 3600; % EUE $/KJ -> $/kWh
    S(:,2,cn) = S(:,2,cn) / 3600; % EPE PJ/(t$)^0.3 / (billion cap)^0.7 -> PWh/(t$)^0.3 / (billion cap)^0.7
    S(:,12,cn) = S(:,12,cn) / 3600; % energy PJ -> PWh
    S(:,21,cn) = S(:,21,cn) / 3600; % cumulative green energy PJ -> PWh
end

npv=sum(sum(S(121:396,36,3:(cn_num+1)),3),1); % global total NPV of utility

end
