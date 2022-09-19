% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.5
tic
clear;

load('files\cndata134.dat','-mat'); cndata=cndata2; clear cndata2;
load('files\cou_iform.dat','-mat'); % 1: id for 222 countries; 2: 2 developing/ 1 developed; 3: 12 region id; 4 OECD; 5 id for 112 countries; 6 pi temperature
cou_iform=cou_iform2; clear cou_iform2;
load('..\mic\output_pop.dat','-mat'); % output_pop=zeros(400,cn_num+1);
cn_num=size(cndata,1);

EndSav=2;

% 0-FudanCCM;
% 1-slow economic growth; 2-fast economic growth
% 3-low climate sensitivity; 4-high climate sensitivity
% 5-learning rate of 0.1; 6-learning rate of 0.3
% 7-optimal temperature of 12C; 8-optimal temperature of 14C
cmap1=[0 0 0;
    51 204 204; 51 204 204; 
    204 0 0; 204 0 0; 
    102 0 255; 102 0 255; 
    204 204 0; 204 204 0;
    ]./255;

% 0-FudanCCM;
% 9-discount 1%;
% 10-discounting 2%
% 11-low elasmu 0.9
% 12-50% abatement of non-CO2 GHGs
% 13-weitzman catastrophe;
% 14-cai catastrophe
% 15-negative emissions
cmap2=[0 0 0;
    0 102 255;
    151 193 255; 
    216 216 216; 
    255 153 204;
    255 102 0;
    234 178 0
    153 204 0 
    ]./255;

% 16-DICE;
% 17-discount 1%;
% 18-discounting 2%
% 19-low elasmu 0.9
% 20-weitzman catastrophe
% 21-cai catastrophe
% 29-negative emissions
cmap3=[0 0 0;
    0 102 255;
    151 193 255; 
    216 216 216; 
    255 102 0;
    234 178 0;
    153 204 0;
    255 153 204;
    ]./255;

% 22-Burke
% 23-Burke-R
% 24-Burke-R with discounting 1%
% 25-Burke-R with low elasmu 0.9
% 26-Burke-R with weitzman catastrophe
% 27-Burke-R with cai catastrophe
% 28-Burke-R with negative emissions
cmap4=[102 0 255;
    0 0 0;
    0 102 255;
    216 216 216; 
    255 102 0;
    234 178 0;
    153 204 0;
    151 193 255; 
    255 153 204;
    ]./255;

%FudanCCM
subplot(2,2,1);
for mc=2:8
    load(strcat('..\mic\output_mic_EndSav',num2str(EndSav),'_MC',num2str(mc+7),'.dat'),'-mat');
    plot(output_mic(cn_num,1:12),'LineStyle','-','LineWidth',1,'Color',cmap2(mc,1:3)); hold on;
end
load(strcat('..\mic\output_mic_EndSav',num2str(EndSav),'_MC0.dat'),'-mat');
plot(output_mic(cn_num,1:12),'LineStyle','-','LineWidth',3,'Color',cmap2(1,1:3)); hold on;

%FudanCCM
subplot(2,2,2);
for mc=2:9
    load(strcat('..\mic\output_mic_EndSav',num2str(EndSav),'_MC',num2str(mc-1),'.dat'),'-mat');
    if mc==3 || mc==4 || mc==7 || mc==9
        plot(output_mic(cn_num,1:12),'LineStyle',':','LineWidth',2,'Color',cmap1(mc,1:3)); hold on;
    else
        plot(output_mic(cn_num,1:12),'LineStyle','-','LineWidth',1,'Color',cmap1(mc,1:3)); hold on;
    end
end
load(strcat('..\mic\output_mic_EndSav',num2str(EndSav),'_MC0.dat'),'-mat');
plot(output_mic(cn_num,1:12),'LineStyle','-','LineWidth',3,'Color',cmap1(1,1:3)); hold on;

% DICE
subplot(2,2,3);
for mc=1:6
    load(strcat('..\mic\output_mic_EndSav',num2str(EndSav),'_MC',num2str(mc+15),'.dat'),'-mat');
    if  mc==1
        plot(output_mic(cn_num,1:12),'LineStyle','-','LineWidth',3,'Color',cmap3(1,1:3)); hold on;
    else
        plot(output_mic(cn_num,1:12),'LineStyle','-','LineWidth',1,'Color',cmap3(mc,1:3)); hold on;
    end
end
% load(strcat('..\mic\output_mic_EndSav',num2str(EndSav),'_MC16.dat'),'-mat');
% plot(output_mic(cn_num,1:12),'LineStyle','-','LineWidth',3,'Color',cmap3(1,1:3)); hold on;
load(strcat('..\mic\output_mic_EndSav',num2str(EndSav),'_MC',num2str(29),'.dat'),'-mat');
plot(output_mic(cn_num,1:12),'LineStyle','-','LineWidth',1,'Color',cmap3(7,1:3)); hold on;
load(strcat('..\mic\output_mic_EndSav',num2str(EndSav),'_MC',num2str(30),'.dat'),'-mat');
plot(output_mic(cn_num,1:12),'LineStyle','-','LineWidth',1,'Color',cmap3(8,1:3)); hold on;

% DICE
subplot(2,2,4);
for mc=1:7
    if mc==4
        continue; % not showing elasmu = 1.45
    end
    load(strcat('..\mic\output_mic_EndSav',num2str(EndSav),'_MC',num2str(mc+21),'.dat'),'-mat');
    if mc==2
        plot(output_mic(cn_num,1:12),'LineStyle','-','LineWidth',3,'Color',cmap4(mc,1:3)); hold on;
    else
        plot(output_mic(cn_num,1:12),'LineStyle','-','LineWidth',1,'Color',cmap4(mc,1:3)); hold on;
    end
end
load(strcat('..\mic\output_mic_EndSav',num2str(EndSav),'_MC',num2str(31),'.dat'),'-mat');
plot(output_mic(cn_num,1:12),'LineStyle','-','LineWidth',1,'Color',cmap4(8,1:3)); hold on;
load(strcat('..\mic\output_mic_EndSav',num2str(EndSav),'_MC',num2str(32),'.dat'),'-mat');
plot(output_mic(cn_num,1:12),'LineStyle','-','LineWidth',1,'Color',cmap4(9,1:3)); hold on;

% subplot(2,3,4);
% load(strcat('..\mic\output_var_EndSav',num2str(EndSav),'_MC',num2str(13),'.dat'),'-mat');
% bar(5-sum(output_var(121+2150-2025,3:12,19:23),3)); % number of tipping events

subplot(2,2,1);
load(strcat('..\mic\output_mic_EndSav1_MC0.dat'),'-mat'); % exo saving
plot(output_mic(cn_num,1:12),'LineStyle','--','LineWidth',2,'Color',cmap2(1,1:3)); hold on;
load(strcat('..\mic\output_mic_EndSav2_MC34.dat'),'-mat'); % neg + tipping
plot(output_mic(cn_num,1:12),'LineStyle',':','LineWidth',2,'Color',cmap2(8,1:3)); hold on;

subplot(2,2,2);
load(strcat('..\mic\output_mic_EndSav1_MC0.dat'),'-mat'); % exo saving
plot(output_mic(cn_num,1:12),'LineStyle','--','LineWidth',2,'Color',cmap2(1,1:3)); hold on;

subplot(2,2,3);
load(strcat('..\mic\output_mic_EndSav1_MC16.dat'),'-mat'); % exo saving
plot(output_mic(cn_num,1:12),'LineStyle','--','LineWidth',2,'Color',cmap3(1,1:3)); hold on;
load(strcat('..\mic\output_mic_EndSav2_MC33.dat'),'-mat'); % neg + tipping
plot(output_mic(cn_num,1:12),'LineStyle',':','LineWidth',2,'Color',cmap3(7,1:3)); hold on;

subplot(2,2,4);
load(strcat('..\mic\output_mic_EndSav1_MC23.dat'),'-mat'); % exo saving
plot(output_mic(cn_num,1:12),'LineStyle','--','LineWidth',2,'Color',cmap4(2,1:3)); hold on;
load(strcat('..\mic\output_mic_EndSav1_MC22.dat'),'-mat'); % exo saving
plot(output_mic(cn_num,1:12),'LineStyle','--','LineWidth',2,'Color',cmap4(1,1:3)); hold on;
load(strcat('..\mic\output_mic_EndSav2_MC35.dat'),'-mat'); % neg + tipping
plot(output_mic(cn_num,1:12),'LineStyle',':','LineWidth',2,'Color',cmap4(7,1:3)); hold on;


micd=zeros(12,3); % rou from 1 to 50 and discounting from 0.02, 0.025 to 0.07
mcss=[0 16 23];
elasmu=1.45;
for simu=1:3
    load(strcat('..\mic\output_mic50_EndSav',num2str(1),'_MC',num2str(mcss(simu)),'.dat'),'-mat');
    load(strcat('..\mic\output_gdp_EndSav',num2str(1),'_MC',num2str(mcss(simu)),'.dat'),'-mat');
    for s=1:12
        % rate of GDP per cap
        clog=log(sum(output_gdp((116+s*5):(126+s*5),1:(cn_num-1),1),2)./sum(output_pop((116+s*5):(126+s*5),1:(cn_num-1)),2));
        gc=mean(clog(2:11,1)-clog(1:10,1),1);
        ri=50-floor(gc*elasmu*1000);
        micd(s,simu)=output_mic50(cn_num,s,ri);
    end
end

subplot(2,2,1);
plot(micd(1:12,1),'LineStyle','-.','LineWidth',2,'Color',[0 0 0]); hold on;

subplot(2,2,2);
plot(micd(1:12,1),'LineStyle','-.','LineWidth',2,'Color',[0 0 0]); hold on;

subplot(2,2,3);
plot(micd(1:12,2),'LineStyle','-.','LineWidth',2,'Color',[0 0 0]); hold on;

subplot(2,2,4);
plot(micd(1:12,3),'LineStyle','-.','LineWidth',2,'Color',[0 0 0]); hold on;   


