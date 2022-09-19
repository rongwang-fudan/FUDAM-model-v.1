% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.5
tic
clear;

load('files\cndata134.dat','-mat'); cndata=cndata2; clear cndata2;
load('files\cou_iform.dat','-mat'); % 1: id for 222 countries; 2: 2 developing/ 1 developed; 3: 12 region id; 4 OECD; 5 id for 112 countries; 6 pi temperature
cou_iform=cou_iform2; clear cou_iform2;
load('..\mic\output_pop.dat','-mat'); % output_pop=zeros(400,cn_num-1);
load('..\mic\output_gdp_EndSav1_MC0.dat','-mat');
cn_num=size(cndata,1);

EndSav=1;
rou=15;
elasmu=1.45;
mc2=1000;
pop=sum(output_pop,2);
pop2=zeros(400,cn_num-1);
for i=121:396
    pop2(i,:)=output_pop(i,:)./pop(i,1);
end
tt=[5:5:25];
tx=[2025:5:2100]; tx=tx';
para0=[3.1	0.2	13	0.166780624	0.24099032	1.243061217];
cmap=jet(12);

rou5=zeros(12,1);
for s=1:12
    % rate of GDP per cap
    clog=log(sum(output_gdp((116+s*5):(126+s*5),1:(cn_num-1),1),2)./sum(output_pop((116+s*5):(126+s*5),1:(cn_num-1)),2));
    gc=mean(clog(2:11,1)-clog(1:10,1),1);
    rou5(s)=min(50,max(1,50-floor(gc*elasmu*1000))); % rou for discounting=5%
end

%FudanCCM
mcmode=2;
gfrac=zeros(12,mc2,50); % global mic
temp=zeros(12,mc2,36);
para=zeros(6,mc2);
for mc=1:mc2
    load(strcat('..\micmonte\output_mic_sce',num2str(mcmode),'_MC',num2str(mc),'.dat'),'-mat');
    load(strcat('..\micmonte\output_temp_sce',num2str(mcmode),'_MC',num2str(mc),'.dat'),'-mat');
    for s=1:12
        temp(s,mc,1:36)=output_temp(1:36,s+2);
        for ri=1:50
            gfrac(s,mc,ri)=output_mic(cn_num,s,ri);
        end
    end
    load(strcat('..\micmonte\parameters_sce',num2str(mcmode),'_MC',num2str(mc),'.dat'),'-mat');
    para(1:6,mc) = parameters(1:6);
end
save(strcat('..\micmonte\gfrac_sce',num2str(mcmode),'.dat'),'gfrac');
save(strcat('..\micmonte\temp_sce',num2str(mcmode),'.dat'),'temp');
save(strcat('..\micmonte\para_sce',num2str(mcmode),'.dat'),'para');

