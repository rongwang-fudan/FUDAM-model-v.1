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

% Abating 50% non-CO2 GHGs
load(strcat('..\micmonte\gfrac_sce3.dat'),'-mat');
load(strcat('..\micmonte\temp_sce3.dat'),'-mat');
gfrac_ghg=gfrac;
temp_ghg=temp;
idx_ghg=find(max(gfrac_ghg(1:12,1:mc2,rou),[],1)<1);
mc3_ghg=size(idx_ghg,2);

% Using 50% negative emissions
load(strcat('..\micmonte\gfrac_sce4.dat'),'-mat');
load(strcat('..\micmonte\temp_sce4.dat'),'-mat');
gfrac_neg=gfrac;
temp_neg=temp;
idx_neg=find(max(gfrac_neg(1:12,1:mc2,rou),[],1)<1);
mc3_neg=size(idx_neg,2);

% Central case
load(strcat('..\micmonte\gfrac_sce2.dat'),'-mat');
load(strcat('..\micmonte\temp_sce2.dat'),'-mat');
load(strcat('..\micmonte\para_sce2.dat'),'-mat');

% identify the outliers
idx=find(max(gfrac(1:12,1:mc2,rou),[],1)<1);
mc3=size(idx,2);

subplot(2,3,1);
% MIC against temp
mave1=zeros(4,19); mave1(1,1:19)=[1.2:0.1:3];
for s=1:12
    for mc=1:mc3
        plot(temp(s,idx(mc),s+2),gfrac(s,idx(mc),rou5(s)),'o','MarkerEdgeColor',cmap(s,1:3),'MarkerFaceColor','none','MarkerSize',6); hold on;
        i=floor(temp(s,idx(mc),s+2)*10-11.5)+1;
        if i>=1 && i<=19
            mave1(2,i)=mave1(2,i)+1;
            mave1(3,i)=mave1(3,i)+gfrac(s,idx(mc),rou5(s));
        end
    end
end
mave1(4,1:19)=mave1(3,1:19)./mave1(2,1:19);
plot(mave1(1,1:19),mave1(4,1:19),'LineStyle','-','LineWidth',3,'Color',[0.8 0.8 0.8]); hold on;
%
mave1_ghg=zeros(4,19); mave1_ghg(1,1:19)=[1.2:0.1:3];
for s=1:12
    for mc=1:mc3_ghg
        i=floor(temp_ghg(s,idx_ghg(mc),s+2)*10-11.5)+1;
        if i>=1 && i<=19
            mave1_ghg(2,i)=mave1_ghg(2,i)+1;
            mave1_ghg(3,i)=mave1_ghg(3,i)+gfrac_ghg(s,idx_ghg(mc),rou5(s));
        end
    end
end
mave1_ghg(4,1:19)=mave1_ghg(3,1:19)./mave1_ghg(2,1:19);
plot(mave1_ghg(1,1:19),mave1_ghg(4,1:19),'LineStyle','-','LineWidth',3,'Color',[0.8 0 0.6]); hold on;
%
mave1_neg=zeros(4,19); mave1_neg(1,1:19)=[1.2:0.1:3];
for s=1:12
    for mc=1:mc3_neg
        i=floor(temp_neg(s,idx_neg(mc),s+2)*10-11.5)+1;
        if i>=1 && i<=19
            mave1_neg(2,i)=mave1_neg(2,i)+1;
            mave1_neg(3,i)=mave1_neg(3,i)+gfrac_neg(s,idx_neg(mc),rou5(s));
        end
    end
end
mave1_neg(4,1:19)=mave1_neg(3,1:19)./mave1_neg(2,1:19);
plot(mave1_neg(1,1:19),mave1_neg(4,1:19),'LineStyle','-','LineWidth',3,'Color',[1 0.4 1]);
axis([1 3.3 -0.25 0.5]);

subplot(2,3,4);
% MIC trend agaist temp
mave2=zeros(4,14); mave2(1,1:14)=[1.4:0.1:2.7];
for s=1:8
    for mc=1:mc3
        A(1:5) = gfrac(s:(s+4),idx(mc),rou5(s+2));
        [rs,ms,bs] = regression(tt,A);
        plot(temp(s+2,idx(mc),s+4),ms,'o','MarkerEdgeColor',cmap(s+2,1:3),'MarkerFaceColor','none','MarkerSize',6); hold on;
        i=floor(temp(s+2,idx(mc),s+4)*10-13.5)+1;
        if i>=1 && i<=14
            mave2(2,i)=mave2(2,i)+1;
            mave2(3,i)=mave2(3,i)+ms;
        end
    end
end
mave2(4,1:14)=mave2(3,1:14)./mave2(2,1:14);
plot(mave2(1,1:14),mave2(4,1:14),'LineStyle','-','LineWidth',3,'Color',[0.8 0.8 0.8]); hold on;
%
mave2_ghg=zeros(4,14); mave2_ghg(1,1:14)=[1.4:0.1:2.7];
for s=1:8
    for mc=1:mc3_ghg
        A(1:5) = gfrac_ghg(s:(s+4),idx_ghg(mc),rou5(s+2));
        [rs,ms,bs] = regression(tt,A);
        i=floor(temp_ghg(s+2,idx_ghg(mc),s+4)*10-13.5)+1;
        if i>=1 && i<=14
            mave2_ghg(2,i)=mave2_ghg(2,i)+1;
            mave2_ghg(3,i)=mave2_ghg(3,i)+ms;
        end
    end
end
mave2_ghg(4,1:14)=mave2_ghg(3,1:14)./mave2_ghg(2,1:14);
plot(mave2_ghg(1,1:14),mave2_ghg(4,1:14),'LineStyle','-','LineWidth',3,'Color',[0.8 0 0.6]); hold on;
%
mave2_neg=zeros(4,14); mave2_neg(1,1:14)=[1.4:0.1:2.7];
for s=1:8
    for mc=1:mc3_neg
        A(1:5) = gfrac_neg(s:(s+4),idx_neg(mc),rou5(s+2));
        [rs,ms,bs] = regression(tt,A);
        i=floor(temp_neg(s+2,idx_neg(mc),s+4)*10-13.5)+1;
        if i>=1 && i<=14
            mave2_neg(2,i)=mave2_neg(2,i)+1;
            mave2_neg(3,i)=mave2_neg(3,i)+ms;
        end
    end
end
mave2_neg(4,1:14)=mave2_neg(3,1:14)./mave2_neg(2,1:14);
plot(mave2_neg(1,1:14),mave2_neg(4,1:14),'LineStyle','-','LineWidth',3,'Color',[1 0.4 1]);
axis([1 3.3 -0.02 0.015]);

subplot(2,3,2);
% violin: mic
fpositivemic=zeros(4,2);
for g=1:4
    x=zeros(1,1); i=0;
    for s=1:12
        for mc=1:mc3
            if temp(s,idx(mc),s+2)<1.5
                g2=1;
            elseif temp(s,idx(mc),s+2)<2
                g2=2;
            elseif temp(s,idx(mc),s+2)<2.5
                g2=3;
            else
                g2=4;
            end
            if g==g2
                i=i+1;
                x(i,1)=gfrac(s,idx(mc),rou5(s));
            end
        end
    end
    [y1,y2,u] = ksdensity(x,'kernel','normal');
    y11=max(y1,[],2);
    plot(g*2.5+y1./y11,y2,'LineStyle','-','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
    plot(g*2.5-y1./y11,y2,'LineStyle','-','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
    x2=[(g*2.5-0.5) median(x,1);(g*2.5+0.5) median(x,1)];
    plot(x2(:,1),x2(:,2),'LineStyle','-','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
    
    idxy=find(y2>0);
    fpositivemic(g,1)=sum(y1(idxy),2)/sum(y1,2);
end
%
for g=1:4
    x=zeros(1,1); i=0;
    for s=1:12
        for mc=1:mc3_ghg
            if temp_ghg(s,idx_ghg(mc),s+2)<1.5
                g2=1;
            elseif temp_ghg(s,idx_ghg(mc),s+2)<2
                g2=2;
            elseif temp_ghg(s,idx_ghg(mc),s+2)<2.5
                g2=3;
            else
                g2=4;
            end
            if g==g2
                i=i+1;
                x(i,1)=gfrac_ghg(s,idx_ghg(mc),rou5(s));
            end
        end
    end
%     [y1,y2,u] = ksdensity(x,'kernel','normal');
%     y11=max(y1,[],2);
%     plot(g*2.5+y1./y11,y2,'LineStyle','-','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
%     plot(g*2.5-y1./y11,y2,'LineStyle','-','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
%     x2=[(g*2.5-0.5) median(x,1);(g*2.5+0.5) median(x,1)];
%     plot(x2(:,1),x2(:,2),'LineStyle','-','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
    plot(g*2.5,median(x,1),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',cmap(g*3-2,1:3),'MarkerSize',3); hold on;
end
%
for g=1:4
    x=zeros(1,1); i=0;
    for s=1:12
        for mc=1:mc3_neg
            if temp_neg(s,idx_neg(mc),s+2)<1.5
                g2=1;
            elseif temp_neg(s,idx_neg(mc),s+2)<2
                g2=2;
            elseif temp_neg(s,idx_neg(mc),s+2)<2.5
                g2=3;
            else
                g2=4;
            end
            if g==g2
                i=i+1;
                x(i,1)=gfrac_neg(s,idx_neg(mc),rou5(s));
            end
        end
    end
%     [y1,y2,u] = ksdensity(x,'kernel','normal');
%     y11=max(y1,[],2);
%     plot(g*2.5+y1./y11,y2,'LineStyle','--','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
%     plot(g*2.5-y1./y11,y2,'LineStyle','--','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
%     x2=[(g*2.5-0.5) median(x,1);(g*2.5+0.5) median(x,1)];
%     plot(x2(:,1),x2(:,2),'LineStyle','--','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
    plot(g*2.5,median(x,1),'^','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',cmap(g*3-2,1:3),'MarkerSize',5); hold on;
end
axis([0 12 -0.25 0.5]);

subplot(2,3,5);
% violin: mic trend by temp
for g=1:4
    x=zeros(1,1); i=0;
    for s=1:8
        for mc=1:mc3
            if temp(s+2,idx(mc),s+4)<1.5
                g2=1;
            elseif temp(s+2,idx(mc),s+4)<2
                g2=2;
            elseif temp(s+2,idx(mc),s+4)<2.5
                g2=3;
            else
                g2=4;
            end
            if g==g2
                i=i+1;
                A(1:5) = gfrac(s:(s+4),idx(mc),rou5(s+2));
                [rs,ms,bs] = regression(tt,A);
                x(i,1)=ms;
            end
        end
    end
    [y1,y2,u] = ksdensity(x,'kernel','normal');
    y11=max(y1,[],2);
    plot(g*2.5+y1./y11,y2,'LineStyle','-','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
    plot(g*2.5-y1./y11,y2,'LineStyle','-','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
    x2=[(g*2.5-0.5) median(x,1);(g*2.5+0.5) median(x,1)];
    plot(x2(:,1),x2(:,2),'LineStyle','-','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
    
    idxy=find(y2>0);
    fpositivemic(g,2)=sum(y1(idxy),2)/sum(y1,2);
end
%
for g=1:4
    x=zeros(1,1); i=0;
    for s=1:8
        for mc=1:mc3_ghg
            if temp_ghg(s+2,idx_ghg(mc),s+4)<1.5
                g2=1;
            elseif temp_ghg(s+2,idx_ghg(mc),s+4)<2
                g2=2;
            elseif temp_ghg(s+2,idx_ghg(mc),s+4)<2.5
                g2=3;
            else
                g2=4;
            end
            if g==g2
                i=i+1;
                A(1:5) = gfrac_ghg(s:(s+4),idx_ghg(mc),rou5(s+2));
                [rs,ms,bs] = regression(tt,A);
                x(i,1)=ms;
            end
        end
    end
%     [y1,y2,u] = ksdensity(x,'kernel','normal');
%     y11=max(y1,[],2);
%     plot(g*2.5+y1./y11,y2,'LineStyle','-','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
%     plot(g*2.5-y1./y11,y2,'LineStyle','-','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
%     x2=[(g*2.5-0.5) median(x,1);(g*2.5+0.5) median(x,1)];
%     plot(x2(:,1),x2(:,2),'LineStyle','-','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
    plot(g*2.5,median(x,1),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',cmap(g*3-2,1:3),'MarkerSize',3); hold on;
end
%
for g=1:4
    x=zeros(1,1); i=0;
    for s=1:8
        for mc=1:mc3_neg
            if temp_neg(s+2,idx_neg(mc),s+4)<1.5
                g2=1;
            elseif temp_neg(s+2,idx_neg(mc),s+4)<2
                g2=2;
            elseif temp_neg(s+2,idx_neg(mc),s+4)<2.5
                g2=3;
            else
                g2=4;
            end
            if g==g2
                i=i+1;
                A(1:5) = gfrac_neg(s:(s+4),idx_neg(mc),rou5(s+2));
                [rs,ms,bs] = regression(tt,A);
                x(i,1)=ms;
            end
        end
    end
%     [y1,y2,u] = ksdensity(x,'kernel','normal');
%     y11=max(y1,[],2);
%     plot(g*2.5+y1./y11,y2,'LineStyle','--','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
%     plot(g*2.5-y1./y11,y2,'LineStyle','--','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
%     x2=[(g*2.5-0.5) median(x,1);(g*2.5+0.5) median(x,1)];
%     plot(x2(:,1),x2(:,2),'LineStyle','--','LineWidth',1,'Color',cmap(g*3-2,1:3)); hold on;
    plot(g*2.5,median(x,1),'^','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',cmap(g*3-2,1:3),'MarkerSize',5); hold on;
end
axis([0 12 -0.02 0.015]);

subplot(2,3,3);
frac2=zeros(12,50+11+2); % rou from 1 to 50 and discounting from 0.02, 0.025 to 0.07
for s=1:12
    for ri=1:50
        idx2=find(gfrac(s,idx,ri)<0);
        frac2(s,ri)=size(idx2,2)/mc3; % probability of mic<0
    end
    % rate of GDP per cap
    clog=log(sum(output_gdp((116+s*5):(126+s*5),1:(cn_num-1),1),2)./sum(output_pop((116+s*5):(126+s*5),1:(cn_num-1)),2));
    gc=mean(clog(2:11,1)-clog(1:10,1),1); 
    % discounting 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05 (7), 0.055, 0.06, 0.065, 0.07
    for di=1:11
        ri=min(50,max(1,(di+3)*5-floor(gc*elasmu*1000)));
        frac2(s,di+50)=min(1,frac2(s,ri));
        if di==7
            frac2(s,62)=min(1,size(find(gfrac_ghg(s,idx_ghg,ri)<0),2)/mc3_ghg);
            frac2(s,63)=min(1,size(find(gfrac_neg(s,idx_neg,ri)<0),2)/mc3_neg);
        end
    end
end
for di=2:10
    plot(tx(3:14),frac2(1:12,di+50),'LineStyle','-','LineWidth',1,'Color',cmap(di,1:3)); hold on;
end
plot(tx(3:14),frac2(1:12,7+50),'LineStyle','--','LineWidth',2,'Color',[0.8 0.8 0.8]); hold on;
plot(tx(3:14),frac2(1:12,62),'LineStyle','--','LineWidth',2,'Color',[0.8 0 0.6]); hold on;
plot(tx(3:14),frac2(1:12,63),'LineStyle','--','LineWidth',2,'Color',[1 0.4 1]); hold on;
axis([2035 2090 0 1]);

subplot(2,3,6);
gdmic2=zeros(8,11+2); % discounting 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05 (7), 0.055, 0.06, 0.065, 0.07
for s=1:8
    for di=2:10
        ri=max(1,min(50,rou5(s+2)+(di-7)*5));
        for mc=1:mc3
            A(1:5) = gfrac(s:(s+4),idx(mc),ri);
            [rs,ms,bs] = regression(tt,A);
            if ms<0
                gdmic2(s,di)=gdmic2(s,di)+1/mc3; % probability of dmic/dt<0
            end
        end
    end
    
    % scenario of abating non-CO2 ghgs
    for mc=1:mc3_ghg
        A(1:5) = gfrac_ghg(s:(s+4),idx_ghg(mc),rou5(s+2));
        [rs,ms,bs] = regression(tt,A);
        if ms<0
            gdmic2(s,12)=gdmic2(s,12)+1/mc3_ghg; % probability of dmic/dt<0
        end
    end
    
    % scenario of abating non-CO2 ghgs + negative emissions\
    for mc=1:mc3_neg
        A(1:5) = gfrac_neg(s:(s+4),idx_neg(mc),rou5(s+2));
        [rs,ms,bs] = regression(tt,A);
        if ms<0
            gdmic2(s,13)=gdmic2(s,13)+1/mc3_neg; % probability of dmic/dt<0
        end
    end
end

for di=2:10
    plot(tx(5:12),gdmic2(1:8,di),'LineStyle','-','LineWidth',1,'Color',cmap(di,1:3)); hold on;
end
plot(tx(5:12),gdmic2(1:8,7),'LineStyle','--','LineWidth',2,'Color',[0.8 0.8 0.8]); hold on;
plot(tx(5:12),gdmic2(1:8,12),'LineStyle','--','LineWidth',2,'Color',[0.8 0 0.6]); hold on;
plot(tx(5:12),gdmic2(1:8,13),'LineStyle','--','LineWidth',2,'Color',[1 0.4 1]); hold on;

axis([2045 2080 0 1]);
