tic
clear;

global  elas ESC INT FFlux deltarf fudanccm_exo

%Fudan University model for climate change mitigation % 1 fudanccm endo; 2 fudanccm exo; 3 dice
fudanccm_exo=1; 
%Elasticity of substitution
elas = 0.4;
%CO2 Forcings of equilibrium CO2 doubling (Wm-2)
deltarf = 3.8;
%Equilibrium sensitivity of climate Sherwood SC, et al.. An Assessment of Earth's Climate Sensitivity Using Multiple Lines of Evidence. Rev Geophys. 2020 Dec;58(4):e2019RG000678.
ESC = 3.1/deltarf;
%Time inertia of climate system to reach equilibirum (year)
INT = 53; % 53 in OSCAR; 38 in DICE
%Climate damage function: dpo for power coefficient on temperature, dcoef for damage as a percentage of GDP for 1 degree warming
dpo = 2;
%Data of climate damage: 1 (high damage, unconfirmed) 1999-Climate change policy_ quantifying uncertainties for damages and optimal carbon taxes; 2 (moderate, used by DICE) 2017-A Survey of Global Impacts of Climate Change: Replication, Survey Methods, and a Statistical Analysis
damagedata = 2;
%Learning rate on the cost curve
LR = 0.2;
%Peak population
LA = 11500;
%Year of COVID-19 outbreak
covidyear = 2020;
%year to initiate mitigation
abtyear = 2025;
%Ultimate fraction of CO2 emission abatement (1 for zero emissions)
abtfrac = 1;
%Time (years) to abate CO2 emissions
abtlen = 10;
% switcher for  C1	C2	S1	S2	S3	S4	S5	T1	T2	T3	T4
switcher = ones(1,10);
%Simulation output in scenarios
output_abt = zeros(396,22*3);

%Historical data of economy
Initialset;

%Historical data of climate
InitialsetC;

%Land-use change emissions GtC and radiative forcing for non-CO2
AerosolsLUC;

%Time series of population: L
L = population( LA );

%regression of regional temperature with global temperature in 11 regions
tempz = tempzoneregress( 11 );

%Calibration of climate damage function
[dcoef, xy_damage] = damage( dpo, damagedata );

%Calibration of induced efficiency change
[iec_cn] = Calibration_IEC( cndata, 0 );

%Calibration of equilibrium sensitivity of climate
[output_esc] = Calibration_ESC( FFlux, 1 );

%Calibration of savings rate by capital, energy and ouput
[calrsav, output_cap] = Calibration_CAP( L, iec_cn(1,:), LR, switcher, 1, 0 );


erate=zeros(45,5);
for i=1:45
    i1=max(1,i-2); i2=min(45,i+2);
    x=realtime(i1:i2,1);
    for j=1:5
        if j<5
            y=log(output_cap(i1:i2,26+j)); % K, E, Y, L
        elseif j==5
            y=log(output_cap(i1:i2,28).*Egreen(i1:i2,8)); % Zero-carbon energy
        end
        [r,m,b]=regression(x',y');
        erate(i,j)=m;
    end
end

%National data
cnK = load('files\cnK.txt'); % 175x71 (1971-2019) repeat China
cnE = load('files\cnE.txt'); % 142x50 (1971-2019) repeat China
cnY = load('files\cnY.txt');
cnL = load('files\cnL.txt');
cnP1 = load('files\cnE_total.txt'); % 215x40 (1980-2018) repeat China/US
cnP2 = load('files\cnE_nuclear.txt');
cnP3 = load('files\cnE_renew.txt');
cndata=zeros(142,1+49*5); % KEYL+Ez
i2=0;
for i=1:142
    if i==24 || i==26
        continue;
    end
    for s=1:5
        if s==1
            idk=0;
            for ik=1:175
                if cnK(ik,1)==cnE(i,1) && ik~=38 && ik~=39 && ik~=40 % repeated country
                    x=cnK(ik,23:71);
                    idk=ik; % we get country ID in cnK the same as that in cnE
                end
            end
            if idk==0
                break; % no country ID in cnK is the same as that in cnE
            end
        elseif s==2
            x=cnE(i,2:50);
        elseif s==3
            x=cnY(i,2:50);
        elseif s==4
            x=cnL(i,2:50);
        elseif s==5
            % energy by carbon-zero
            idk=0;
            for ik=1:215
                if cnP1(ik,1)==cnE(i,1) && ik~=41 && ik~=42 && ik~=43 && ik~=203 && ik~=204 && ik~=205 % repeated country
                    x=zeros(1,49);
                    for j=1:49
                        j2=min(40,max(j-8,2));
                        if cnP1(ik,j2)>0 && (cnP2(ik,j2)+cnP3(ik,j2))>0
                            idk=ik;
                            x(j)=(cnP2(ik,j2)+cnP3(ik,j2))/cnP1(ik,j2)*cndata(i2,2*49-48+j);
                        else
                            x(j)=-999;
                        end
                    end
                end
            end
            if idk==0
                for j=1:49
                    cndata(i2,s*49-48+j)=Egreen(min(j,45),8)*cndata(i2,2*49-48+j);
                end
                continue;
            end
        end
        % find the longest time series
        len=zeros(49,1);
        for j=1:49
            for k=1:(50-j)
                if x(j+k-1)>0
                    len(j)=len(j)+1;
                else
                    break;
                end
            end
        end
        [IB, IX] = sort(len, 1); % rank of per capita GDP
        if IB(end)<10 && s~=5
            if s>1
                i2=i2-1; % erase this country from the list
            end
            break;
        end
        if s==1
            i2=i2+1; % add this country to the list
            cndata(i2,1)=cnE(i,1);
        end
        % fill in missing data
        if IX(end)>1
            for j=(IX(end)-1):-1:1
                x(j)=x(j+1)/(1+erate(min(45,j),s));
            end
        end
        if (IX(end)+IB(end)-1)<49
            for j=(IX(end)+IB(end)):49
                x(j)=x(j-1)/(1+erate(min(45,j),s));
            end
        end
        cndata(i2,(s*49-47):(s*49+1))=x(1,1:49);
    end
end
cndata=cndata(1:(i2+1),1:end);

% calibration by the global data
cndata(i2+1,1:end)=sum(cndata(1:i2,:),1);
for s=1:5
    for j=1:49
        if j<=45
            if s==5
                ratio=output_cap(j,28)*Egreen(j,8)/cndata(i2+1,s*49-48+j);
            else
                ratio=output_cap(j,26+s)/cndata(i2+1,s*49-48+j);
            end
            if s==2 || s==5
                ratio=ratio * 3600; % energy PWh -> PJ
            end
        end
        cndata(1:i2,s*49-48+j)=cndata(1:i2,s*49-48+j).*ratio;
    end
end
cndata(i2+1,1:end)=sum(cndata(1:i2,:),1);
cndata(i2+1,1)=0; % global
cndata2=cndata;
cndata2(1,:)=cndata(i2+1,:); % KEYLG
cndata2(2:(i2+1),:)=cndata(1:i2,:);

% update of energy data from IEA new dataset https://www.iea.org/data-and-statistics/data-product/greenhouse-gas-emissions-from-energy-highlights
energycongo=[93.5999787	96.38997894	97.0922657	98.56483491	100.7537395	109.9112091	104.0025379	107.9222619	107.6303053	111.5669977	114.7948796	121.0104848	134.6663953	130.1116442	134.5323028	132.1136328	134.7807292	131.2955377	135.6375716	140.2436674	145.8427562	146.3516027	151.2762946	130.8114325	151.9352572	140.7037156	135.2421911	135.8617515	132.5141736	130.2660362	142.5448338	150.7215342	175.6684566	186.7321153	203.1967998	256.5812199	234.5882455	250.42153	280.6648933	307.8777032	414.1170264	476.3349733	504.5905841	507.8739359	536.1703228	567.072711	569.2932604	597.5931686	617.1134665];
cndata2(27,51:99)=energycongo(1,1:49);

save('files\cndata134.dat','cndata2');







