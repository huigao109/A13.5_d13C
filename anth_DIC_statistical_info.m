close all;clear;clc;
% load data
% load 2010 data
load('F:\carbon_study\South Atlantic\GLOPADv2_2019\A13.5_2010.mat')

% MLR fit equation DIC=a+b*temperature+c*salinity+d*AOU+e*silicate+f*phosphate
in=find(temperature~=-9999  & salinity~=-9999 & gamma1~=-9999 & aou~=-9999 ...
        & nitrate~=-9999 & silicate~=-9999 & phosphate~=-9999 & tco2~=-9999 ...
       & latitudedegrees_north>=-42 & latitudedegrees_north<=-32 & depthm>=50);  

x10(:,1)=temperature(in); 
x10(:,2)=salinity(in); 
x10(:,3)=gamma1(in);
x10(:,4)=aou(in);    
x10(:,5)=nitrate(in);   
x10(:,6)=silicate(in);    
x10(:,7)=phosphate(in);
% x(:,8)=talk(in);
y10=tco2(in);

density=gamma1(in);
in1=find(density<26.8);
in2=find(density>=26.8 & density<27.23);
in3=find(density>=27.23 & density<27.5);
in4=find(density>=27.5 & density<28);
% in5=find(density>=28 & density<28.27);
% in6=find(density>=28.27);
in5=find(density>=28);

x2010{1}=x10(in1,1:7);y2010{1}=y10(in1);
x2010{2}=x10(in2,1:7);y2010{2}=y10(in2);
x2010{3}=x10(in3,1:7);y2010{3}=y10(in3);
x2010{4}=x10(in4,1:7);y2010{4}=y10(in4);
x2010{5}=x10(in5,1:7);y2010{5}=y10(in5);
% x2010{6}=x(in6,1:7);y2010{6}=y(in6);


clear in1 in2 in3 in4 in5 in6 x y pres depth
% load 2020 data
% load('F:\carbon_study\South Atlantic\d13Cv2\A13.5_2020v6.mat')
load('F:\carbon_study\South Atlantic\C13\A13.5_2020v5.mat')

in=find(CTDTMP~=-999 & CTDSAL~=-999 & AOU~=-999 & NITRAT~=-999 & SILCAT~=-999 ...
        & PHSPHT~=-999  & TCARBN~=-999 & cal_depth>=50); 
x20(:,1)=CTDTMP(in);
x20(:,2)=CTDSAL(in);
x20(:,3)=cal_gamma(in);
x20(:,4)=AOU(in);
x20(:,5)=NITRAT(in);
x20(:,6)=SILCAT(in);
x20(:,7)=PHSPHT(in);
% x(:,8)=ALKALI(in);

DIC=TCARBN(in);
pres=CTDPRS(in);
depth=cal_depth(in);

y20=TCARBN(in);   % offset


density=cal_gamma(in);
in1=find(density<26.8);
in2=find(density>=26.8 & density<27.23);
in3=find(density>=27.23 & density<27.5);
in4=find(density>=27.5 & density<28);
% in5=find(density>=28 & density<28.27);
% in6=find(density>=28.27);
in5=find(density>=28);

x2020{1}=x20(in1,1:7);y2020{1}=y20(in1);
x2020{2}=x20(in2,1:7);y2020{2}=y20(in2);
x2020{3}=x20(in3,1:7);y2020{3}=y20(in3);
x2020{4}=x20(in4,1:7);y2020{4}=y20(in4);
x2020{5}=x20(in5,1:7);y2020{5}=y20(in5);
% x2020{6}=x(in6,1:7);y2020{6}=y(in6);


VarVec=[1 1 0 1 1 1 0
        1 1 0 1 0 1 1
        1 1 0 1 1 0 0
        1 1 0 1 0 0 1
        1 1 0 1 0 1 0
        0 1 1 1 1 1 0
        0 1 1 1 0 1 1];

%% five layers    
for slab=1:5   
    x1=x2010{1,slab};
    y1=y2010{1,slab};
    x2=x2020{1,slab};
    y2=y2020{1,slab};

for eq=1:1
   UseVars1=x1.*VarVec(eq,:);
   [b1,stats1] = robustfit(UseVars1,y1);  
   coe1(eq,:)=b1; 
   UseVars2=x2.*VarVec(eq,:);
   [b2,stats2] = robustfit(UseVars2,y2);  
   coe2(eq,:)=b2; 
   coediff=coe2-coe1;
   
   coe11{slab}=coe1; 
   coe22{slab}=coe2; 
   
   resi1=stats1.resid;
   resi2=stats2.resid;
   rmse1(slab,eq)=stats1.ols_s;
   rmse2(slab,eq)=stats2.ols_s;
   mresi1(slab,eq)=nanmean(resi1);
   stdresi1(slab,eq)=nanstd(resi1);
   mresi2(slab,eq)=nanmean(resi2);
   stdresi2(slab,eq)=nanstd(resi2);
   
    pd13C1{eq}=coe1(eq,1)+sum(coe1(eq,2:8).*UseVars1,2);
    pd13C2{eq}=coe2(eq,1)+sum(coe2(eq,2:8).*UseVars2,2);

    [rho1(slab,eq),pval1(slab,eq)] = corr(y1,pd13C1{eq});
    [rho2(slab,eq),pval2(slab,eq)] = corr(y2,pd13C2{eq});
    
    NSE1(slab)=1-(sum((y1-pd13C1{eq}).^2))/(sum((y1-mean(y1)).^2));
    NSE2(slab)=1-(sum((y2-pd13C2{eq}).^2))/(sum((y2-mean(y2)).^2));

   preddC(:,eq)=coe1(eq,1)+sum(coe1(eq,2:8).*UseVars2,2);
   mlrd13C(:,eq)=y2-preddC(:,eq);   % MLR method   % measured-predicted
   delC(:,eq)=coediff(eq,1)+sum(coediff(eq,2:8).*UseVars2,2);
 
end    
     md=mean(delC(:,eq),2);
     md_mlr=mean(mlrd13C(:,eq),2);
     mprd=mean(preddC(:,eq),2);
     clear delC  mlrd13C preddC
     dC{slab}=md;
     dC_mlr{slab}=md_mlr;
     preddelC{slab}=mprd;         
     
end

%% all samples
for eq=1:1
   UseVars1=x10.*VarVec(eq,:);
   [b1,stats1] = robustfit(UseVars1,y10);  
   coe1(eq,:)=b1; 
   UseVars2=x20.*VarVec(eq,:);
   [b2,stats2] = robustfit(UseVars2,y20);  
   coe2(eq,:)=b2; 
   coediff=coe2-coe1;
   
   resi1=stats1.resid;
   resi2=stats2.resid;
   rmse1(eq)=stats1.ols_s;
   rmse2(eq)=stats2.ols_s;
   mresi1(eq)=nanmean(resi1);
   stdresi1(eq)=nanstd(resi1);
   mresi2(eq)=nanmean(resi2);
   stdresi2(eq)=nanstd(resi2);
   
    pd13C1{eq}=coe1(eq,1)+sum(coe1(eq,2:8).*UseVars1,2);
    pd13C2{eq}=coe2(eq,1)+sum(coe2(eq,2:8).*UseVars2,2);

    [rho1(eq),pval1(eq)] = corr(y10,pd13C1{eq});
    [rho2(eq),pval2(eq)] = corr(y20,pd13C2{eq});
    
    NSE1=1-(sum((y10-pd13C1{eq}).^2))/(sum((y10-mean(y10)).^2));
    NSE2=1-(sum((y20-pd13C2{eq}).^2))/(sum((y20-mean(y20)).^2));

   preddC(:,eq)=coe1(eq,1)+sum(coe1(eq,2:8).*UseVars2,2);
   mlrd13C(:,eq)=y20-preddC(:,eq);   % MLR method   % measured-predicted
   delC(:,eq)=coediff(eq,1)+sum(coediff(eq,2:8).*UseVars2,2);
end