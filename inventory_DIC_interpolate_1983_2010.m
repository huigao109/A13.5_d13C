close all;clear;clc;
% load data
load('A13.5_1983.mat')

in=find(theta~=-9999 & salinity~=-9999 & depthm>=50 & ...
    sigma0~=-9999 & gamma1~=-9999 & aou~=-9999 ...
    & nitrate~=-9999 & silicate~=-9999 & phosphate~=-9999 ...
    & talk~=-9999 & tco2~=-9999);
Cstar=tco2(in)+117/170*oxygen(in)-0.5*(talk(in)-16/170*oxygen(in));


x(:,1)=theta(in); 
x(:,2)=salinity(in); 
x(:,3)=gamma1(in);
x(:,4)=aou(in);    
x(:,5)=nitrate(in);   
x(:,6)=silicate(in);    
x(:,7)=phosphate(in);
x(:,8)=talk(in);

y=Cstar;

gamma=gamma1;

density=gamma(in);
in1=find(density<26.8);
in2=find(density>=26.8 & density<27.23);
in3=find(density>=27.23 & density<27.5);
in4=find(density>=27.5 & density<28);
in5=find(density>=28 & density<28.27);
in6=find(density>=28.27);

x1989{1}=x(in1,1:8);y1989{1}=y(in1);
x1989{2}=x(in2,1:8);y1989{2}=y(in2);
x1989{3}=x(in3,1:8);y1989{3}=y(in3);
x1989{4}=x(in4,1:8);y1989{4}=y(in4);
x1989{5}=x(in5,1:8);y1989{5}=y(in5);
% x1989{6}=x(in6,1:8);y1989{6}=y(in6);


clear in Cstar theta salinity aou nitrate silicate x y gamma
clear in1 in2 in3 in4 in5 in6 in7 in8 in9 in10 in11 in12 in13 in14
% load data
load('A13.5_2010.mat')

in=find(temperature~=-9999 & theta~=-9999 & salinity~=-9999 & depthm>=50 & ...
    sigma0~=-9999 & gamma1~=-9999 & aou~=-9999 ...
    & nitrate~=-9999 & silicate~=-9999 & phosphate~=-9999 ...
    & talk~=-9999 & tco2~=-9999);

Cstar=tco2(in)+117/170*oxygen(in)-0.5*(talk(in)+16/170*oxygen(in));

% MLR fit equation DIC=a+b*theta+c*salinity+d*AOU+e*silicate+f*phosphate
x(:,1)=theta(in); 
x(:,2)=salinity(in); 
x(:,3)=gamma1(in);
x(:,4)=aou(in);    
x(:,5)=nitrate(in);   
x(:,6)=silicate(in);    
x(:,7)=phosphate(in);
x(:,8)=talk(in);

y=Cstar;
% y=Chat;
% y=tco2(in);
gamma=gamma1;

density=gamma(in);
in1=find(density<26.8);
in2=find(density>=26.8 & density<27.23);
in3=find(density>=27.23 & density<27.5);
in4=find(density>=27.5 & density<28);
in5=find(density>=28 & density<28.27);
in6=find(density>=28.27);

x2013{1}=x(in1,1:8);y2013{1}=y(in1);
x2013{2}=x(in2,1:8);y2013{2}=y(in2);
x2013{3}=x(in3,1:8);y2013{3}=y(in3);
x2013{4}=x(in4,1:8);y2013{4}=y(in4);
x2013{5}=x(in5,1:8);y2013{5}=y(in5);
% x2013{6}=x(in6,1:8);y2013{6}=y(in6);

VarVec=[1 1 0 1 1 1 0 0
        1 1 0 1 0 1 1 0
        1 1 0 1 1 0 0 0
        1 1 0 1 0 0 1 0
        1 1 0 1 0 1 0 0
        0 1 1 1 1 1 0 0
        0 1 1 1 0 1 1 0
        1 0 0 1 1 0 0 1
        1 0 0 1 0 0 1 1
        0 0 1 1 1 0 0 1
        0 0 1 1 0 0 1 1];

% VarVec=[1 1 0 1 1 1 0 0
%         1 1 0 1 0 1 1 0
%         0 1 1 1 1 1 0 0
%         0 1 1 1 0 1 1 0
%         1 0 0 1 1 0 0 1
%         1 0 0 1 0 0 1 1
%         0 0 1 1 1 0 0 1
%         0 0 1 1 0 0 1 1
%         1 1 0 1 1 0 0 0
%         1 1 0 1 0 0 1 0
%         1 1 0 1 0 1 0 0]; 

    
for slab=1:5   
    x1=x1989{1,slab};
    y1=y1989{1,slab};
    x2=x2013{1,slab};
    y2=y2013{1,slab};

for eq=1:11
   UseVars1=x1.*VarVec(eq,:);
   [b1,stats1] = robustfit(UseVars1,y1);  
   coe1(eq,:)=b1; 
   UseVars2=x2.*VarVec(eq,:);
   [b2,stats2] = robustfit(UseVars2,y2);  
   coe2(eq,:)=b2; 
   coediff=coe2-coe1;
   delC(:,eq)=coediff(eq,1)+sum(coediff(eq,2:9).*UseVars2,2);
end    
     md=mean(delC,2);
     clear delC
     dC{slab}=md;
end

d1=dC{1};d2=dC{2};d3=dC{3};
d4=dC{4};d5=dC{5};
% d6=dC{6};

dic(1:length(d1))=d1;
dic(1+length(d1):length(d1)+length(d2))=d2;
dic(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=d3;
dic(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=d4;
dic(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=d5;
% dic(1+length(d1)+length(d2)+length(d3)+...
%     length(d4)+length(d5):length(d1)+length(d2)+length(d3)+length(d4)+length(d5)+length(d6))=d6;

lat=latitudedegrees_north(in);
slat(1:length(d1))=lat(in1);
slat(1+length(d1):length(d1)+length(d2))=lat(in2);
slat(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=lat(in3);
slat(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=lat(in4);
slat(1+length(d1)+length(d2)+length(d3)+length(d4):...
    length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=lat(in5);
% slat(1+length(d1)+length(d2)+length(d3)+...
%     length(d4)+length(d5):length(d1)+length(d2)+length(d3)+length(d4)+length(d5)+length(d6))=lat(in6);


dep=depthm(in);
sdep(1:length(d1))=dep(in1);
sdep(1+length(d1):length(d1)+length(d2))=dep(in2);
sdep(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=dep(in3);
sdep(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=dep(in4);
sdep(1+length(d1)+length(d2)+length(d3)+length(d4):...
    length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=dep(in5);
% sdep(1+length(d1)+length(d2)+length(d3)+...
%     length(d4)+length(d5):length(d1)+length(d2)+length(d3)+length(d4)+length(d5)+length(d6))=dep(in6);

gamma=gamma1(in);
sgamma(1:length(d1))=gamma(in1);
sgamma(1+length(d1):length(d1)+length(d2))=gamma(in2);
sgamma(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=gamma(in3);
sgamma(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=gamma(in4);
sgamma(1+length(d1)+length(d2)+length(d3)+length(d4):...
    length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=gamma(in5);

bd=BotDepthm(in);

clear in x y z
in=find(tco2~=-9999 & talk~=-9999 & phosphate~=-9999 & depthm<=50 ...
        & temperature~=-9999 & salinity~=-9999 &  silicate~=-9999);
par1type =    1; % The first parameter supplied is of type "1", which is "alkalinity"
par1     =   talk(in); % value of the first parameter
par2type =    2; % The first parameter supplied is of type "2", which is "DIC"
par2     =   tco2(in); % value of the second parameter
sal      =   salinity(in); % Salinity of the sample
tempin   =   temperature(in); % Temperature at input conditions
presin   =    pressure(in); % Pressure    at input conditions
tempout  =    0; % Temperature at output conditions 
presout  =    0; % Pressure    at output conditions 
sil      =   silicate(in); % Concentration of silicate  in the sample (in umol/kg)
po4      =   phosphate(in); % Concentration of phosphate in the sample (in umol/kg)
pHscale  =    1; % pH scale at which the input pH is reported ("1" means "Total Scale")  - doesn't matter in this example
k1k2c    =    10; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 )
kso4c    =    3; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

% Do the calculation. See CO2SYS's help for syntax and output format
A=CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);
pco22=A(:,4); % '04' means pco2 input   '19' means pco2 output
load('co2anmmlov2.mat')
pco2diff=co2anmmlov2(find(co2anmmlov2(:,1)==2010),2)-co2anmmlov2(find(co2anmmlov2(:,1)==1983),2);
pco21=pco22-pco2diff;
B=CO2SYS(par1,pco21,1,4,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);
DICe=B(:,2);
Canth=tco2(in)-DICe;

x=latitudedegrees_north(in);
y=depthm(in);
z=Canth;

sg=gamma1(in);

% figure
% scatter(slat,sdep,S,dic,'filled')
% hold on
% scatter(x,y,S,z,'filled')
% colormap(jet)
% caxis([-5 30]);
% colorbar('Fontsize',28,'Fontname','Times New Roman');
% basevalue = 6000;area(lat,bd,basevalue,'FaceColor',[0.5 0.5 0.5])
% set(gca,'ydir','reverse');
% grid on
% set(gca,'Fontsize',28,'XLim',[-55 5],'YLim',[0 6000],'Fontname','Times New Roman')
% xlabel('Latitude','Fontsize',28,'Fontname','Times New Roman');
% ylabel('Depth (m)','Fontsize',28,'Fontname','Times New Roman');
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 6])
%  saveas(gcf,['./A13.56watermassCanth20101983Chat.tif']);


m=20;n=60;
XI=linspace(-42,-32,m);  
YI=linspace(0,120,n); 
sdep(1+length(sdep):length(sdep)+length(y))=y;
slat(1+length(slat):length(slat)+length(x))=x;
dic(1+length(dic):length(dic)+length(z))=z;
sgamma(1+length(sgamma):length(sgamma)+length(sg))=sg;

dept=sdep/50;
% ZI=griddata(slat,sdep,dic,XI,YI.','v4');
ZI=IDW2(slat,dept,dic,XI,YI,-2,'ng',5);

inn=find(gamma1~=-9999);
dlat=latitudedegrees_north(inn);
ddep=depthm(inn)/50;
den=gamma1(inn);
nden=IDW2(dlat,ddep,den,XI,YI,-2,'ng',5);

% bottom
for st=1:max(station)
   sbd=BotDepthm(find(station==st));
   stlat=latitudedegrees_north(find(station==st));
    stbd(st)=sbd(1);
    sttlat(st)=stlat(1);    
end
interpstbd=interp1(sttlat,stbd,XI,'linear','extrap');

sz = size(ZI);
ma=ones(sz);
YI=YI*50;

for aa=1:20
    mas=(ma(:,aa));
    mas(YI>=interpstbd(aa))=0;
    mask(:,aa)=mas;
end

Canth=ZI/(2010-1984)*10;

Cant2=Canth.*mask;

