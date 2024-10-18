close all;clear;clc;
% load data
% load 2010 data

load('A13.5_2010.mat')

in=find(temperature~=-9999  & salinity~=-9999 & gamma1~=-9999 & aou~=-9999 ...
        & nitrate~=-9999 & silicate~=-9999 & phosphate~=-9999 & c13~=-9999 ...
       & latitudedegrees_north>=-42 & latitudedegrees_north<=-32);  

x(:,1)=temperature(in); 
x(:,2)=salinity(in); 
x(:,3)=gamma1(in);
x(:,4)=aou(in);    
x(:,5)=nitrate(in);   
x(:,6)=silicate(in);    
x(:,7)=phosphate(in);
% x(:,8)=talk(in);
y=c13(in);

density=gamma1(in);
in1=find(density<26.8);
in2=find(density>=26.8 & density<27.23);
in3=find(density>=27.23 & density<27.5);
in4=find(density>=27.5 & density<28);
% in5=find(density>=28 & density<28.27);
% in6=find(density>=28.27);
in5=find(density>=28);

x2010{1}=x(in1,1:7);y2010{1}=y(in1);
x2010{2}=x(in2,1:7);y2010{2}=y(in2);
x2010{3}=x(in3,1:7);y2010{3}=y(in3);
x2010{4}=x(in4,1:7);y2010{4}=y(in4);
x2010{5}=x(in5,1:7);y2010{5}=y(in5);
% x2010{6}=x(in6,1:7);y2010{6}=y(in6);


clear in1 in2 in3 in4 in5 in6 x y pres depth
% load 2020 data
load('A13.5_2020v6.mat')

in=find(CTDTMP~=-999 & CTDSAL~=-999 & AOU~=-999 & NITRAT~=-999 & SILCAT~=-999 ...
        & PHSPHT~=-999  & DELC13~=-999); 

x(:,1)=CTDTMP(in);
x(:,2)=CTDSAL(in);
x(:,3)=cal_gamma(in);
x(:,4)=AOU(in);
x(:,5)=NITRAT(in);
x(:,6)=SILCAT(in);
x(:,7)=PHSPHT(in);
% x(:,8)=ALKALI(in);

DIC=TCARBN(in);
pres=CTDPRS(in);
depth=cal_depth(in);

y=DELC13(in)+0.07;   % offset
% y=DELC13(in);  


density=cal_gamma(in);
in1=find(density<26.8);
in2=find(density>=26.8 & density<27.23);
in3=find(density>=27.23 & density<27.5);
in4=find(density>=27.5 & density<28);
% in5=find(density>=28 & density<28.27);
% in6=find(density>=28.27);
in5=find(density>=28);

x2020{1}=x(in1,1:7);y2020{1}=y(in1);
x2020{2}=x(in2,1:7);y2020{2}=y(in2);
x2020{3}=x(in3,1:7);y2020{3}=y(in3);
x2020{4}=x(in4,1:7);y2020{4}=y(in4);
x2020{5}=x(in5,1:7);y2020{5}=y(in5);
% x2020{6}=x(in6,1:7);y2020{6}=y(in6);


VarVec=[1 1 0 1 1 1 0
        1 1 0 1 0 1 1
        1 1 0 1 1 0 0
        1 1 0 1 0 0 1
        1 1 0 1 0 1 0
        0 1 1 1 1 1 0
        0 1 1 1 0 1 1];

    
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
   
   preddC(:,eq)=coe1(eq,1)+sum(coe1(eq,2:8).*UseVars2,2);
   mlrd13C(:,eq)=y2-preddC(:,eq);   % MLR method   % measured-predicted
   delC(:,eq)=coediff(eq,1)+sum(coediff(eq,2:8).*UseVars2,2);
 
end    
     md=mean(delC(:,eq),2);
     md_mlr=mean(mlrd13C(:,eq),2);
     mprd=mean(preddC(:,eq),2);
%      stdC=nanstd(delC(:,eq),2);
%      std_mlr=nanstd(mlrd13C(:,eq),2);
%      stdprd=nanstd(preddC(:,eq),2);
     clear delC  mlrd13C preddC
     dC{slab}=md;
     dC_mlr{slab}=md_mlr;
     preddelC{slab}=mprd;     
     
%      stddC{slab}=stdC;
%      stdC_mlr{slab}=std_mlr;
%      stdpreddelC{slab}=stdprd;  
     
end

% eMLR anthropogenic d13C change
d1=dC{1};d2=dC{2};d3=dC{3};d4=dC{4};d5=dC{5};
d13C_eMLR(1:length(d1))=d1;
d13C_eMLR(1+length(d1):length(d1)+length(d2))=d2;
d13C_eMLR(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=d3;
d13C_eMLR(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=d4;
d13C_eMLR(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=d5;

% MLR anthropogenic d13C change
dd1=dC_mlr{1};dd2=dC_mlr{2};dd3=dC_mlr{3};dd4=dC_mlr{4};dd5=dC_mlr{5};
d13C_MLR(1:length(d1))=dd1;
d13C_MLR(1+length(d1):length(d1)+length(d2))=dd2;
d13C_MLR(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=dd3;
d13C_MLR(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=dd4;
d13C_MLR(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=dd5;


%  predicted d13C 2020
pd1=preddelC{1};pd2=preddelC{2};pd3=preddelC{3};pd4=preddelC{4};pd5=preddelC{5};
pd13C(1:length(d1))=pd1;
pd13C(1+length(d1):length(d1)+length(d2))=pd2;
pd13C(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=pd3;
pd13C(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=pd4;
pd13C(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=pd5;


% % eMLR anthropogenic d13C change std
% sd1=stddC{1};sd2=stddC{2};sd3=stddC{3};sd4=stddC{4};sd5=stddC{5};
% stdd13C_eMLR(1:length(d1))=sd1;
% stdd13C_eMLR(1+length(d1):length(d1)+length(d2))=sd2;
% stdd13C_eMLR(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=sd3;
% stdd13C_eMLR(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=sd4;
% stdd13C_eMLR(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=sd5;
% 
% % MLR anthropogenic d13C change std
% sdd1=stdC_mlr{1};sdd2=stdC_mlr{2};sdd3=stdC_mlr{3};sdd4=stdC_mlr{4};sdd5=stdC_mlr{5};
% stdd13C_MLR(1:length(d1))=sdd1;
% stdd13C_MLR(1+length(d1):length(d1)+length(d2))=sdd2;
% stdd13C_MLR(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=sdd3;
% stdd13C_MLR(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=sdd4;
% stdd13C_MLR(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=sdd5;
% 
% 
% %  predicted d13C 2020 std
% spd1=stdpreddelC{1};spd2=stdpreddelC{2};spd3=stdpreddelC{3};spd4=stdpreddelC{4};spd5=stdpreddelC{5};
% stdpd13C(1:length(d1))=spd1;
% stdpd13C(1+length(d1):length(d1)+length(d2))=spd2;
% stdpd13C(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=spd3;
% stdpd13C(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=spd4;
% stdpd13C(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=spd5;

% arrange by station number
sn=STNNBR(in);
ssn1=sn(in1);ssn2=sn(in2);ssn3=sn(in3);
ssn4=sn(in4);ssn5=sn(in5);

ssn(1:length(d1))=ssn1;
ssn(1+length(d1):length(d1)+length(d2))=ssn2;
ssn(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=ssn3;
ssn(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=ssn4;
ssn(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=ssn5;

ins1=find(ssn==1);
ins2=find(ssn==2);
ins3=find(ssn==3);
ins4=find(ssn==4);
ins5=find(ssn==5);
ins6=find(ssn==6);
ins7=find(ssn==7);
ins8=find(ssn==8);

% predicted 2020 d13C for each station 
pd13Cs1=pd13C(ins1);
pd13Cs2=pd13C(ins2);
pd13Cs3=pd13C(ins3);
pd13Cs4=pd13C(ins4);
pd13Cs5=pd13C(ins5);
pd13Cs6=pd13C(ins6);
pd13Cs7=pd13C(ins7);
pd13Cs8=pd13C(ins8);

% each station depth
sdep1=depth(in1);sdep2=depth(in2);sdep3=depth(in3);
sdep4=depth(in4);sdep5=depth(in5);
sdep(1:length(d1))=sdep1;
sdep(1+length(d1):length(d1)+length(d2))=sdep2;
sdep(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=sdep3;
sdep(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=sdep4;
sdep(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=sdep5;

deps1=sdep(ins1);
deps2=sdep(ins2);
deps3=sdep(ins3);
deps4=sdep(ins4);
deps5=sdep(ins5);
deps6=sdep(ins6);
deps7=sdep(ins7);
deps8=sdep(ins8);

% each station pressure
spres1=pres(in1);spres2=pres(in2);spres3=pres(in3);
spres4=pres(in4);spres5=pres(in5);
spres(1:length(d1))=spres1;
spres(1+length(d1):length(d1)+length(d2))=spres2;
spres(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=spres3;
spres(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=spres4;
spres(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=spres5;

pres1=spres(ins1);
pres2=spres(ins2);
pres3=spres(ins3);
pres4=spres(ins4);
pres5=spres(ins5);
pres6=spres(ins6);
pres7=spres(ins7);
pres8=spres(ins8);

% measured 2020 d13C 
md13C(1:length(d1))=y2020{1};
md13C(1+length(d1):length(d1)+length(d2))=y2020{2};
md13C(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=y2020{3};
md13C(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=y2020{4};
md13C(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=y2020{5};

md13Cs1=md13C(ins1);
md13Cs2=md13C(ins2);
md13Cs3=md13C(ins3);
md13Cs4=md13C(ins4);
md13Cs5=md13C(ins5);
md13Cs6=md13C(ins6);
md13Cs7=md13C(ins7);
md13Cs8=md13C(ins8);

%  anthropogenic d13C change for each station
d13Cs1_eMLR=d13C_eMLR(ins1);
d13Cs2_eMLR=d13C_eMLR(ins2);
d13Cs3_eMLR=d13C_eMLR(ins3);
d13Cs4_eMLR=d13C_eMLR(ins4);
d13Cs5_eMLR=d13C_eMLR(ins5);
d13Cs6_eMLR=d13C_eMLR(ins6);
d13Cs7_eMLR=d13C_eMLR(ins7);
d13Cs8_eMLR=d13C_eMLR(ins8);

d13Cs1_MLR=d13C_MLR(ins1);
d13Cs2_MLR=d13C_MLR(ins2);
d13Cs3_MLR=d13C_MLR(ins3);
d13Cs4_MLR=d13C_MLR(ins4);
d13Cs5_MLR=d13C_MLR(ins5);
d13Cs6_MLR=d13C_MLR(ins6);
d13Cs7_MLR=d13C_MLR(ins7);
d13Cs8_MLR=d13C_MLR(ins8);


% sort by depth
[sdeps1,Is1]=sort(deps1);
spd13Cs1=pd13Cs4(Is1);
smd13Cs1=md13Cs4(Is1);
sd13Cs1_eMLR=d13Cs1_eMLR(Is1);
sd13Cs1_MLR=d13Cs1_MLR(Is1);
% std13Cs1_eMLR=stdd13Cs1_eMLR(Is1);
% std13Cs1_MLR=stdd13Cs1_MLR(Is1);

[sdeps2,Is2]=sort(deps2);
spd13Cs2=pd13Cs2(Is2);
smd13Cs2=md13Cs2(Is2);
sd13Cs2_eMLR=d13Cs2_eMLR(Is2);
sd13Cs2_MLR=d13Cs2_MLR(Is2);
% std13Cs2_eMLR=stdd13Cs2_eMLR(Is2);
% std13Cs2_MLR=stdd13Cs2_MLR(Is2);

[sdeps3,Is3]=sort(deps3);
spd13Cs3=pd13Cs3(Is3);
smd13Cs3=md13Cs3(Is3);
sd13Cs3_eMLR=d13Cs3_eMLR(Is3);
sd13Cs3_MLR=d13Cs3_MLR(Is3);
% std13Cs3_eMLR=stdd13Cs3_eMLR(Is3);
% std13Cs3_MLR=stdd13Cs3_MLR(Is3);

[sdeps4,Is4]=sort(deps4);
spd13Cs4=pd13Cs4(Is4);
smd13Cs4=md13Cs4(Is4);
sd13Cs4_eMLR=d13Cs4_eMLR(Is4);
sd13Cs4_MLR=d13Cs4_MLR(Is4);
% std13Cs4_eMLR=stdd13Cs4_eMLR(Is4);
% std13Cs4_MLR=stdd13Cs4_MLR(Is4);

[sdeps5,Is5]=sort(deps5);
spd13Cs5=pd13Cs5(Is5);
smd13Cs5=md13Cs5(Is5);
sd13Cs5_eMLR=d13Cs5_eMLR(Is5);
sd13Cs5_MLR=d13Cs5_MLR(Is5);
% std13Cs5_eMLR=stdd13Cs5_eMLR(Is5);
% std13Cs5_MLR=stdd13Cs5_MLR(Is5);

[sdeps6,Is6]=sort(deps6);
spd13Cs6=pd13Cs6(Is6);
smd13Cs6=md13Cs6(Is6);
sd13Cs6_eMLR=d13Cs6_eMLR(Is6);
sd13Cs6_MLR=d13Cs6_MLR(Is6);
% std13Cs6_eMLR=stdd13Cs6_eMLR(Is6);
% std13Cs6_MLR=stdd13Cs6_MLR(Is6);

[sdeps7,Is7]=sort(deps7);
spd13Cs7=pd13Cs7(Is7);
smd13Cs7=md13Cs7(Is7);
sd13Cs7_eMLR=d13Cs7_eMLR(Is7);
sd13Cs7_MLR=d13Cs7_MLR(Is7);
% std13Cs7_eMLR=stdd13Cs7_eMLR(Is7);
% std13Cs7_MLR=stdd13Cs7_MLR(Is7);

[sdeps8,Is8]=sort(deps8);
spd13Cs8=pd13Cs8(Is8);
smd13Cs8=md13Cs8(Is8);
sd13Cs8_eMLR=d13Cs8_eMLR(Is8);
sd13Cs8_MLR=d13Cs8_MLR(Is8);
% std13Cs8_eMLR=stdd13Cs8_eMLR(Is8);
% std13Cs8_MLR=stdd13Cs8_MLR(Is8);

%%
interp_depth=[0:50:6000];         

z2_eachsta(1,:)=jininterp1(sdeps1,sd13Cs1_eMLR,interp_depth);  
z2_eachsta(2,:)=jininterp1(sdeps2,sd13Cs2_eMLR,interp_depth);  
z2_eachsta(3,:)=jininterp1(sdeps3,sd13Cs3_eMLR,interp_depth);  
z2_eachsta(4,:)=jininterp1(sdeps4,sd13Cs4_eMLR,interp_depth);  
z2_eachsta(5,:)=jininterp1(sdeps5,sd13Cs5_eMLR,interp_depth);  
z2_eachsta(6,:)=jininterp1(sdeps6,sd13Cs6_eMLR,interp_depth);  
z2_eachsta(7,:)=jininterp1(sdeps7,sd13Cs7_eMLR,interp_depth);  
z2_eachsta(8,:)=jininterp1(sdeps8,sd13Cs8_eMLR,interp_depth);  

% inve=z2_eachsta.*(1000+den_eachsta);
inve=z2_eachsta;
inve(isnan(inve))=0;

for sta=1:8
inventory(sta,:)=trapz(interp_depth,inve(sta,:));
inventory3(sta,:)=trapz(interp_depth(1:21),squeeze(inve(sta,1:21))); % depth 0-1000 m
inventory2(sta,:)=trapz(interp_depth(1:41),squeeze(inve(sta,1:41))); % depth 0-2000 m
end

mean(inventory2)
std(inventory2)

invens1_eMLR=cumtrapz(sdeps1,sd13Cs1_eMLR);
invens2_eMLR=cumtrapz(sdeps2,sd13Cs2_eMLR);
invens3_eMLR=cumtrapz(sdeps3,sd13Cs3_eMLR);
invens4_eMLR=cumtrapz(sdeps4,sd13Cs4_eMLR);
invens5_eMLR=cumtrapz(sdeps5,sd13Cs5_eMLR);
invens6_eMLR=cumtrapz(sdeps6,sd13Cs6_eMLR);
invens7_eMLR=cumtrapz(sdeps7,sd13Cs7_eMLR);
invens8_eMLR=cumtrapz(sdeps8,sd13Cs8_eMLR);

invens1_MLR=cumtrapz(sdeps1,sd13Cs1_MLR);
invens2_MLR=cumtrapz(sdeps2,sd13Cs2_MLR);
invens3_MLR=cumtrapz(sdeps3,sd13Cs3_MLR);
invens4_MLR=cumtrapz(sdeps4,sd13Cs4_MLR);
invens5_MLR=cumtrapz(sdeps5,sd13Cs5_MLR);
invens6_MLR=cumtrapz(sdeps6,sd13Cs6_MLR);
invens7_MLR=cumtrapz(sdeps7,sd13Cs7_MLR);
invens8_MLR=cumtrapz(sdeps8,sd13Cs8_MLR);


% %%  figure
% 
%  eMLR Anthropogenic \delta^1^3C change (?)
figure
axes('position', [0.18 0.18 0.70 0.70])       
hold on
plot(d13Cs1_eMLR,deps1,'k*','MarkerSize',12);
plot(d13Cs2_eMLR,deps2,'r^','MarkerSize',7,'MarkerEdgeColor','r',...
    'MarkerFaceColor','r');
plot(d13Cs3_eMLR,deps3,'gs','MarkerSize',7,'MarkerEdgeColor','g',...
    'MarkerFaceColor','g');
plot(d13Cs4_eMLR,deps4,'bd','MarkerSize',7,'MarkerEdgeColor','b',...
    'MarkerFaceColor','b');
plot(d13Cs5_eMLR,deps5,'p','MarkerSize',7,'MarkerEdgeColor',[0.4940 0.1840 0.5560],...
    'MarkerFaceColor',[0.4940 0.1840 0.5560]);
plot(d13Cs6_eMLR,deps6,'co','MarkerSize',6,'MarkerEdgeColor','c',...
    'MarkerFaceColor','c');
plot(d13Cs7_eMLR,deps7,'mv','MarkerSize',6,'MarkerEdgeColor','m',...
    'MarkerFaceColor','m');
plot(d13Cs8_eMLR,deps8,'h','MarkerSize',6,'color',[0.8500 0.3250 0.0980],...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980],...
    'MarkerFaceColor',[0.8500 0.3250 0.0980]);
line([0,0],[0,6000],'Color','k');
ax = gca;ax.YDir= 'reverse'; 
legend('2020 s1','2020 s2',...
    '2020 s3','2020 s4',...
    '2020 s5','2020 s6',...
    '2020 s7','2020 s8','Location','southwest',...
    'NumColumns',1,'Fontsize',16,'Fontname','Times New Roman')
grid on
box on
xticks([-0.6 -0.4 -0.2 0 0.2])
set(gca,'XLim',[-0.6 0.2],'YLim',[0 6000],'Fontsize',18,'Fontname','Times New Roman')
xlabel('Anthropogenic \delta^1^3C change (¡ë)','Fontsize',20,'Fontname','Times New Roman');
ylabel('Depth (m)','Fontsize',20,'Fontname','Times New Roman');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 6])
% saveas(gcf,['./eMLR_Anthropogenic_d13C_change_A135_all.tif']);
%  print(gcf,'-dtiff','-r300','eMLR_anth_d13C_change_no_offset'); 
% 
% %  MLR Anthropogenic \delta^1^3C change (?)
% figure
% axes('position', [0.18 0.18 0.70 0.70])       
% hold on
% plot(d13Cs1_MLR,deps1,'k*','MarkerSize',12);
% plot(d13Cs2_MLR,deps2,'r^','MarkerSize',7,'MarkerEdgeColor','r',...
%     'MarkerFaceColor','r');
% plot(d13Cs3_MLR,deps3,'gs','MarkerSize',7,'MarkerEdgeColor','g',...
%     'MarkerFaceColor','g');
% plot(d13Cs4_MLR,deps4,'bd','MarkerSize',7,'MarkerEdgeColor','b',...
%     'MarkerFaceColor','b');
% plot(d13Cs5_MLR,deps5,'p','MarkerSize',7,'MarkerEdgeColor',[0.4940 0.1840 0.5560],...
%     'MarkerFaceColor',[0.4940 0.1840 0.5560]);
% plot(d13Cs6_MLR,deps6,'co','MarkerSize',6,'MarkerEdgeColor','c',...
%     'MarkerFaceColor','c');
% plot(d13Cs7_MLR,deps7,'mv','MarkerSize',6,'MarkerEdgeColor','m',...
%     'MarkerFaceColor','m');
% plot(d13Cs8_MLR,deps8,'h','MarkerSize',6,'color',[0.8500 0.3250 0.0980],...
%     'MarkerEdgeColor',[0.8500 0.3250 0.0980],...
%     'MarkerFaceColor',[0.8500 0.3250 0.0980]);
% line([0,0],[0,6000],'Color','k');
% ax = gca;ax.YDir= 'reverse'; 
% legend('2020 s1','2020 s2',...
%     '2020 s3','2020 s4',...
%     '2020 s5','2020 s6',...
%     '2020 s7','2020 s8','Location','southwest',...
%     'NumColumns',1,'Fontsize',16,'Fontname','Times New Roman')
% grid on
% box on
% xticks([-0.6 -0.4 -0.2 0 0.2])
% set(gca,'XLim',[-0.6 0.2],'YLim',[0 6000],'Fontsize',18,'Fontname','Times New Roman')
% xlabel('Anthropogenic \delta^1^3C change (¡ë)','Fontsize',20,'Fontname','Times New Roman');
% ylabel('Depth (m)','Fontsize',20,'Fontname','Times New Roman');
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 6])
% % saveas(gcf,['./eMLR_Anthropogenic_d13C_change_A135_all.tif']);
% print(gcf,'-dtiff','-r300','MLR_anth_d13C_change_no_offset');
% 
% %%  inventory
% % eMLR
% figure
% axes('position', [0.18 0.18 0.70 0.70])       
% hold on
% plot(invens1_eMLR,sdeps1,'k*','MarkerSize',12);
% plot(invens2_eMLR,sdeps2,'r^','MarkerSize',8,'MarkerEdgeColor','r',...
%     'MarkerFaceColor','r');
% plot(invens3_eMLR,sdeps3,'gs','MarkerSize',8,'MarkerEdgeColor','g',...
%     'MarkerFaceColor','g');
% plot(invens4_eMLR,sdeps4,'bd','MarkerSize',8,'MarkerEdgeColor','b',...
%     'MarkerFaceColor','b');
% plot(invens5_eMLR,sdeps5,'p','MarkerSize',8,'MarkerEdgeColor',[0.4940 0.1840 0.5560],...
%     'MarkerFaceColor',[0.4940 0.1840 0.5560]);
% plot(invens6_eMLR,sdeps6,'co','MarkerSize',8,'MarkerEdgeColor','c',...
%     'MarkerFaceColor','c');
% plot(invens7_eMLR,sdeps7,'mv','MarkerSize',8,'MarkerEdgeColor','m',...
%     'MarkerFaceColor','m');
% plot(invens8_eMLR,sdeps8,'h','MarkerSize',8,'color',[0.8500 0.3250 0.0980],...
%     'MarkerEdgeColor',[0.8500 0.3250 0.0980],...
%     'MarkerFaceColor',[0.8500 0.3250 0.0980]);
% ax = gca;ax.YDir= 'reverse'; 
% legend('station 1','station 2',...
%     'station 3','station 4',...
%     'station 5','station 6',...
%     'station 7','station 8','Location','southeast',...
%     'NumColumns',2,'Fontsize',16,'Fontname','Times New Roman')
% grid on
% box on
% set(gca,'XLim',[-300 50],'YLim',[0 6000],'Fontsize',22,'Fontname','Times New Roman')
% % set(gca,'XLim',[-600 0],'YLim',[0 6000],'Fontsize',22,'Fontname','Times New Roman')
% xlabel('Anthropogenic \delta^1^3C inventory change (¡ë)','Fontsize',24,'Fontname','Times New Roman');
% ylabel('Depth (m)','Fontsize',24,'Fontname','Times New Roman');
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
% print(gcf,'-dtiff','-r300','eMLR_d13C_inventory_change_offset');
% 
% % MLR
% figure
% axes('position', [0.18 0.18 0.70 0.70])       
% hold on
% plot(invens1_MLR,sdeps1,'k*','MarkerSize',12);
% plot(invens2_MLR,sdeps2,'r^','MarkerSize',8,'MarkerEdgeColor','r',...
%     'MarkerFaceColor','r');
% plot(invens3_MLR,sdeps3,'gs','MarkerSize',8,'MarkerEdgeColor','g',...
%     'MarkerFaceColor','g');
% plot(invens4_MLR,sdeps4,'bd','MarkerSize',8,'MarkerEdgeColor','b',...
%     'MarkerFaceColor','b');
% plot(invens5_MLR,sdeps5,'p','MarkerSize',8,'MarkerEdgeColor',[0.4940 0.1840 0.5560],...
%     'MarkerFaceColor',[0.4940 0.1840 0.5560]);
% plot(invens6_MLR,sdeps6,'co','MarkerSize',8,'MarkerEdgeColor','c',...
%     'MarkerFaceColor','c');
% plot(invens7_MLR,sdeps7,'mv','MarkerSize',8,'MarkerEdgeColor','m',...
%     'MarkerFaceColor','m');
% plot(invens8_MLR,sdeps8,'h','MarkerSize',8,'color',[0.8500 0.3250 0.0980],...
%     'MarkerEdgeColor',[0.8500 0.3250 0.0980],...
%     'MarkerFaceColor',[0.8500 0.3250 0.0980]);
% ax = gca;ax.YDir= 'reverse'; 
% legend('station 1','station 2',...
%     'station 3','station 4',...
%     'station 5','station 6',...
%     'station 7','station 8','Location','southeast',...
%     'NumColumns',2,'Fontsize',16,'Fontname','Times New Roman')
% grid on
% box on
% set(gca,'XLim',[-300 50],'YLim',[0 6000],'Fontsize',22,'Fontname','Times New Roman')
% % set(gca,'XLim',[-600 0],'YLim',[0 6000],'Fontsize',22,'Fontname','Times New Roman')
% xlabel('Anthropogenic \delta^1^3C inventory change (¡ë)','Fontsize',24,'Fontname','Times New Roman');
% ylabel('Depth (m)','Fontsize',24,'Fontname','Times New Roman');
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
% print(gcf,'-dtiff','-r300','MLR_d13C_inventory_change_offset');
% 
% % %% std
% % %  eMLR Anthropogenic \delta^1^3C change (?)
% % figure
% % axes('position', [0.18 0.18 0.70 0.70])       
% % hold on
% % plot(std13Cs1_eMLR,deps1,'k*','MarkerSize',12);
% % plot(std13Cs2_eMLR,deps2,'r^','MarkerSize',7,'MarkerEdgeColor','r',...
% %     'MarkerFaceColor','r');
% % plot(std13Cs3_eMLR,deps3,'gs','MarkerSize',7,'MarkerEdgeColor','g',...
% %     'MarkerFaceColor','g');
% % plot(std13Cs4_eMLR,deps4,'bd','MarkerSize',7,'MarkerEdgeColor','b',...
% %     'MarkerFaceColor','b');
% % plot(std13Cs5_eMLR,deps5,'p','MarkerSize',7,'MarkerEdgeColor',[0.4940 0.1840 0.5560],...
% %     'MarkerFaceColor',[0.4940 0.1840 0.5560]);
% % plot(std13Cs6_eMLR,deps6,'co','MarkerSize',6,'MarkerEdgeColor','c',...
% %     'MarkerFaceColor','c');
% % plot(std13Cs7_eMLR,deps7,'mv','MarkerSize',6,'MarkerEdgeColor','m',...
% %     'MarkerFaceColor','m');
% % plot(std13Cs8_eMLR,deps8,'h','MarkerSize',6,'color',[0.8500 0.3250 0.0980],...
% %     'MarkerEdgeColor',[0.8500 0.3250 0.0980],...
% %     'MarkerFaceColor',[0.8500 0.3250 0.0980]);
% % line([0,0],[0,6000],'Color','k');
% % ax = gca;ax.YDir= 'reverse'; 
% % legend('2020 s1','2020 s2',...
% %     '2020 s3','2020 s4',...
% %     '2020 s5','2020 s6',...
% %     '2020 s7','2020 s8','Location','southwest',...
% %     'NumColumns',1,'Fontsize',16,'Fontname','Times New Roman')
% % grid on
% % box on
% % % xticks([-0.6 -0.4 -0.2 0 0.2])
% % % set(gca,'XLim',[-0.6 0.2],'YLim',[0 6000],'Fontsize',18,'Fontname','Times New Roman')
% % xlabel('Anthropogenic \delta^1^3C change (¡ë)','Fontsize',20,'Fontname','Times New Roman');
% % ylabel('Depth (m)','Fontsize',20,'Fontname','Times New Roman');
% % set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 6])
% % % saveas(gcf,['./eMLR_Anthropogenic_d13C_change_A135_all.tif']);
% % print(gcf,'-dtiff','-r300','eMLR_anth_d13C_change_reg1_7_std'); 
% % 
% % %  MLR Anthropogenic \delta^1^3C change (?)
% % figure
% % axes('position', [0.18 0.18 0.70 0.70])       
% % hold on
% % plot(std13Cs1_MLR,deps1,'k*','MarkerSize',12);
% % plot(std13Cs2_MLR,deps2,'r^','MarkerSize',7,'MarkerEdgeColor','r',...
% %     'MarkerFaceColor','r');
% % plot(std13Cs3_MLR,deps3,'gs','MarkerSize',7,'MarkerEdgeColor','g',...
% %     'MarkerFaceColor','g');
% % plot(std13Cs4_MLR,deps4,'bd','MarkerSize',7,'MarkerEdgeColor','b',...
% %     'MarkerFaceColor','b');
% % plot(std13Cs5_MLR,deps5,'p','MarkerSize',7,'MarkerEdgeColor',[0.4940 0.1840 0.5560],...
% %     'MarkerFaceColor',[0.4940 0.1840 0.5560]);
% % plot(std13Cs6_MLR,deps6,'co','MarkerSize',6,'MarkerEdgeColor','c',...
% %     'MarkerFaceColor','c');
% % plot(std13Cs7_MLR,deps7,'mv','MarkerSize',6,'MarkerEdgeColor','m',...
% %     'MarkerFaceColor','m');
% % plot(std13Cs8_MLR,deps8,'h','MarkerSize',6,'color',[0.8500 0.3250 0.0980],...
% %     'MarkerEdgeColor',[0.8500 0.3250 0.0980],...
% %     'MarkerFaceColor',[0.8500 0.3250 0.0980]);
% % line([0,0],[0,6000],'Color','k');
% % ax = gca;ax.YDir= 'reverse'; 
% % legend('2020 s1','2020 s2',...
% %     '2020 s3','2020 s4',...
% %     '2020 s5','2020 s6',...
% %     '2020 s7','2020 s8','Location','southwest',...
% %     'NumColumns',1,'Fontsize',16,'Fontname','Times New Roman')
% % grid on
% % box on
% % % xticks([-0.6 -0.4 -0.2 0 0.2])
% % % set(gca,'XLim',[-0.6 0.2],'YLim',[0 6000],'Fontsize',18,'Fontname','Times New Roman')
% % xlabel('Anthropogenic \delta^1^3C change (¡ë)','Fontsize',20,'Fontname','Times New Roman');
% % ylabel('Depth (m)','Fontsize',20,'Fontname','Times New Roman');
% % set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5.5 6])
% % % saveas(gcf,['./eMLR_Anthropogenic_d13C_change_A135_all.tif']);
% % print(gcf,'-dtiff','-r300','MLR_anth_d13C_change_reg1_7_std');
% 
% 
