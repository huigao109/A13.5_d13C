close all;clear;clc;
% load data
% load 2010 data
load('A13.5_2010.mat')

% in=find(temperature~=-9999  & salinity~=-9999 & gamma1~=-9999 & aou~=-9999 ...
%         & nitrate~=-9999 & silicate~=-9999 & phosphate~=-9999 & c13~=-9999 ...
%        & latitudedegrees_north>=-42 & latitudedegrees_north<=-32 & depthm>=50); 

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

% in=find(CTDTMP~=-999 & CTDSAL~=-999 & AOU~=-999 & NITRAT~=-999 & SILCAT~=-999 ...
%         & PHSPHT~=-999  & DELC13~=-999 & cal_depth>=50); 
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

% y=d13C_lab_Najid(in)+0.07;   % offset
% y=DELC13(in);  
y=DELC13(in)+0.07;   % offset
ssn=STNNBR(in);
lats=LATITUDE(in);


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

%  depth
sdep1=depth(in1);sdep2=depth(in2);sdep3=depth(in3);
sdep4=depth(in4);sdep5=depth(in5);
sdep(1:length(d1))=sdep1;
sdep(1+length(d1):length(d1)+length(d2))=sdep2;
sdep(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=sdep3;
sdep(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=sdep4;
sdep(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=sdep5;

%  denisty
sden1=density(in1);sden2=density(in2);sden3=density(in3);
sden4=density(in4);sden5=density(in5);
sden(1:length(d1))=sden1;
sden(1+length(d1):length(d1)+length(d2))=sden2;
sden(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=sden3;
sden(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=sden4;
sden(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=sden5;

%  latitude
slat1=lats(in1);slat2=lats(in2);slat3=lats(in3);
slat4=lats(in4);slat5=lats(in5);
slat(1:length(d1))=slat1;
slat(1+length(d1):length(d1)+length(d2))=slat2;
slat(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=slat3;
slat(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=slat4;
slat(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=slat5;

% arrange by station number
ssn1=ssn(in1);ssn2=ssn(in2);ssn3=ssn(in3);
ssn4=ssn(in4);ssn5=ssn(in5);
stn(1:length(d1))=ssn1;
stn(1+length(d1):length(d1)+length(d2))=ssn2;
stn(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=ssn3;
stn(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=ssn4;
stn(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=ssn5;


%%%% interpolate into standard depth at each station
XI=[-42:0.5:-41 -37:0.5:-32]; 
YI=[0:50:6000]; 

[stns,Is]=sort(stn);
sdeps=sdep(Is);
sdens=sden(Is);
slats=slat(Is);
pd13Canths=d13C_eMLR(Is);

minsta=min(stns);maxsta=max(stns);
    sta=minsta;
% for i=1:maxsta-minsta+1
i=0;
for st=[2 1 3 4 5 6 7 8]
    i=i+1;
    inn=find(stns==st);
deps=sdeps(inn);
dens=sdens(inn);

pd13Canthse=pd13Canths(inn);
slatss(i)=slats(inn(1));
[sdepss,Isd]=sort(deps);
data=pd13Canthse(Isd);
sdenss=dens(Isd);
sdata=jininterp1(sdepss,data,YI); 
ssdenss=jininterp1(sdepss,sdenss,YI); 
z_eachsta(i,:)=sdata; 
z_eachden(i,:)=ssdenss;
%     sta=sta+1;
end

for k=1:length(YI)    
data=z_eachsta(:,k);
eachden=z_eachden(:,k);
sdata=jininterp1(slatss,data,XI);  
s_den=jininterp1(slatss,eachden,XI); 
szdata(:,k)=sdata; 
szden(:,k)=s_den;

end

inn1=find(szden<26.8);
inn2=find(szden>=26.8 & szden<27.23);
inn3=find(szden>=27.23 & szden<27.5);
inn4=find(szden>=27.5 & szden<28);
inn5=find(szden>=28);

sd=[25.8 27 27.35 27.75 28.15];

yy{1}=szdata(inn1);yy{2}=szdata(inn2);yy{3}=szdata(inn3);
yy{4}=szdata(inn4);yy{5}=szdata(inn5);

for i=1:5
md13C_eMLR(i)=nanmean(yy{i});
stdd13C_eMLR(i)=nanstd(yy{i});
end


figure
sz=60;

top_margin = 0.08; % top margin
btm_margin = 0.20; % bottom margin
left_margin = 0.20;% left margin
right_margin = 0.08;% right margin
fig_margin = 0.001; % margin beween figures(sub)  
row = 2; % rows
col = 1; % cols 

% Calculate figure height and width according to rows and cols 
fig_h = (1- top_margin - btm_margin - (row-1) * fig_margin)/row;
fig_w = (1 - left_margin - right_margin - (col-1) * fig_margin) / col;

% figure position: you can refer to 'help axes' to review the        
% parameter meaning, note that original point is lower-left        
position = [left_margin + (1-1)*(fig_margin+fig_w), ...           
    1- (top_margin + 1 * fig_h + (1-1) * fig_margin)+fig_h/2.4, ...           
    fig_w, fig_h/8*5];       
axes('position', position)       
% draw colorful pictures...    
hold on;
scatter(d13C_eMLR,sden,sz,'ks');
errorbar(md13C_eMLR,sd,stdd13C_eMLR,'horizontal','ks','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',2);
line([0,0],[24,28.5],'Color',[0.5 0.5 0.5]);
line([-0.2,-0.2],[24,28.5],'Color',[0.8 0.8 0.8]);
line([-0.4,-0.4],[24,28.5],'Color',[0.8 0.8 0.8]);

set(gca,'Fontsize',20,'XLim',[-0.6 0.2],'YLim',[24.8 26.8],'Fontname','Times New Roman')
xticks([ ])
yticks([26 26.8])
box on
grid on
ax = gca;ax.YDir= 'reverse'; 


position = [left_margin + (1-1)*(fig_margin+fig_w), ...           
    1- (top_margin + 2 * fig_h + (2-1) * fig_margin)-fig_h/10, ...           
    fig_w, fig_h/8*12];       
axes('position', position)
hold on
scatter(d13C_eMLR,sden,sz,'ks');
errorbar(md13C_eMLR,sd,stdd13C_eMLR,'horizontal','ks','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',2);
line([0,0],[24,28.5],'Color',[0.5 0.5 0.5]);

set(gca,'Fontsize',20,'XLim',[-0.6 0.2],'YLim',[26.8 28.3],'Fontname','Times New Roman')
xticks([-0.6 -0.4 -0.2 0 0.2])
yticks([27.23 27.5 28 28.27])
box on
grid on
% hold off
ax = gca;ax.YDir= 'reverse'; 
% legend('2010 measured','2020 measured','2020 MLR-predicted', ...
%    '2010 measured average','2010 measured STD',...
%    '2020 measured average','2020 measured STD',...
%    '2020 MLR-predicted average','2020 MLR-predicted STD', ...
%    'Location','southeast',...
%     'NumColumns',1,'Fontsize',16,'Fontname','Times New Roman')
xlabel('Anthropogenic \delta^1^3C change (‰)','Fontsize',22,'Fontname','Times New Roman');
% ylabel('Depth (m)','Fontsize',26,'Fontname','Times New Roman');
[ax,h2]=suplabel('Neutral density (kg m^{-3})','y');
set(h2,'Fontsize',22,'Fontname','Times New Roman','position',[-0.1 0.5])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6])
% print(gcf,'-dtiff','-r300','eMLR_d13C_change_density_F3d.tif'); 

