close all;clear;clc;
%%
% load 2010 data
load('A13.5_2010.mat')

% d13C
in=find(temperature~=-9999 & salinity~=-9999 & aou~=-9999 & nitrate~=-9999 ...
       & c13~=-9999 & silicate~=-9999 & latitudedegrees_north>=-42 & latitudedegrees_north<=-32);  
x(:,1)=temperature(in);
x(:,2)=salinity(in);
x(:,3)=aou(in);
x(:,4)=nitrate(in);
x(:,5)=silicate(in);
dep=depthm(in);
% y=tco2(in);
y=c13(in);
den=gamma1(in);

density=gamma1(in);
in1=find(density<26.8);
in2=find(density>=26.8 & density<27.23);
in3=find(density>=27.23 & density<27.5);
in4=find(density>=27.5 & density<28);
% in5=find(density>=28 & density<28.27);
% in6=find(density>=28.27);
in5=find(density>=28);

x2010{1}=x(in1,1:5);y2010{1}=y(in1);
x2010{2}=x(in2,1:5);y2010{2}=y(in2);
x2010{3}=x(in3,1:5);y2010{3}=y(in3);
x2010{4}=x(in4,1:5);y2010{4}=y(in4);
x2010{5}=x(in5,1:5);y2010{5}=y(in5);
% x2010{6}=x(in6,1:5);y2010{6}=y(in6);

sd=[25.8 27 27.35 27.75 28.15];

for i=1:5
md13C_2010(i)=nanmean(y2010{i});
stdd13C_2010(i)=nanstd(y2010{i});
end



%% load 2020 data
clear in in1 in2 in3 in4 in5 density
load('A13.5_2020v5.mat')

% MLR fit equation DIC=a+b*theta+c*salinity+d*AOU+e*silicate+f*phosphate
in=find(CTDTMP~=-999 & CTDSAL~=-999 & AOU~=-999 & NITRAT~=-999 & SILCAT~=-999 & d13C_lab_Najid~=-999); 
xx(:,1)=CTDTMP(in);
xx(:,2)=CTDSAL(in);
xx(:,3)=AOU4_S(in);
xx(:,4)=NITRAT(in);
xx(:,5)=SILCAT(in);
pres=CTDPRS(in);
dept=cal_depth(in);
% yy=d13C_lab_Najid(in);
% yy=TCARBN(in);
yy=d13C_lab_Najid(in); 
dens=cal_gamma(in);

density=cal_gamma(in);
in1=find(density<26.8);
in2=find(density>=26.8 & density<27.23);
in3=find(density>=27.23 & density<27.5);
in4=find(density>=27.5 & density<28);
% in5=find(density>=28 & density<28.27);
% in6=find(density>=28.27);
in5=find(density>=28);

x2020{1}=xx(in1,1:5);y2020{1}=yy(in1);
x2020{2}=xx(in2,1:5);y2020{2}=yy(in2);
x2020{3}=xx(in3,1:5);y2020{3}=yy(in3);
x2020{4}=xx(in4,1:5);y2020{4}=yy(in4);
x2020{5}=xx(in5,1:5);y2020{5}=yy(in5);
% x2020{6}=x(in6,1:5);y2020{6}=y(in6);

for i=1:5
md13C_2020(i)=nanmean(y2020{i});
stdd13C_2020(i)=nanstd(y2020{i});
end


%%
figure
sz=60;

top_margin = 0.08; % top margin
btm_margin = 0.18; % bottom margin
left_margin = 0.19;% left margin
right_margin = 0.05;% right margin
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
line([0.5,0.5],[24,28.5],'Color',[0.9 0.9 0.9]);
line([1,1],[24,28.5],'Color',[0.9 0.9 0.9]);
line([1.5,1.5],[24,28.5],'Color',[0.9 0.9 0.9]);
line([2,2],[24,28.5],'Color',[0.9 0.9 0.9]);

scatter(y,den,sz,'r^','MarkerEdgeAlpha',0.5);
scatter(yy,dens,sz,'k*','MarkerEdgeAlpha',0.5);
errorbar(md13C_2010,sd,stdd13C_2010,'horizontal','r^','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','r','linewidth',2);
errorbar(md13C_2020,sd,stdd13C_2020,'horizontal','ko','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',2);


set(gca,'Fontsize',16,'XLim',[0.45 2.2],'YLim',[24.8 26.8],'Fontname','Times New Roman')
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
scatter(y,den,sz,'r^','MarkerEdgeAlpha',0.5);
scatter(yy,dens,sz,'k*','MarkerEdgeAlpha',0.5);
errorbar(md13C_2010,sd,stdd13C_2010,'horizontal','r^','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','r','linewidth',2);
errorbar(md13C_2020,sd,stdd13C_2020,'horizontal','ko','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',2);

line([0,0],[24,28.5],'Color',[0.8 0.8 0.8]);

set(gca,'Fontsize',16,'XLim',[0.45 2.2],'YLim',[26.8 28.3],'Fontname','Times New Roman')
% xticks([-0.6 -0.4 -0.2 0 0.2])
yticks([27.23 27.5 28 28.27])
box on
grid on
% hold off
ax = gca;ax.YDir= 'reverse'; 
% legend('Measured DIC in 2010','Measured DIC in 2020',...
%    'Mean & Standard deviation in 2010','Mean & Standard deviation in 2020',...
%    'Location','southwest',...
%     'NumColumns',1,'Fontsize',13.5,'Fontname','Times New Roman')
legend('2010','2020',...
   '','',...
   'Location','southeast',...
    'NumColumns',1,'Fontsize',13.5,'Fontname','Times New Roman')


xlabel('\delta^1^3C (‰)','Fontsize',18,'Fontname','Times New Roman');
% ylabel('Depth (m)','Fontsize',26,'Fontname','Times New Roman');
[ax,h2]=suplabel('Neutral density (kg m^{-3})','y');
set(h2,'Fontsize',18,'Fontname','Times New Roman','position',[-0.1 0.5])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 6])
% print(gcf,'-dtiff','-r300','F2d_d13C_density.tif'); 

