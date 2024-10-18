close all;clear;clc;
%%
% load 2010 data
load('A13.5_2010.mat')

% d13C
in=find(temperature~=-9999 & salinity~=-9999 & aou~=-9999 & nitrate~=-9999 ...
       & tco2~=-9999 & silicate~=-9999 & latitudedegrees_north>=-42 & latitudedegrees_north<=-32);  
x(:,1)=temperature(in);
x(:,2)=salinity(in);
x(:,3)=aou(in);
x(:,4)=nitrate(in);
x(:,5)=silicate(in);
dep=depthm(in);
y=tco2(in);

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

interp_depth=[10 50 100 150 200 250 300 400 500 600 700 800 900 1000 1100 1200 ...
    1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400 3600 3800 4000 4200 ...
    4400 4600 4800 4900];

sn=station(in);
nn=1;
for n=min(sn):max(sn)
ins{nn}=find(sn==n);
ds=dep(ins{nn});

[l,h]=size(ds);
for i=2:l
    if ds(i-1)==ds(i)
        ds(i)=ds(i) +0.00001;      
    end
end

dsn{nn}=ds;

ysn{nn}=y(ins{nn});
ysin(:,nn)=interp1(dsn{nn},ysn{nn},interp_depth);
nn=nn+1;
end


for i=1:35
my1(i)=nanmean(ysin(i,:));
stdy1(i)=nanstd(ysin(i,:));
end

sd1=interp_depth;
sd1(isnan(my1))=[];
my1(isnan(my1))=[];
stdy1(isnan(stdy1))=[];


x1_low=my1-stdy1;
x1_up=flip(my1+stdy1);
x1_er = [x1_low,x1_up,x1_low(1)];
y1_er = [sd1,flip(sd1),sd1(1)];


%% load 2020 data

load('A13.5_2020v5.mat')

% MLR fit equation DIC=a+b*theta+c*salinity+d*AOU+e*silicate+f*phosphate
in=find(CTDTMP~=-999 & CTDSAL~=-999 & AOU~=-999 & NITRAT~=-999 & SILCAT~=-999 & TCARBN~=-999); 
xx(:,1)=CTDTMP(in);
xx(:,2)=CTDSAL(in);
xx(:,3)=AOU4_S(in);
xx(:,4)=NITRAT(in);
xx(:,5)=SILCAT(in);
pres=CTDPRS(in);
dept=cal_depth(in);
% yy=d13C_lab_Najid(in);
yy=TCARBN(in);

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

ssn=STNNBR(in);
for n=1:8
inssn{n}=find(ssn==n);
dssn{n}=dept(inssn{n});
yssn{n}=yy(inssn{n});
yysin(:,n)=interp1(dssn{n},yssn{n},interp_depth);

end

for i=1:35
my2(i)=nanmean(yysin(i,:));
stdy2(i)=nanstd(yysin(i,:));
end

sd2=interp_depth;
sd2(isnan(my2))=[];
my2(isnan(my2))=[];
stdy2(isnan(stdy2))=[];

x2_low=my2-stdy2;
x2_up=flip(my2+stdy2);
x2_er = [x2_low,x2_up,x2_low(1)];
y2_er = [sd2,flip(sd2),sd2(1)];


%%  anthropogenic d13C change
for slab=1:5   
    x11=x2010{1,slab};
    y11=y2010{1,slab};
    x22=x2020{1,slab};
    y22=y2020{1,slab};

   [b1,stats1] = robustfit(x11,y11);  
   coe1(slab,:)=b1'; 
   [b2,stats2] = robustfit(x22,y22);  
   coe2(slab,:)=b2'; 
   
%    resi1=stats1.resid;
%    resi2=stats2.resid;
%    rmse1(slab)=stats1.ols_s;
%    rmse2(slab)=stats2.ols_s;
%    mresi1(slab)=nanmean(resi1);
%    stdresi1(slab)=nanstd(resi1);
%    mresi2(slab)=nanmean(resi2);
%    stdresi2(slab)=nanstd(resi2);
   
    pd13C1{slab}=coe1(slab,1)+sum(coe1(slab,2:6).*x11,2);
    pd13C2{slab}=coe2(slab,1)+sum(coe2(slab,2:6).*x22,2);

    [rho1(:,slab),pval1(:,slab)] = corr(y11,pd13C1{slab});
    [rho2(:,slab),pval2(:,slab)] = corr(y22,pd13C2{slab});
   
   coediff(slab,:)=b2'-b1';
   preddelC{slab}=coe1(slab,1)+sum(coe1(slab,2:6).*x22,2);  % MLR predicted 
   chd13C{slab}=y22-preddelC{slab};   % MLR method   % measured-predicted
   delC{slab}=coediff(slab,1)+sum(coediff(slab,2:6).*x22,2);  % eMLR method

      clear  x11 x22 y11 y22
end

%  predicted d13C 2020
pd1=preddelC{1};pd2=preddelC{2};pd3=preddelC{3};pd4=preddelC{4};pd5=preddelC{5};
pd13C(1:length(pd1))=pd1;
pd13C(1+length(pd1):length(pd1)+length(pd2))=pd2;
pd13C(1+length(pd1)+length(pd2):length(pd1)+length(pd2)+length(pd3))=pd3;
pd13C(1+length(pd1)+length(pd2)+length(pd3):length(pd1)+length(pd2)+length(pd3)+length(pd4))=pd4;
pd13C(1+length(pd1)+length(pd2)+length(pd3)+length(pd4):length(pd1)+length(pd2)+length(pd3)+length(pd4)+length(pd5))=pd5;

%  depth
sdep1=dept(in1);sdep2=dept(in2);sdep3=dept(in3);
sdep4=dept(in4);sdep5=dept(in5);
sdep(1:length(pd1))=sdep1;
sdep(1+length(pd1):length(pd1)+length(pd2))=sdep2;
sdep(1+length(pd1)+length(pd2):length(pd1)+length(pd2)+length(pd3))=sdep3;
sdep(1+length(pd1)+length(pd2)+length(pd3):length(pd1)+length(pd2)+length(pd3)+length(pd4))=sdep4;
sdep(1+length(pd1)+length(pd2)+length(pd3)+length(pd4):length(pd1)+length(pd2)+length(pd3)+length(pd4)+length(pd5))=sdep5;

% arrange by station number

ssn1=ssn(in1);ssn2=ssn(in2);ssn3=ssn(in3);
ssn4=ssn(in4);ssn5=ssn(in5);

stn(1:length(pd1))=ssn1;
stn(1+length(pd1):length(pd1)+length(pd2))=ssn2;
stn(1+length(pd1)+length(pd2):length(pd1)+length(pd2)+length(pd3))=ssn3;
stn(1+length(pd1)+length(pd2)+length(pd3):length(pd1)+length(pd2)+length(pd3)+length(pd4))=ssn4;
stn(1+length(pd1)+length(pd2)+length(pd3)+length(pd4):length(pd1)+length(pd2)+length(pd3)+length(pd4)+length(pd5))=ssn5;

for n=1:8
instn{n}=find(stn==n);
dstn{n}=sdep(instn{n});
ystn{n}=pd13C(instn{n});
yyysin(:,n)=interp1(dstn{n},ystn{n},interp_depth);

end

for i=1:35
my3(i)=nanmean(yyysin(i,:));
stdy3(i)=nanstd(yyysin(i,:));
end

sd3=interp_depth;
sd3(isnan(my3))=[];
my3(isnan(my3))=[];
stdy3(isnan(stdy3))=[];

x3_low=my3-stdy3;
x3_up=flip(my3+stdy3);
x3_er = [x3_low,x3_up,x3_low(1)];
y3_er = [sd3,flip(sd3),sd3(1)];


%%
figure

sz=40;
top_margin = 0.08; % top margin
btm_margin = 0.18; % bottom margin
left_margin = 0.16;% left margin
right_margin = 0.08;% right margin
fig_margin = 0.001; % margin beween figures(sub)  
row = 2; % rows
col = 1; % cols 

% Calculate figure height and width according to rows and cols 
fig_h = (1- top_margin - btm_margin - (row-1) * fig_margin) / row + 0.04;
fig_w = (1 - left_margin - right_margin - (col-1) * fig_margin) / col;

% figure position: you can refer to 'help axes' to review the        
% parameter meaning, note that original point is lower-left        
position = [left_margin + (1-1)*(fig_margin+fig_w), ...           
    1- (top_margin + 1 * fig_h + (1-1) * fig_margin) + 0.02, ...           
    fig_w, fig_h];       
axes('position', position)       
% draw colorful pictures...    

hold on;
line([2050,2050],[0 2000],'Color',[0.9 0.9 0.9]);
line([2100,2100],[0 2000],'Color',[0.9 0.9 0.9]);
line([2150,2150],[0 2000],'Color',[0.9 0.9 0.9]);
line([2200,2200],[0 2000],'Color',[0.9 0.9 0.9]);
line([2250,2250],[0 2000],'Color',[0.9 0.9 0.9]);
line([2000,2300],[500 500],'Color',[0.9 0.9 0.9]);
line([2000,2300],[1000 1000],'Color',[0.9 0.9 0.9]);
line([2000,2300],[1500 1500],'Color',[0.9 0.9 0.9]);


% plot(y,dep,'r^','MarkerSize',8);
% plot(yy,dept,'k*','MarkerSize',10);
% plot(pd13C,sdep,'bo','MarkerSize',8);
scatter(y,dep,sz,'r^','MarkerEdgeAlpha',0.5);
scatter(yy,dept,sz,'k*','MarkerEdgeAlpha',0.5);
% scatter(pd13C,sdep,sz,'bo','MarkerEdgeAlpha',0.5);

% p1=errorbar(my1,sd,stdy1,'horizontal','rs-','MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','r');
% p2=errorbar(my2,sd,stdy2,'horizontal','ks','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
% p3=errorbar(my3,sd,stdy3,'horizontal','bs-','MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor','b');
plot(my1,sd1,'r-','linewidth',1.5)
h1=fill(x1_er,y1_er,'r');
h1.FaceColor=[1 0 0];
h1.FaceAlpha=0.15;
% h1.EdgeColor=[1 1 1];
h1.EdgeAlpha=0;

plot(my2,sd2,'k-','linewidth',1.5)
h2=fill(x2_er,y2_er,'k');
h2.FaceColor=[0 0 0];
h2.FaceAlpha=0.15;
% h2.EdgeColor=[1 1 1];
h2.EdgeAlpha=0;

% plot(my3,sd3,'b-','linewidth',1.5)
% h3=fill(x3_er,y3_er,'b');
% h3.FaceColor=[0 0 1];
% h3.FaceAlpha=0.15;
% % h3.EdgeColor=[1 1 1];
% h3.EdgeAlpha=0;


set(gca,'Fontsize',16,'XLim',[2040 2260],'YLim',[0 2000],'Fontname','Times New Roman')
xlabel([ ])
xticks([ ])
box on
ax = gca;ax.YDir= 'reverse'; 


position = [left_margin + (1-1)*(fig_margin+fig_w), ...           
    1- (top_margin + 2 * fig_h + (2-1) * fig_margin) + 0.02, ...           
    fig_w, fig_h];       
axes('position', position)
hold on
% plot(y,dep,'r^','MarkerSize',8);
% plot(yy,dept,'k*','MarkerSize',10);
% plot(pd13C,sdep,'bo','MarkerSize',8);
scatter(y,dep,sz,'r^','MarkerEdgeAlpha',0.5);
scatter(yy,dept,sz,'k*','MarkerEdgeAlpha',0.5);
% scatter(pd13C,sdep,sz,'bo','MarkerEdgeAlpha',0.5);

% p4=errorbar(my1,sd,stdy1,'horizontal','rs-','MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','r');
% p5=errorbar(my2,sd,stdy2,'horizontal','ks','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');
% p6=errorbar(my3,sd,stdy3,'horizontal','bs-','MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor','b');
plot(my1,sd1,'r-','linewidth',1.5)
h1=fill(x1_er,y1_er,'r');
h1.FaceColor=[1 0 0];
h1.FaceAlpha=0.15;
% h1.EdgeColor=[1 1 1];
h1.EdgeAlpha=0;

plot(my2,sd2,'k-','linewidth',1.5)
h2=fill(x2_er,y2_er,'k');
h2.FaceColor=[0 0 0];
h2.FaceAlpha=0.15;
% h2.EdgeColor=[1 1 1];
h2.EdgeAlpha=0;

% plot(my3,sd3,'b-','linewidth',1.5)
% h3=fill(x3_er,y3_er,'b');
% h3.FaceColor=[0 0 1];
% h3.FaceAlpha=0.15;
% % h3.EdgeColor=[1 1 1];
% h3.EdgeAlpha=0;


set(gca,'Fontsize',16,'XLim',[2040 2260],'YLim',[2000 6000],'Fontname','Times New Roman')
xlabel('Latitude','Fontsize',18,'Fontname','Times New Roman');
yticks([3000 4000 5000 6000])
box on
grid on
% hold off
ax = gca;ax.YDir= 'reverse'; 
% legend('2010 measured','2020 measured','2020 MLR predicted', ...
%    '2010 meas average','2010 meas STD',...
%    '2020 meas average','2020 meas STD',...
%    '2020 MLR pred average','2020 MLR pred STD', ...
%    'Location','southwest',...
%     'NumColumns',1,'Fontsize',13.5,'Fontname','Times New Roman')
legend('Measured DIC in 2010','Measured DIC in 2020',...
   'Mean value in 2010','Standard deviation in 2010',...
   'Mean value in 2020','Standard deviation in 2020',...
   'Location','southwest',...
    'NumColumns',1,'Fontsize',13.5,'Fontname','Times New Roman')
% set(gca,'XLim',[0.4 1.8],'YLim',[0 6000],'Fontsize',24,'Fontname','Times New Roman')
xlabel('DIC (\mumol kg^{-1})','Fontsize',18,'Fontname','Times New Roman');
% ylabel('Depth (m)','Fontsize',26,'Fontname','Times New Roman');
[ax,h2]=suplabel('Depth (m)','y');
set(h2,'Fontsize',18,'Fontname','Times New Roman','position',[-0.08 0.5])

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 6])
print(gcf,'-dtiff','-r300','F2a_compare_DIC_M_mean.tif'); 

