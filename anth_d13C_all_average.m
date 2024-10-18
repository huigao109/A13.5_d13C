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
pd13Canthse=pd13Canths(inn);
slatss(i)=slats(inn(1));
[sdepss,Isd]=sort(deps);
data=pd13Canthse(Isd);
sdata=jininterp1(sdepss,data,YI);    
z_eachsta(i,:)=sdata; 
%     sta=sta+1;
end

for k=1:length(YI)    
data=z_eachsta(:,k);
sdata=jininterp1(slatss,data,XI);    
szdata(:,k)=sdata; 
end

sz = 50;
[sx,sy]=meshgrid(XI,YI);
Canth_grid=reshape(szdata',[],1);
sslat=reshape(sx,[],1);
ssdep=reshape(sy,[],1);
basevalue = 6000;

figure
% scatter(sslat,ssdep,sz,ssz1,'filled')
% scatter(sslat,ssdep,sz,ssz2,'filled')
scatter(sslat,ssdep,sz,Canth_grid,'filled')
colormap(jet)
colorbar
caxis([-0.5 0.1]);
hold on
area(LATITUDE,DEPTH,basevalue,'FaceColor',[0.5 0.5 0.5])
% area(x2,bd2,basevalue,'FaceColor',[0.5 0.5 0.5])
set(gca,'ydir','reverse');
grid on
set(gca,'Fontsize',20,'XLim',[-42 -32],'YLim',[0 6000],'Fontname','Times New Roman')
xlabel('Latitude (¡ãS)','Fontsize',20,'Fontname','Times New Roman');
ylabel('Depth (m)','Fontsize',20,'Fontname','Times New Roman');
box on

sel_depth=[50 150 300 500 700 900 1100 1300 1500 1700 1900 2250 2750 3250 3750 4250 4750];
ind{1}=find(YI>=0 & YI<=100);    % 50
ind{2}=find(YI>=100 & YI<=200);   % 150
ind{3}=find(YI>=200 & YI<=400);   % 300
ind{4}=find(YI>=400 & YI<=600);    % 500
ind{5}=find(YI>=600 & YI<=800);    % 700
ind{6}=find(YI>=800 & YI<=1000);   % 900
ind{7}=find(YI>=1000 & YI<=1200);  % 1100
ind{8}=find(YI>=1200 & YI<=1400);    % 1300
ind{9}=find(YI>=1400 & YI<=1600);   % 1500
ind{10}=find(YI>=1600 & YI<=1800);  % 1700
ind{11}=find(YI>=1800 & YI<=2000);    % 1900
ind{12}=find(YI>=2000 & YI<=2500);   %2250
ind{13}=find(YI>=2500 & YI<=3000);  % 2750
ind{14}=find(YI>=3000 & YI<=3500);    % 3250
ind{15}=find(YI>=3500 & YI<=4000);    % 3750
ind{16}=find(YI>=4000 & YI<=4500);   % 4250
ind{17}=find(YI>=4500 & YI<=5000);  % 4750

for m=1:length(sel_depth)
        pcanth=szdata(:,ind{1,m});
        mspcanth{m}=reshape(pcanth,[],1);  
        my(m)=nanmean(mspcanth{m});
        stdy(m)=nanstd(mspcanth{m});
       
end



figure

top_margin = 0.08; % top margin
btm_margin = 0.20; % bottom margin
left_margin = 0.20;% left margin
right_margin = 0.08;% right margin
fig_margin = 0.001; % margin beween figures(sub)  
row = 2; % rows
col = 1; % cols 

% Calculate figure height and width according to rows and cols 
fig_h = (1- top_margin - btm_margin - (row-1) * fig_margin) / row;
fig_w = (1 - left_margin - right_margin - (col-1) * fig_margin) / col;

% figure position: you can refer to 'help axes' to review the        
% parameter meaning, note that original point is lower-left        
position = [left_margin + (1-1)*(fig_margin+fig_w), ...           
    1- (top_margin + 1 * fig_h + (1-1) * fig_margin), ...           
    fig_w, fig_h];       
axes('position', position)       
% draw colorful pictures...    
hold on;
line([-0.2,-0.2],[0,6000],'Color',[0.8 0.8 0.8]);
line([-0.4,-0.4],[0,6000],'Color',[0.8 0.8 0.8]);
scatter(d13C_eMLR,sdep,sz,'ks');
errorbar(my,sel_depth,stdy,'horizontal','ks-','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',2);
line([0,0],[0,6000],'Color',[0.5 0.5 0.5]);

set(gca,'Fontsize',20,'XLim',[-0.6 0.2],'YLim',[0 2000],'Fontname','Times New Roman')
xticks([ ])
yticks([0 500 1000 1500 2000])

box on
grid on
ax = gca;ax.YDir= 'reverse'; 


position = [left_margin + (1-1)*(fig_margin+fig_w), ...           
    1- (top_margin + 2 * fig_h + (2-1) * fig_margin), ...           
    fig_w, fig_h];       
axes('position', position)
hold on
scatter(d13C_eMLR,sdep,sz,'ks');
errorbar(my,sel_depth,stdy,'horizontal','ks-','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',2);
line([0,0],[0,6000],'Color',[0.5 0.5 0.5]);

set(gca,'Fontsize',20,'XLim',[-0.6 0.2],'YLim',[2000 6000],'Fontname','Times New Roman')
xticks([-0.6 -0.4 -0.2 0 0.2])
yticks([3000 4000 5000 6000])
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
xlabel('Anthropogenic \delta^1^3C change (¡ë)','Fontsize',22,'Fontname','Times New Roman');
% ylabel('Depth (m)','Fontsize',26,'Fontname','Times New Roman');
[ax,h2]=suplabel('Depth (m)','y');
set(h2,'Fontsize',22,'Fontname','Times New Roman','position',[-0.1 0.5])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6])
% print(gcf,'-dtiff','-r300','F3b_eMLR_d13C_change.tif'); 

%%
% % % arrange by station number
% % ssn=STNNBR(in);
% % ssn1=ssn(in1);ssn2=ssn(in2);ssn3=ssn(in3);
% % ssn4=ssn(in4);ssn5=ssn(in5);
% % 
% % stn(1:length(d1))=ssn1;
% % stn(1+length(d1):length(d1)+length(d2))=ssn2;
% % stn(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=ssn3;
% % stn(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=ssn4;
% % stn(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=ssn5;
% % 
% % interp_depth=[15 50 100 150 200 250 300 400 500 600 700 800 900 1000 1100 1200 ...
% %     1400 1600 1800 2000 2200 2400 2600 2800 3000 3200 3400 3600 3800 4000 4200 ...
% %     4400 4600 4800 4900];
% % 
% % for n=1:8
% % instn{n}=find(stn==n);
% % dstn{n}=sdep(instn{n});
% % ystn{n}=d13C_eMLR(instn{n});
% % yyysin(:,n)=interp1(dstn{n},ystn{n},interp_depth);
% % 
% % end
% % 
% % for i=1:35
% % md13C_eMLR(i)=nanmean(yyysin(i,:));
% % stdd13C_eMLR(i)=nanstd(yyysin(i,:));
% % end
% % 
% % sd=interp_depth;
% % sd(isnan(md13C_eMLR))=[];
% % md13C_eMLR(isnan(md13C_eMLR))=[];
% % stdd13C_eMLR(isnan(stdd13C_eMLR))=[];
% % 
% % x_low=md13C_eMLR-stdd13C_eMLR;
% % x_up=flip(md13C_eMLR+stdd13C_eMLR);
% % x_er = [x_low,x_up,x_low(1)];
% % y_er = [sd,flip(sd),sd(1)];
% % 
% % % figure
% % % axes('position', [0.21 0.20 0.72 0.76])  
% % % hold on
% % % p1=plot(d13C_eMLR,sdep,'ks','MarkerSize',8);
% % % p2=errorbar(md13C_eMLR,sd,stdd13C_eMLR,'horizontal','ks-','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',2);
% % % ax = gca;ax.YDir= 'reverse'; 
% % % % legend('2020 s7','2020 s8','Location','southeast',...
% % % %     'NumColumns',1,'Fontsize',16,'Fontname','Times New Roman')
% % % grid on
% % % box on
% % % xticks([-0.6 -0.4 -0.2 0 0.2])
% % % set(gca,'XLim',[-0.6 0.2],'YLim',[0 6000],'Fontsize',22,'Fontname','Times New Roman')
% % % xlabel('Anthropogenic \delta^1^3C change (?)','Fontsize',23,'Fontname','Times New Roman');
% % % ylabel('Depth (m)','Fontsize',23,'Fontname','Times New Roman');
% % % set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6])
% % % % saveas(gcf,['./eMLR_Anthropogenic_d13C_change_A135_all.tif']);
% % % % print(gcf,'-dtiff','-r300','eMLR_d13C_change_A135_mean'); 
% % 
% % figure
% % sz=60;
% % 
% % top_margin = 0.08; % top margin
% % btm_margin = 0.18; % bottom margin
% % left_margin = 0.20;% left margin
% % right_margin = 0.08;% right margin
% % fig_margin = 0.001; % margin beween figures(sub)  
% % row = 2; % rows
% % col = 1; % cols 
% % 
% % % Calculate figure height and width according to rows and cols 
% % fig_h = (1- top_margin - btm_margin - (row-1) * fig_margin) / row;
% % fig_w = (1 - left_margin - right_margin - (col-1) * fig_margin) / col;
% % 
% % % figure position: you can refer to 'help axes' to review the        
% % % parameter meaning, note that original point is lower-left        
% % position = [left_margin + (1-1)*(fig_margin+fig_w), ...           
% %     1- (top_margin + 1 * fig_h + (1-1) * fig_margin), ...           
% %     fig_w, fig_h];       
% % axes('position', position)       
% % % draw colorful pictures...    
% % hold on;
% % scatter(d13C_eMLR,sdep,sz,'ks');
% % % errorbar(md13C_eMLR,sd,stdd13C_eMLR,'horizontal','ks-','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',2);
% % plot(md13C_eMLR,sd,'k-','linewidth',2)
% % h1=fill(x_er,y_er,'k');
% % h1.FaceColor=[0 0 0];
% % h1.FaceAlpha=0.3;
% % % h1.EdgeColor=[1 1 1];
% % h1.EdgeAlpha=0;
% % 
% % 
% % line([0,0],[0,6000],'Color','k');
% % 
% % 
% % set(gca,'Fontsize',22,'XLim',[-0.6 0.2],'YLim',[0 2000],'Fontname','Times New Roman')
% % xlabel([ ])
% % xticks([ ])
% % box on
% % ax = gca;ax.YDir= 'reverse'; 
% % 
% % 
% % position = [left_margin + (1-1)*(fig_margin+fig_w), ...           
% %     1- (top_margin + 2 * fig_h + (2-1) * fig_margin), ...           
% %     fig_w, fig_h];       
% % axes('position', position)
% % hold on
% % scatter(d13C_eMLR,sdep,sz,'ks');
% % % errorbar(md13C_eMLR,sd,stdd13C_eMLR,'horizontal','ks-','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',2);
% % plot(md13C_eMLR,sd,'k-','linewidth',2)
% % h1=fill(x_er,y_er,'k');
% % h1.FaceColor=[0 0 0];
% % h1.FaceAlpha=0.3;
% % % h1.EdgeColor=[1 1 1];
% % h1.EdgeAlpha=0;
% % line([0,0],[0,6000],'Color','k');
% % 
% % set(gca,'Fontsize',22,'XLim',[-0.6 0.2],'YLim',[2000 6000],'Fontname','Times New Roman')
% % xlabel('Latitude','Fontsize',24,'Fontname','Times New Roman');
% % xticks([-0.6 -0.4 -0.2 0 0.2])
% % yticks([3000 4000 5000 6000])
% % box on
% % % hold off
% % ax = gca;ax.YDir= 'reverse'; 
% % % legend('2010 measured','2020 measured','2020 MLR-predicted', ...
% % %    '2010 measured average','2010 measured STD',...
% % %    '2020 measured average','2020 measured STD',...
% % %    '2020 MLR-predicted average','2020 MLR-predicted STD', ...
% % %    'Location','southeast',...
% % %     'NumColumns',1,'Fontsize',16,'Fontname','Times New Roman')
% % xlabel('\delta^1^3C (?)','Fontsize',24,'Fontname','Times New Roman');
% % % ylabel('Depth (m)','Fontsize',26,'Fontname','Times New Roman');
% % [ax,h2]=suplabel('Depth (m)','y');
% % set(h2,'Fontsize',24,'Fontname','Times New Roman','position',[-0.1 0.5])
% % 
% % set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6])
% % % print(gcf,'-dtiff','-r300','eMLR_d13C_change_A135_meanv4.tif'); 
