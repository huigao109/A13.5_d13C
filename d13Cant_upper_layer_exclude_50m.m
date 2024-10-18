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
y=c13(in);
dep=depthm(in);

density=gamma1(in);
% in1=find(density<26.8 & dep>=50);
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

%%
clear in x y in1 pres depth
% load 2020 data

load('A13.5_2020v5.mat')

in=find(CTDTMP~=-999 & CTDSAL~=-999 & AOU~=-999 & NITRAT~=-999 & SILCAT~=-999 & d13C_lab_Najid~=-999); 
x(:,1)=CTDTMP(in);
x(:,2)=CTDSAL(in);
x(:,3)=AOU4_S(in);
x(:,4)=NITRAT(in);
x(:,5)=SILCAT(in);
DIC=TCARBN(in);
pres=CTDPRS(in);
depth=cal_depth(in);
% d13C=d13C_lab(in);
% y=d13C_lab(in);
% y=d13C_lab(in)+0.05;
% y=d13C_lab_Najid(in);
y=d13C_lab_Najid(in)+0.07;   % offset
% y=d13C_lab_Najid(in)+0.05;   % offset

ssn=STNNBR(in);
lats=LATITUDE(in);

density=cal_gamma(in);
% in1=find(density<26.8 & depth>=50);
in1=find(density<26.8);
in2=find(density>=26.8 & density<27.23);
in3=find(density>=27.23 & density<27.5);
in4=find(density>=27.5 & density<28);
% in5=find(density>=28 & density<28.27);
% in6=find(density>=28.27);
in5=find(density>=28);

x2020{1}=x(in1,1:5);y2020{1}=y(in1);
x2020{2}=x(in2,1:5);y2020{2}=y(in2);
x2020{3}=x(in3,1:5);y2020{3}=y(in3);
x2020{4}=x(in4,1:5);y2020{4}=y(in4);
x2020{5}=x(in5,1:5);y2020{5}=y(in5);
% x2020{6}=x(in6,1:5);y2020{6}=y(in6);


%%  anthropogenic d13C change
for slab=1:5   
    x1=x2010{1,slab};
    y1=y2010{1,slab};
    x2=x2020{1,slab};
    y2=y2020{1,slab};

   [b1,stats1] = robustfit(x1,y1);  
   coe1(slab,:)=b1'; 
   [b2,stats2] = robustfit(x2,y2);  
   coe2(slab,:)=b2'; 
   
   resi1=stats1.resid;
   resi2=stats2.resid;
   rmse1(slab)=stats1.ols_s;
   rmse2(slab)=stats2.ols_s;
   mresi1(slab)=nanmean(resi1);
   stdresi1(slab)=nanstd(resi1);
   mresi2(slab)=nanmean(resi2);
   stdresi2(slab)=nanstd(resi2);
   
    pd13C1{slab}=coe1(slab,1)+sum(coe1(slab,2:6).*x1,2);
    pd13C2{slab}=coe2(slab,1)+sum(coe2(slab,2:6).*x2,2);

    [rho1(:,slab),pval1(:,slab)] = corr(y1,pd13C1{slab});
    [rho2(:,slab),pval2(:,slab)] = corr(y2,pd13C2{slab});
   
   coediff(slab,:)=b2'-b1';
   preddelC{slab}=coe1(slab,1)+sum(coe1(slab,2:6).*x2,2);
   chd13C{slab}=y2-preddelC{slab};   % MLR method   % measured-predicted
   delC{slab}=coediff(slab,1)+sum(coediff(slab,2:6).*x2,2);  % eMLR method

      clear  x1 x2 y1 y2
end

% eMLR anthropogenic d13C change
d1=delC{1};d2=delC{2};d3=delC{3};d4=delC{4};d5=delC{5};
d13C_eMLR(1:length(d1))=d1;
d13C_eMLR(1+length(d1):length(d1)+length(d2))=d2;
d13C_eMLR(1+length(d1)+length(d2):length(d1)+length(d2)+length(d3))=d3;
d13C_eMLR(1+length(d1)+length(d2)+length(d3):length(d1)+length(d2)+length(d3)+length(d4))=d4;
d13C_eMLR(1+length(d1)+length(d2)+length(d3)+length(d4):length(d1)+length(d2)+length(d3)+length(d4)+length(d5))=d5;

% each station depth
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


% % ind1=find(sdep>=0 & sdep<=100);    % 50
% % ind2=find(sdep>=100 & sdep<=200);   % 150
% % ind3=find(sdep>=200 & sdep<=400);   % 300
% % ind4=find(sdep>=400 & sdep<=600);    % 500
% % ind5=find(sdep>=600 & sdep<=800);    % 700
% % ind6=find(sdep>=800 & sdep<=1000);   % 900
% % ind7=find(sdep>=1000 & sdep<=1200);  % 1100
% % ind8=find(sdep>=1200 & sdep<=1400);    % 1300
% % ind9=find(sdep>=1400 & sdep<=1600);   % 1500
% % ind10=find(sdep>=1600 & sdep<=1800);  % 1700
% % ind11=find(sdep>=1800 & sdep<=2000);    % 1900
% % ind12=find(sdep>=2000 & sdep<=2500);   %2250
% % ind13=find(sdep>=2500 & sdep<=3000);  % 2750
% % ind14=find(sdep>=3000 & sdep<=3500);    % 3250
% % ind15=find(sdep>=3500 & sdep<=4000);    % 3750
% % ind16=find(sdep>=4000 & sdep<=4500);   % 4250
% % ind17=find(sdep>=4500 & sdep<=5000);  % 4750
% % 
% % sd=[50 150 300 500 700 900 1100 1300 1500 1700 1900 2250 2750 3250 3750 4250 4750];
% % 
% % yy{1}=d13C_eMLR(ind1);yy{2}=d13C_eMLR(ind2);yy{3}=d13C_eMLR(ind3);
% % yy{4}=d13C_eMLR(ind4);yy{5}=d13C_eMLR(ind5);yy{6}=d13C_eMLR(ind6);
% % yy{7}=d13C_eMLR(ind7);yy{8}=d13C_eMLR(ind8);yy{9}=d13C_eMLR(ind9);
% % yy{10}=d13C_eMLR(ind10);yy{11}=d13C_eMLR(ind11);yy{12}=d13C_eMLR(ind12);
% % yy{13}=d13C_eMLR(ind13);yy{14}=d13C_eMLR(ind14);yy{15}=d13C_eMLR(ind15);
% % yy{16}=d13C_eMLR(ind16);yy{17}=d13C_eMLR(ind17);
% % 
% % for i=1:17
% % md13C_eMLR(i)=nanmean(yy{i});
% % stdd13C_eMLR(i)=nanstd(yy{i});
% % end

figure
axes('position', [0.21 0.23 0.72 0.72])  
hold on
% p1=plot(d13C_eMLR,sdep,'ks','MarkerSize',8);
% p2=errorbar(md13C_eMLR,sd,stdd13C_eMLR,'horizontal','ks-','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',2);
scatter(d13C_eMLR,sdep,sz,'ks');
errorbar(my,sel_depth,stdy,'horizontal','ks-','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',2);

ax = gca;ax.YDir= 'reverse'; 
% legend('2020 s7','2020 s8','Location','southeast',...
%     'NumColumns',1,'Fontsize',16,'Fontname','Times New Roman')
grid on
box on
xticks([-0.6 -0.4 -0.2 0 0.2])
set(gca,'XLim',[-0.6 0.2],'YLim',[0 1000],'Fontsize',20,'Fontname','Times New Roman')
xlabel('Anthropogenic \delta^1^3C change (¡ë)','Fontsize',20,'Fontname','Times New Roman');
ylabel('Depth (m)','Fontsize',23,'Fontname','Times New Roman');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 4])
% saveas(gcf,['./eMLR_Anthropogenic_d13C_change_A135_all.tif']);
% print(gcf,'-dtiff','-r300','eMLR_d13C_change_exclude_50m'); 
print(gcf,'-dtiff','-r300','eMLR_d13C_change_extrapolate'); 
