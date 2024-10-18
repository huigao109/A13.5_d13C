close all;clear;clc;
% % load 2010 data
load('A13.5_2010.mat')

% inn=find(depthm>=2000 & c13~=-9999 & phosphate~=-9999);  
inn=find(depthm>=2000 & c13~=-9999 & phosphate~=-9999 & latitudedegrees_north>=-42 & latitudedegrees_north<=-32);  

density=gamma1(inn);
d13C=c13(inn);
PO4=phosphate(inn);

% [rho,pval] = corr(PO4,d13C);

% ig1=find(density<26.8);
% ig2=find(density>=26.8 & density<27.23);
% ig3=find(density>=27.23 & density<27.5);
% ig4=find(density>=27.5 & density<28);
% ig5=find(density>=28);
% 
% [p1,s1] = polyfit(PO4(ig1),d13C(ig1),1);
% [p2,s2] = polyfit(PO4(ig2),d13C(ig2),1);
% [p3,s3] = polyfit(PO4(ig3),d13C(ig3),1);
% [p4,s4] = polyfit(PO4(ig4),d13C(ig4),1);
% [p5,s5] = polyfit(PO4(ig5),d13C(ig5),1);

[p,s] = polyfit(PO4,d13C,1);

% load 2020 data
load('A13.5_2020v5.mat')

in=find(cal_depth>=2000 & DELC13~=-999 & PHSPHT~=-999);   

density2=cal_gamma(in);
d13C2=DELC13(in);
PO42=PHSPHT(in);

% [rho,pval] = corr(PO4,d13C);

% id1=find(density2<26.8);
% id2=find(density2>=26.8 & density2<27.23);
% id3=find(density2>=27.23 & density2<27.5);
% id4=find(density2>=27.5 & density2<28);
% id5=find(density2>=28);


% % d13C_pre1=p1(1)*(PO42(id1))+p1(2);
% % d13C_pre2=p2(1)*(PO42(id2))+p2(2);
% % d13C_pre3=p3(1)*(PO42(id3))+p3(2);
% % d13C_pre4=p4(1)*(PO42(id4))+p4(2);
% % d13C_pre5=p5(1)*(PO42(id5))+p5(2);
% 
% d13C_pre1=polyval(p1,PO42(id1));
% d13C_pre2=polyval(p2,PO42(id2));
% d13C_pre3=polyval(p3,PO42(id3));
% d13C_pre4=polyval(p4,PO42(id4));
% d13C_pre5=polyval(p5,PO42(id5));

d13C_pre=polyval(p,PO42);
[rho,pval] = corr(d13C_pre,d13C2)

% [rho1,pval1] = corr(d13C_pre1,d13C2(id1))
% [rho2,pval2] = corr(d13C_pre2,d13C2(id2))
% [rho3,pval3] = corr(d13C_pre3,d13C2(id3))
% [rho4,pval4] = corr(d13C_pre4,d13C2(id4))
% [rho5,pval5] = corr(d13C_pre5,d13C2(id5))

figure
hold on
% plot(d13C_pre1,d13C2(id1),'ko')
% plot(d13C_pre2,d13C2(id2),'ro')
% plot(d13C_pre3,d13C2(id3),'bo')
% plot(d13C_pre4,d13C2(id4),'go')
% plot(d13C_pre5,d13C2(id5),'mo')
plot(d13C_pre,d13C2,'bs')
line([0 2],[0 2],'Color','k','LineStyle','-')
xlim([0,1.6]);ylim([0,1.6])
xticks(0:0.2:1.6);yticks(0:0.2:1.6);
box on
set(gca,'Fontsize',20,'Fontname','Times New Roman')
xlabel('MLR_{predicted} \delta^1^3C','Fontsize',20,'Fontname','Times New Roman');
ylabel('\delta^1^3C','Fontsize',20,'Fontname','Times New Roman');
% legend('SW: \gamma<26.8','SAMW: 26.8<=\gamma<27.23','AAIW: 27.23<=\gamma<27.5', ...
%     'UCDW: 27.5<=\gamma<28','LCDW+AABW: \gamma>=28','', ...
%    'Location','southeast','Fontsize',12,'Fontname','Times New Roman');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 6])
print(gcf,'-dtiff','-r300','d13C_pred13C_2000_32_42'); 


