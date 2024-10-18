close all;clear;clc;

load('DCanth_from_d13C_RC_density.mat') 
mCantly=md13C_eMLR;
stdCantly=stdd13C_eMLR;

load('DCanth_from_DIC_density.mat') 

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
scatter(d13C_eMLR,sden,sz,'ko');
errorbar(md13C_eMLR,sd,stdd13C_eMLR,'horizontal','ko','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',3);
scatter(Cant_eMLR,den_gamma,sz,'ro');
errorbar(mCantly,sd,stdCantly,'horizontal','ro','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','r','linewidth',2);
line([0,0],[24,28.5],'Color',[0.5 0.5 0.5]);
line([10,10],[24,28.5],'Color',[0.8 0.8 0.8]);
line([20,20],[24,28.5],'Color',[0.8 0.8 0.8]);
line([15,15],[24,28.5],'Color',[0.8 0.8 0.8]);
line([5,5],[24,28.5],'Color',[0.8 0.8 0.8]);
line([-5,-5],[24,28.5],'Color',[0.8 0.8 0.8]);

set(gca,'Fontsize',20,'XLim',[-10 25],'YLim',[24.8 26.8],'Fontname','Times New Roman')
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
scatter(d13C_eMLR,sden,sz,'ko');
errorbar(md13C_eMLR,sd,stdd13C_eMLR,'horizontal','ko','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',3);

scatter(Cant_eMLR,den_gamma,sz,'ro');
errorbar(mCantly,sd,stdCantly,'horizontal','ro','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','r','linewidth',2);

line([0,0],[24,28.5],'Color',[0.5 0.5 0.5]);

set(gca,'Fontsize',20,'XLim',[-10 25],'YLim',[26.8 28.3],'Fontname','Times New Roman')
% xticks([-0.6 -0.4 -0.2 0 0.2])
yticks([27.23 27.5 28 28.27])
xticks([-10 -5 0 5 10 15 20 25])
box on
grid on
% hold off
ax = gca;ax.YDir= 'reverse'; 
xlabel('Anthropogenic CO_2 change (\mumol kg^-^1)','Fontsize',22,'Fontname','Times New Roman');
% ylabel('Depth (m)','Fontsize',26,'Fontname','Times New Roman');
[ax,h2]=suplabel('Neutral density (kg m^{-3})','y');
set(h2,'Fontsize',22,'Fontname','Times New Roman','position',[-0.1 0.5])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6])
print(gcf,'-dtiff','-r300','Canth_change_density_compare_F3c.tif'); 

