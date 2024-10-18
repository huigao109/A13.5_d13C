close all;clear;clc;

load('DCanth_from_d13C_RC_depth_v2.mat') 

load('DCanth_from_DIC_depth_v2.mat') 


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
line([-5,-5],[0,6000],'Color',[0.8 0.8 0.8]);
line([5,5],[0,6000],'Color',[0.8 0.8 0.8]);
line([10,10],[0,6000],'Color',[0.8 0.8 0.8]);
line([15,15],[0,6000],'Color',[0.8 0.8 0.8]);
line([20,20],[0,6000],'Color',[0.8 0.8 0.8]);


scatter(d13C_eMLR,sdep,sz,'ko');
errorbar(my,sel_depth,stdy,'horizontal','ko-','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',3);
scatter(Cant_eMLR,dep,sz,'ro');
errorbar(my2,sel_depth,stdy2,'horizontal','ro-','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','r','linewidth',2);
line([0,0],[0,6000],'Color',[0.5 0.5 0.5]);

set(gca,'Fontsize',20,'XLim',[-10 25],'YLim',[0 2000],'Fontname','Times New Roman')
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
scatter(d13C_eMLR,sdep,sz,'ko');
errorbar(my,sel_depth,stdy,'horizontal','ko-','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k','linewidth',3);
scatter(Cant_eMLR,dep,sz,'ro');
errorbar(my2,sel_depth,stdy2,'horizontal','ro-','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','r','linewidth',2);
line([0,0],[0,6000],'Color',[0.5 0.5 0.5]);

set(gca,'Fontsize',20,'XLim',[-10 25],'YLim',[2000 6000],'Fontname','Times New Roman')
% xticks([-0.6 -0.4 -0.2 0 0.2])
xticks([-10 -5 0 5 10 15 20 25])
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
xlabel('Anthropogenic CO_2 change (\mumol kg^-^1)','Fontsize',22,'Fontname','Times New Roman');
% ylabel('Depth (m)','Fontsize',26,'Fontname','Times New Roman');
[ax,h2]=suplabel('Depth (m)','y');
set(h2,'Fontsize',22,'Fontname','Times New Roman','position',[-0.1 0.5])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6])
print(gcf,'-dtiff','-r300','Canth_change_depth_compare_F3a.tif'); 




