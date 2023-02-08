clearvars
close all
clc



x_disc = 0.1;
X_ratio = 0.1:x_disc:2;

T_total = 100*10^-9;
% tL = 100*10^-9;
% tR = 1*tL;
nL = 10;
nR = 1*nL;
OA = 0;
OB = 100*pi;


R_ratio = 0.1:x_disc:2;

c = 3*10^8;
omegal =4*pi*c/(T_total*4*nL);
omega = (0.1:0.0001:3)*omegal;





for i = 1:length(X_ratio)
    display(i)
    tL = T_total/(X_ratio(i)+1);
    tR = T_total*X_ratio(i)/(X_ratio(i)+1);
    locate_bloch2_X = Band_function(omega,omegal,nL,nR,tL,tR,OA,OB);
    [y_width_X(i),y_center_X(i)] = band_width(locate_bloch2_X,omega,omegal);
end


figure()
for i = 1:length(R_ratio)
    display(i)
    tL = T_total/2;
    tR = T_total/2;
    locate_bloch2_R = Band_function(omega,omegal,nL,R_ratio(i)*nR,tL,tR,OA,OB);
    [y_width_R(i),y_center_R(i)] = band_width(locate_bloch2_R,omega,omegal);
    check{i} = locate_bloch2_R;
    plot(omega/omegal, locate_bloch2_R(2,:)/pi,'.')
    hold on
end

%%
figure()
subplot(2,2,1)
scatter(X_ratio,y_width_X,'MarkerFaceColor',[0 0.4470 0.7410])
xlabel('T_{ratio}')
ylabel('Band gap width')
set(gca,'FontSize',20)
subplot(2,2,2)
scatter(X_ratio,y_center_X,'*','MarkerEdgeColor',[0 0.4470 0.7410])
xlabel('T_{ratio}')
ylabel('Band gap center')
ylim([0.8,1.2])
set(gca,'FontSize',20)
subplot(2,2,3)
scatter(R_ratio,y_width_R,'MarkerFaceColor',[0.9290 0.6940 0.1250])
xlabel('R_{ratio}')
ylabel('Band gap width')
set(gca,'FontSize',20)
subplot(2,2,4)
scatter(R_ratio,y_center_R,'*','MarkerEdgeColor',[0.9290 0.6940 0.1250])
hold on
xlabel('R_{ratio}')
ylabel('Band gap center')
set(gca,'FontSize',20)

% 
% figure()
% subplot(2,1,1)
% scatter(R_ratio,y_width_R,'MarkerFaceColor',[0 0.4470 0.7410])
% ylabel('Band width')
% set(gca,'FontSize',20)
% subplot(2,1,2)
% scatter(R_ratio,y_center_R,'MarkerFaceColor',[0.9290 0.6940 0.1250])
% xlabel('R_{ratio}')
% ylabel('Band gap center')
% set(gca,'FontSize',20)


% close all
% figure(1)
% % plot(1,result_X_1,'Color',"#0072BD")
% scatter(1,result_X_1,'filled','MarkerEdgeColor',[0 0.4470 0.7410])
% hold on
% for i = 2:length(X_ratio)
%     [y_width, y_center] = band_width(nL,nR,tL,X_ratio(i)*tL,OA,OB);
%     figure(1)
%     scatter(X_ratio(i),y_width,'MarkerEdgeColor',[0 0.4470 0.7410])
% end
% set(gca,'filled')