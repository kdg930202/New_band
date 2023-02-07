clearvars
close all
clc

tL = 100*10^-9;
% tR = 1*tL;
nL = 10;
nR = 1*nL;
OA = 0;
OB = 100*pi;


X_ratio = 0.1:0.1:2;
R_ratio = 0.1:0.1:2;

c = 3*10^8;
omegal =4*pi*c/(tL*4*nL);
omega = (0:0.001:3)*omegal;





for i = 1:length(X_ratio)
    locate_bloch2_X = Band_function(omega,omegal,nL,nR,tL,X_ratio(i)*tL,OA,OB);
    [y_width_X(i),y_center_X(i),a] = band_width(locate_bloch2_X,omega,omegal);
end



for i = 1:length(R_ratio)
    locate_bloch2_R = Band_function(omega,omegal,nL,R_ratio(i)*nR,tL,tL,OA,OB);
    [y_width_R(i),y_center_R(i),b] = band_width(locate_bloch2_R,omega,omegal);
end

%%
figure()
subplot(2,2,1)
scatter(X_ratio,y_width_X,'MarkerFaceColor',[0 0.4470 0.7410])
xlabel('T_{ratio}')
ylabel('Band width')
set(gca,'FontSize',20)
subplot(2,2,2)
scatter(X_ratio,y_center_X,'*','MarkerFaceColor',[0 0.4470 0.7410])
xlabel('T_{ratio}')
ylabel('Band gap center')
set(gca,'FontSize',20)
subplot(2,2,3)
scatter(R_ratio,y_width_R,'MarkerFaceColor',[0.9290 0.6940 0.1250])
xlabel('R_{ratio}')
ylabel('Band width')
set(gca,'FontSize',20)
subplot(2,2,4)
scatter(R_ratio,y_center_R,'*','MarkerFaceColor',[0.9290 0.6940 0.1250])
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