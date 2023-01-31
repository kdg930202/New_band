clearvars
% close all
clc

%update2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%T=1, R=1, theta_variation%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tL = 100*10^-9;
nL = 10;
c = 3*10^8;
% omegal =4*pi*c/(tL*4*nL);
% fre = omegal/(2*pi);
% % % lambda = c/fre*10^9;
lambda = 2000*10^-9;
fre = c/lambda;
omegal = 2*pi*fre;


[locate_bloch, omega] = Band_function(omegal,nL,nL,tL,tL,0,0);
[locate_bloch_2, omega] = Band_function(omegal,nL,nL,tL,tL,0,100*pi);
[locate_bloch_3, omgea] = Band_function(omegal,nL,nL,tL,tL,0,200*pi);



omegal_tera = 10^-12*omegal; %Angular frequency 
f_tera = omegal_tera/(2*pi); %375 THZ : Frequency

lambda_omegal = 2*pi*c/omegal*10^9; %800nm

figure()
subplot(1,3,1)
plot(locate_bloch(2,:)/pi, omega/omegal, '.')
hold on
plot(locate_bloch_2(2,:)/pi, omega/omegal, '.')
plot(locate_bloch_3(2,:)/pi, omega/omegal, '.')
% legend({'$\theta_R = 0$','$\theta_R = 100\pi$','$\theta_R = 200\pi$'},'Interpreter','latex')
title('T = 1, R = 1')
xlabel('Normalized Bloch wave vector')
ylabel('Reduced frequency \omega/\omega_0')
set(gca,'FontSize',20)
xlim([0.01,1])
ylim([0,2])
% Create a legend with 3 entries
[h,icons] = legend('\theta_R = 0','\theta_R = 100\pi','\theta_R = 200\pi');
% Find the 'line' objects
icons = findobj(icons,'Type','line');
% Find lines that use a marker
icons = findobj(icons,'Marker','none','-xor');
% Resize the marker in the legend
set(icons,'MarkerSize',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%T=1, R=1, theta_variation%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%T=1, R=2, theta_variation%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[locate_bloch_4, omega] = Band_function(omegal, nL, 1.5*nL, tL, tL,0,pi);
[locate_bloch_5, omega] = Band_function(omegal, nL, 1.5*nL, tL, tL,0,500*pi);
[locate_bloch_6, omega] = Band_function(omegal, nL, 1.5*nL, tL, tL,0,1000*pi);


subplot(1,3,2)
plot(locate_bloch_4(2,:)/pi, omega/omegal, '.')
hold on
plot(locate_bloch_5(2,:)/pi, omega/omegal, '.')
plot(locate_bloch_6(2,:)/pi, omega/omegal, '.')
% legend({'$\theta_R = 0$','$\theta_R = 100\pi$','$\theta_R = 200\pi$'},'Interpreter','latex')
title('T = 1, R = 1.5')
xlabel('Normalized Bloch wave vector')
ylabel('Reduced frequency \omega/\omega_0')
set(gca,'FontSize',20)
xlim([0.01,1])
ylim([0,2])
% Create a legend with 3 entries
[h,icons] = legend('\theta_R = 0','\theta_R = 100\pi','\theta_R = 200\pi');
% Find the 'line' objects
icons = findobj(icons,'Type','line');
% Find lines that use a marker
icons = findobj(icons,'Marker','none','-xor');
% Resize the marker in the legend
set(icons,'MarkerSize',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%T=1, R=2, theta_variation%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





[locate_bloch_4, omega] = Band_function(omegal, nL, 1.5*nL,tL,2*tL,0,0);
[locate_bloch_5, omega] = Band_function(omegal, nL, 1.5*nL,tL,2*tL,0,100*pi);
[locate_bloch_6, omega] = Band_function(omegal, nL, 1.5*nL,tL,2*tL,0,200*pi);
subplot(1,3,3)
plot(locate_bloch_4(2,:)/pi, omega/omegal, '.')
hold on
plot(locate_bloch_5(2,:)/pi, omega/omegal, '.')
plot(locate_bloch_6(2,:)/pi, omega/omegal, '.')
% legend({'$\theta_R = 0$','$\theta_R = 100\pi$','$\theta_R = 200\pi$'},'Interpreter','latex')
title('T = 2, R = 1.5')
xlabel('Normalized Bloch wave vector')
ylabel('Reduced frequency \omega/\omega_0')
set(gca,'FontSize',20)
xlim([0.01,1])
ylim([0,2])
% Create a legend with 3 entries
[h,icons] = legend('\theta_R = 0','\theta_R = 100\pi','\theta_R = 200\pi');
% Find the 'line' objects
icons = findobj(icons,'Type','line');
% Find lines that use a marker
icons = findobj(icons,'Marker','none','-xor');
% Resize the marker in the legend
set(icons,'MarkerSize',20);

