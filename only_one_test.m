clearvars
close all
clc

%update2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%T=1, R=1, theta_variation%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tL = 100*10^-9;
tR = 1*tL;
nL = 10;
nR = 1.2*nL;
OA = 0;
OB = 100*pi;



c = 3*10^8;

omegal =4*pi*c/(tL*4*nL);
% omegal = c/(nL/2+nR/2) * pi/(tL/2+tR/2);


% omegal = c/(2*pi*1500*10^-9)
% wave = c/(omegal/(2*pi));
% fre = omegal/(2*pi);
% lambda = c/fre*10^9;
% tri = 1;
% omegal =1/tri*4*pi*c/(tL*4*nL);



% 
% lambda = 1500*10^-9;
% omegal = 2*pi*c/(lambda);
kL = nL*omegal/c;
kR = nR*omegal/c;


[locate_bloch, omega] = Band_function(omegal,nL,nR,tL,tR,OA,0);
[locate_bloch_22, omega] = Band_function(omegal,nL,nR,tL,tR,OA,OB);

figure()
plot(locate_bloch(2,:)/pi, omega/omegal, '.')
hold on
plot(locate_bloch_22(2,:)/pi, omega/omegal, '.')
% ylim([0,2])
% xlabel('Normalized Bloch wave vector')
% ylabel('Reduced frequency \omega/\omega_0')
set(gca,'FontSize',20)
% xlim([0.01,1])
% ylim([0,2])
% Create a legend with 3 entries
[h,icons] = legend('\theta_R = 0','\theta_R = 300\pi');
% Find the 'line' objects
icons = findobj(icons,'Type','line');
% Find lines that use a marker
icons = findobj(icons,'Marker','none','-xor');
% Resize the marker in the legend
set(icons,'MarkerSize',20);





band_gap = find(locate_bloch_22(2,:)==0);

for i = 1:length(band_gap)
%     band_0_list(1) = band_gap(2);

    if band_gap(1) ~= 1
        if band_gap(i+1) - band_gap(i) == 1
    %         display(i)
           band_0_list(i) = band_gap(i); 
    %         band_0_list(i) = band_gap(i+1);
        else
           band_0_list(i) = band_gap(i);
           break
        end
    else
        band_0_list(1) = band_gap(2);
        if band_gap(i+2) - band_gap(i+1) == 1
            band_0_list(i+1) = band_gap(i+2);     
        else
%             band_0_list(i+1) = band_gap(i);
            break 
        end
    end
end


ome_list = omega/omegal;
result = ome_list(band_0_list(end)) - ome_list(band_0_list(1));


% 
% figure()
% plot(locate_bloch(2,:)/pi, omega/omegal, '.')
% hold on
% plot(locate_bloch_22(2,:)/pi, omega/omegal, '.')
% ylim([0,2])
% xlabel('Normalized Bloch wave vector')
% ylabel('Reduced frequency \omega/\omega_0')
% set(gca,'FontSize',20)
% % xlim([0.01,1])
% ylim([0,2])
% % Create a legend with 3 entries
% [h,icons] = legend('\theta_R = 0','\theta_R = 300\pi');
% % Find the 'line' objects
% icons = findobj(icons,'Type','line');
% % Find lines that use a marker
% icons = findobj(icons,'Marker','none','-xor');
% % Resize the marker in the legend
% set(icons,'MarkerSize',20);

% figure()
% plot(omega/omegal, locate_bloch(2,:)/pi, '.')