clc
clearvars
close all

%

p = 20;
theta_vac = 0;


c = 3*10^8;


type = 1;

if type == 1
    theta_m = pi*(1:1:2*p);
    OA = theta_m(1);
    OB = theta_m(2);
elseif type == 2
    theta_m = pi*(1:2:(4*p-1));
    OA = theta_m(1);
    OB = theta_m(2);
elseif type == 3
    theta_m = repmat([pi,0],1,p);
    OA = theta_m(1);
    OB = theta_m(2);
end



nL = 10;
nR = 0.5*nL;

T_total = 100*10^-9;
T_ratio = 1;
tL = 1/(T_ratio+1)*T_total;
tR = T_ratio/(T_ratio+1)*T_total;
omegal = 4*pi*c/(T_total*4*nL);
% omegal = c/(nL) * pi/(tL);
% omegal = c/(nL/2+nR/2) * pi/(tL/2+tR/2);
fre = omegal/(2*pi);
lambda = c/fre*10^9;


% omegal = 2*pi*c/(tL*4*nL);
omega = (0:0.001:10)*omegal;


n_sub = sqrt(13);



ref = zeros(1,length(theta_m));
for i = 1:length(ref)
    if mod(i,2) == 1
        ref(i) = nL;
    elseif mod(i,2) == 0
        ref(i) = nR;
    end
end





%%%%%%%%%%%%%%%%%%%%%%   

% thick = l/(length(theta_m)-1);
for i = 1:(length(theta_m)-1)
    M_interface(:,:,i) = interface(theta_m(i),theta_m(i+1),ref(i),ref(i+1));
end

M_ti_va = interface(theta_vac,theta_m(1),1,ref(1));
% M_sub_ti = interface(theta_m(end),theta_sub,ref(end),n_sub);

for i=1:length(omega)
    arguL = exp(1i*nL*omega/c*tL);
    arguR = exp(1i*nR*omega/c*tR);
    
    ML(:,:,i) = diag([arguL(i)',arguL(i)',arguL(i),arguL(i)]);
    MR(:,:,i) = diag([arguR(i)',arguR(i)',arguR(i),arguR(i)]);
      
    
    M_total(:,:,i) = MR(:,:,i) * M_interface(:,:,1) * ML(:,:,i);
    for j= 1:p-1
        M_total(:,:,i) = MR(:,:,i) * M_interface(:,:,2*j+1) * ML(:,:,i) * M_interface(:,:,2*j) * M_total(:,:,i);
    end
        
    if type == 3
        M_total(:,:,i) = ML(:,:,i)* interface(theta_m(end),theta_m(end)+pi,ref(end),ref(1)) * M_total(:,:,i);
        M_sub_ti = interface(theta_m(end)+pi,0,ref(1),n_sub);
        M_total(:,:,i) = M_sub_ti * M_total(:,:,i) * M_ti_va;
    else
%         M_total(:,:,i) = ML(:,:,i)* interface(theta_m(end),theta_m(end)+pi,ref(end),ref(1)) * M_total(:,:,i);
        M_sub_ti = interface(theta_m(end),theta_m(end)+pi,ref(end),n_sub);
        M_total(:,:,i) = M_sub_ti * M_total(:,:,i) * M_ti_va;
    end
%     
    
    % ??? ????????? ?????? ????????? ??? ???????????????????
    

% 
%     M_total(:,:,i) = MR(:,:,i) * M_interface(:,:,1) * ML(:,:,i);


%     if type == 1 || type ==3
%         if p>=2
%             for j= 1:p-1
%                 M_total(:,:,i) = MR(:,:,i) * M_interface(:,:,2*j+1) * ML(:,:,i) * M_interface(:,:,2*j) * M_total(:,:,i);
%             end
%         end
%         M_total(:,:,i) = ML(:,:,i)* interface(theta_m(end),theta_m(end)+pi,ref(end),ref(1)) * M_total(:,:,i);
%         M_total(:,:,i) = interface(theta_m(end)+pi,theta_m(end)+2*pi,ref(1),n_sub) * M_total(:,:,i) * M_ti_va;
% 
%     else
%         if p>=2
%             for j= 1:p-1
%                 M_total(:,:,i) = MR(:,:,i) * M_interface(:,:,2*j+1) * ML(:,:,i) * M_interface(:,:,2*j) * M_total(:,:,i);
%             end
%         end
%         M_total(:,:,i) = M_sub_ti * M_total(:,:,i) * M_ti_va;
%     end
    
    
    
    
    
    
        


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    m = M_total(:,:,i);
    M11 = m(1:2, 1:2);
    M12 = m(1:2, 3:4);
    M21 = m(3:4, 1:2);
    M22 = m(3:4, 3:4);
    
    S11 = M11 - M12*inv(M22)*M21;
    S12 = M12*inv(M22);   
    
    tt1(i) = S11(1,1);
    tt2(i) = S11(1,2);
    tt3(i) = S11(2,1);
    tt4(i) = S11(2,2);
    
    rr1(i) = S12(1,1);
    rr2(i) = S12(1,2);
    rr3(i) = S12(2,1);
    rr4(i) = S12(2,2);

    TT1(i) = tt1(i)*conj(tt1(i));
    TT2(i) = tt2(i)*conj(tt2(i));
    TT3(i) = tt3(i)*conj(tt3(i));
    TT4(i) = tt4(i)*conj(tt4(i));
    
    RR1(i) = rr1(i)*conj(rr1(i));
    RR2(i) = rr2(i)*conj(rr2(i));
    RR3(i) = rr3(i)*conj(rr3(i));
    RR4(i) = rr4(i)*conj(rr4(i));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



TT1 = n_sub*TT1;
TT2 = n_sub*TT2;
TT3 = n_sub*TT3;
TT4 = n_sub*TT4;


figure()
plot(omega/omegal,TT1+TT2+RR1+RR2,'LineWidth',2,'Color','#0072BD')
ylim([0.9,1.1])


figure()
subplot(5,1,1)
plot(omega/omegal,TT1,'LineWidth',2,'Color','#0072BD')
title('Transmission(TM->TM)')
set(gca,'FontSize',20)
% title(strcat(num2str(lambda),'nm'))
subplot(5,1,2)
plot(omega/omegal,TT2,'LineWidth',2,'Color','#0072BD')
title('Transmission(TM->TE)')
set(gca,'FontSize',20)
% title(strcat(num2str(lambda),'nm'))
subplot(5,1,3)
plot(omega/omegal,RR1,'LineWidth',2,'Color',"#D95319")
title('Reflection(TM->TM)')
set(gca,'FontSize',20)
% title(strcat(num2str(lambda),'nm'))
subplot(5,1,4)
plot(omega/omegal,RR2,'LineWidth',2,'Color',"#D95319")
title('Reflection(TM->TE)')
set(gca,'FontSize',20)
% figure()
% subplot(1,5,1)
% plot(omega/omegal,TT1,'LineWidth',2)
% ylabel('Transmission(TM->TM)')
% set(gca,'FontSize',20)
% % xlim([0,omega(end)/omegal])
% 
% 
% subplot(1,5,2)
% plot(omega/omegal,RR1,'LineWidth',2,'Color',"#D95319")
% ylabel('Reflection(TM->TM)')
% set(gca,'FontSize',20)
% % xlim([0,omega(end)/omegal])
% 
% subplot(1,5,3)
% plot(omega/omegal,TT2,'LineWidth',2)
% ylabel('Transmission(TM->TE)')
% set(gca,'FontSize',20)
% % xlim([0,omega(end)/omegal])
% 
% subplot(1,5,4)
% plot(omega/omegal,RR2,'LineWidth',2,'Color',"#D95319")
% ylabel('Reflection(TM->TE)')
% set(gca,'FontSize',20)
% % xlim([0,omega(end)/omegal])
% 
% subplot(1,5,5)
% plot(omega/omegal,RR1+RR2,'LineWidth',2,'Color',"#D95319")
% ylabel('Reflection(TM->TM + TE)')
% set(gca,'FontSize',20)
% % xlim([0,omega(end)/omegal])

% subplot(1,3,2)
% plot(omega/omegal,TT2,'LineWidth',2)
% hold on
% plot(omega/omegal,RR2,'LineWidth',2)
% set(gca,'FontSize',30)
% subplot(1,3,3)
% plot(omega/omegal,TT1+RR1+TT2+RR2,'LineWidth',2)
% ylim([0.7,1.3])
% set(gca,'FontSize',30)



clearvars -except nL nR tL tR OA OB omegal omega theta_m lambda uL uR




locate_bloch = Band_function(omega, omegal,nL,nR,tL,tR,OA,0);
locate_bloch_22 = Band_function(omega, omegal,nL,nR,tL,tR,OA,OB);




subplot(5,1,5)
plot(omega/omegal, locate_bloch(2,:)/pi,'k.')
hold on
% plot(omega/omegal, locate_bloch_22(2,:)/pi,'.')
% plot(omega/omegal, locate_bloch_22(2,:)/pi,'Color',[0 0.4470 0.7410])
% legend({'n_R = 0',strcat('n_R = ', num2str(OB/pi),'\pi')})
set(gca,'FontSize',20)
title(strcat(num2str(lambda),'nm'))
% strcat('n_R = ', num2str(OB/pi),'\pi')



function MRL = interface(thetaL,thetaR,nL,nR) % interface matrix for left to right
    
    c = 1; %????????? ???????????? divided ????????? ?????? ??? ?????? ????????????
    
    alpha = 1/137;
    u0 = 1;
    

    uL = 1;
    uR = 1;
    aL = nL/(u0*uL*c);
    bL = alpha*thetaL/(u0*c*pi);
    
    aR = nR/(u0*uR*c);
    bR = alpha*thetaR/(u0*c*pi);
    
    MRL = 1/(2*aR)*...
      [aR+aL, -bR+bL, aR-aL, -bR+bL;
       bR-bL, aR+aL, bR-bL, aR-aL;
       aR-aL, bR-bL, aR+aL, bR-bL;
       -bR+bL, aR-aL, -bR+bL, aR+aL
      ];

end