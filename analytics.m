clc
clearvars
close all

%git hub challenge
%please


nL = 5;
nR = 5;


tL = 100*10^-9; % Fix
tR = 100*10^-9;
c = 3*10^8;



% % 
% omegal = c/(nL/2+nR/2) * pi/(tL/2+tR/2);
% omegal = 0.5*c*pi/((nL/2+nR/2)*(tL/2+tR/2));
% omegal = 2*pi*c/(tL*4*nL);
omegal = 2356*10^12;

%%
% lambda_omegal = 2*pi*c/omegal;
% lambda_omegal_nm = lambda_omegal*10^9;
% Ultraviolet 750 ~ 30,000 THz 400 ~ 10 nm
% Visible     400 ~ 750 THZ    750 ~ 400nm


% omegal_tera = 10^-12*omegal;

dis = 0.001;
omega = (-2:dis:2)*omegal;
bloch = (-1:dis:1)*pi;


nA = nL;
nB = nR;
ZA = 1/nA;
ZB = 1/nB;
dA = tL;
dB = tR;


OA = 0;
OB = 0*pi;


alpha = 1/137;

kA = nA*omega/c;
kB = nB*omega/c;


% delta = alpha*(OA-OB)^2/(4*pi^2);

delta = alpha^2*(OA-OB)^2/(pi^2); %With this definition, it corresponds to the numerical result. 
% delta = pi^2*(OA-OB)^2/alpha^2;


%%
% delta = 20;
DELTA = 0.5*(ZA/ZB + ZB/ZA + delta*ZA*ZB);

LHS = cos(bloch);
RHS = cos(kA*dA).*cos(kB*dB) - DELTA.*sin(kA*dA).*sin(kB*dB);



locate_bloch = zeros(2,length(omega));
locate_bloch2 = zeros(2,length(omega));
for i = 1:length(omega)
%     display(i)
    [number_zeros,locate_zero] = crossing(LHS-RHS(i),bloch); % Zero finding algorithm
    locate_bloch(:,i) = locate_zero;
end


for i = 1:length(omega)
%     display(i)
    if locate_bloch(:,i) == [0,0]
        locate_bloch2(:,i) = [0,0];
    else
        locate_bloch2(:,i) = [ bloch(locate_bloch(1,i)), bloch(locate_bloch(2,i)) ];
    end
end

% band_data = locate_bloch2(2,:)/pi;
% load('data_type3.mat')
% figure(1)
% subplot(2,1,1)
% plot(omega/omegal, TT1,'LineWidth',2)
% set(gca,'FontSize',30)
% subplot(2,1,2)
% plot(omega/omegal, band_data,'LineWidth',2)
% set(gca, 'YDir','reverse')
% set(gca, 'XAxisLocation', 'top')
% ylabel('Normalized Bloch wave-vector')
% set(gca,'FontSize',15)
% set(gca,'FontSize',30)

figure()
scatter(locate_bloch2(2,:)/pi, omega/omegal, '.')
xlim([0.01,1])
ylim([0,1])




function [number_zeros,zero_crossings] = crossing(array,samplerate)
%FINDZEROS finds zerocrossings
%Finds the zeros or the nearest values to zero in a function and gives back
%as result the number of zerocrossings and an array containing median of the
%array with the positions of the value that are zero or nearst to zero in
%a zero crossing area, so its the middlest value of the zero crossing area
z = find(diff(sign(array)));   

if isempty(z) == 1
    number_zeros = 0;
    zero_crossings = [0,0];

else
    a = 1;
    b = 1;
    for i=2:1:length(z)
        if z(i) > z(i-1)+round((samplerate/10000))+1
            a = 1;
            if i == 2
                zci(b,a) = z(i-1);
            end
            zci(end+1,a) = z(i); 
            b = b+1;
        else
            zci(b,a) = z(i); 
            a = a+1;
        end
    end
    number_zeros = b; %output1
    zci2 = [];
    zb = [];
    zc = [];
    zero_crossings = [];
    for b = 1:1:number_zeros
        zci2 = zci(b,:);
        for j=1:1:length(zci2)
            if zci2(j) == 0 && j ~=1
                break
            end
            zci3(b,j) = array(zci2(j));
        end
        zb = find(abs(zci3(b,:)) == min(abs(zci3(b,:))));
        zb = zci2(zb);
        if length(zb) <= 1
            zero_crossings = [zero_crossings zb]; %output2
        else
            zero_crossings(end+1) = zb(floor(length(zb)/2)); %outpu2
        end
    end
end

end


function MRL = interface(thetaL,thetaR,nL,nR,pol) % interface matrix for left to right
    
    c = 1;
    alpha = 1/137;
    uL = 1;
    uR = 1;
    u0 = 1;


    aL = nL/(u0*uL*c);
    bL = alpha*thetaL/(u0*c*pi);
    
    aR = nR/(u0*uR*c);
    bR = alpha*thetaR/(u0*c*pi);
    
    if pol == 0
    
        MRL = 1/(2*aR)*...
          [aR+aL, -bR+bL, aR-aL, -bR+bL;
           bR-bL, aR+aL, bR-bL, aR-aL;
           aR-aL, bR-bL, aR+aL, bR-bL;
           -bR+bL, aR-aL, -bR+bL, aR+aL
          ];
    else
        MRL = 1/(2*aR)*...
      [aR+aL, 1i*pol*(-bR+bL), aR-aL, 1i*pol*(-bR+bL);
       -1i*bR/pol + 1i*bL/pol, aR+aL, -1i*bR/pol + 1i*bL/pol, aR-aL;
       aR-aL, 1i*pol*(bR-bL), aR+aL, 1i*pol*(bR-bL);
       1i*bR/pol - 1i*bL/pol, aR-aL, 1i*bR/pol - 1i*bL/pol, aR+aL
      ];
    end
  
  

end