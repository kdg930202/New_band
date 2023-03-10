function locate_bloch2 = Band_function(omega,omegal,nL,nR,tL,tR,OA,OB)

% nL = 5;
% nR = 5;


% tL = 100*10^-9; % Fix
% tR = 100*10^-9;
c = 3*10^8;




% omegal = c/(nL/2+nR/2) * pi/(tL/2+tR/2);
% omegal = 2*pi*c/(tL*4*nL);


% lambda_omegal = 2*pi*c/omegal*10^9;
% lambda_omegal_nm = lambda_omegal*10^9;
%Ultraviolet 750~30,000 THz 400~10 nm
%Visible     400~750 THZ    750~400nm


% omegal_tera = 10^-12*omegal;

dis = 0.0001;
% omega = (0:dis:2)*omegal;
bloch = (-1:dis:1)*pi;

uA = 1;
uB = 1;

nA = nL;
nB = nR;
ZA = uA/nA;
ZB = uB/nB;
dA = tL;
dB = tR;


% OA = 0;
% OB = 1000;


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

    [number_zeros,locate_zero] = crossing(LHS-RHS(i),bloch); % Zero finding algorithm
    locate_bloch(:,i) = locate_zero;
end


for j = 1:length(omega)

    if locate_bloch(:,j) == [0,0]
        locate_bloch2(:,j) = [0,0];
    else
        locate_bloch2(:,j) = [ bloch(locate_bloch(1,j)), bloch(locate_bloch(2,j)) ];
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

% figure()
% scatter(locate_bloch2(2,:)/pi, omega/omegal, '.')
% xlim([0.01,1])
% ylim([0,1])




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
    for l=2:1:length(z)
        if z(l) > z(l-1)+round((samplerate/10000))+1
            a = 1;
            if l == 2
                zci(b,a) = z(l-1);
            end
            zci(end+1,a) = z(l); 
            b = b+1;
        else
            zci(b,a) = z(l); 
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
        for k=1:1:length(zci2)
            if zci2(k) == 0 && k ~=1
                break
            end
            zci3(b,k) = array(zci2(k));
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




end