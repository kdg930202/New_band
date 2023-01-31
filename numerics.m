clc
clearvars
% close all


nL = 5;
nR = 5;


tL = 100*10^-9; % Fix
tR = 300*10^-9;
c = 3*10^8;
alpha = 1/137;

omegal = c/(nL/2+nR/2) * pi/(tL/2+tR/2);

dis = 0.001;
omega = (-2:dis:2)*omegal;
bloch = (-1:dis:1)*pi;

OA = 0;
OB = 0;

gA = alpha/(c*pi)*OA - 1i*nL/c;
gB = alpha/(c*pi)*OB - 1i*nR/c;

kA = nL*omega/c;
kB = nR*omega/c;


MAB = 1/(2*1i*imag(gB))*[gA - gB', gA' - gB';
                         -gA + gB, -gA' + gB];
MBA = 1/(2*1i*imag(gA))*[gB - gA', gB' - gA';
                         -gB + gA, -gB' + gA];
                     
for i=1:length(omega)
    arguL = exp(1i*nL*omega/c*tL);
    arguR = exp(1i*nR*omega/c*tR);
    
    MA(:,:,i) = diag([arguL(i)',arguL(i)]);
    MB(:,:,i) = diag([arguR(i)',arguR(i)]);
    M_total(:,:,i) = MB(:,:,i)*MBA*MA(:,:,i)*MAB;
    M_index = M_total(:,:,i);
    DELTA(i) = 0.5*(M_index(1,1) + M_index(2,2));
    
end


LHS = cos(bloch);
RHS = DELTA;

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
