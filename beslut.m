%%plot för samband mellan lamellmassa och fjäderkraft för koppling 
clear all
clc
close all
alfa = 35; %vinkel på medbringare [deg]
rtp = 0.0352; %tyngdpunktsradie lamell [m] 
R = 0.060; %trummans innerradie [m]
n1 = 1265; %varvtalet så kopplingen börjar slira [rpm]
my = 0.28; %friktionskoefficient lamell-trumma []
n2 = n1:3000; %undersökt varvtalsområde [rpm]
Nk=2260;%varvtalet då kopplingen ska koppla
N_b=2; %antal lameller [st]

omega2_rad = n2 *2*pi/60; %konverterar till rad/s
MM = motormoment(omega2_rad); %[Nm]
massa = (MM/(2*my*R)*(1-my*tand(alfa)))./(rtp*((omega2_rad).^2-(omega2_rad(1)).^2))*10^3;
Ffj = (MM/(4*my*R)*(1-my*tand(alfa))) ./ ((omega2_rad).^2/(omega2_rad(1)).^2-1);
%Plottar graferna
yyaxis right
plot(n2, Ffj);
axis([1250 2600 10 800])
yyaxis left
plot(n2,massa)
axis([1250 2600 0 1600])
yyaxis left
title(['Lamellmassa och Fjäderkraft för ' num2str(n1) 'rpm'])
xlabel('Kopplingsvarvtal [rpm]')
ylabel('massa [g]')
yyaxis right
ylabel('Kraft [N]')
grid on
%tar reda på vilken massa som krävs för att koppla vid Nk 
massa = massa(Nk-n1)
Ffj = Ffj(Nk-n1)
k = 12*1000; %Fjaderkonstant [N/m]
Delta = 0.001; %Delta [m]
L0=30.5;
F0 = 14.07; 
CC=(((61.7699-2*k*Delta-F0)/k)*1000)+L0 %försträckt längd