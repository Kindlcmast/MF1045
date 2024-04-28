% function MK = kopplingsmoment(omega1)
% % Funktionsfil för beräkning av maximalt överförbart kopplingsmoment som
% % funktion av vinkelhastigheten hos medbringaren.

function MK = kopplingsmoment(omega1)
N_b = 2; %Antal block
R = 0.060; %Trum radie [m]
my = 0.28; %Friktionstal mellan block och trumma
mb =0.199;%Block massa [kg]
rtp = 0.0352; %Tyngdpunktsradie [m]
k = 12*1000; %Fjaderkonstant [N/m]
Delta = 0.001; %Delta [m]
alfa=35; %vinkel på medbringare i grader
L0=30.5;
CCfjader=32.4750;
Forspanning_Delta=(CCfjader-L0)*(10^-3); %[m]
F0 = 14.07;
Forspanningskraft=F0+k*Forspanning_Delta; %Forspanningskraft [N]
Ffj = 2*k*Delta + Forspanningskraft; %Fjaderkraft [N]
MK = N_b * R * my * ((mb * rtp *omega1.^2 - 2*Ffj)./(1-(my*tand(alfa)))); %Kopplingsmoment
if MK<0 MK=0; end
end


