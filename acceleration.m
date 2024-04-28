% Huvudprogram för simulering av accelerationsförlopp av gokart med c-koppling och växel.
% 
% Programmet anropar funktionsfilerna:
%   1. motormoment.m
%   2. kopplingsmoment.m
%   3. lastmoment.m
%   4. derivataacceleration.m
% 
% Anders Söderberg, KTH - Maskinkonstruktion, 2023-08-21

%% Städar arbetsminnet och command window, samt stänger alla öppna fönster
clear all
clc
close all

%% Sätter värden på systemparametrar

% Definera systemparametar som globala
global eta u r m  JM JK1 JK2 JV1 JV2 JL g f c rho A

% Data för gokart
m   = 131;                  % Total massa [kg]
r   = 0.135;                % Hjulradie [m]
A   = 0.5;                  % Projicerad area för luftmotstånd [m^2]
f   = 0.012;                % Rullmotståndskoefficient [-]
g   = 9.81;                 % Tyngdacceleration [m/s^2]
c   = 0.6;                  % Luftmotståndskoefficient [-]
rho = 1.22;                 % Denstitet hos luft [kg/m^3]

% Masströghetsmoment i systemet:
JM  = 0.02;                 % Tröghetsmoment för motorn [kg*m^2]
JK1 = 0.001;                % Tröghetsmoment för medbringar med block [kg*m^2]
JK2 = 0.002;                % Tröghetsmoment för trumma [kg*m^2]
JV1 = 0.0004;                    % Tröghetsmoment för växel drev [kg*m^2]
JV2 = 0.008;                    % Tröghetsmoment för växel hjul [kg*m^2]
JA  = 0.070;                % Tröghetsmoment för bakaxel [kg*m^2]
JB  = 0.015;                    % Tröghetsmoment för bromsskiva [kg*m^2]
JH  = 0.03;                    % Tröghetsmoment per hjul [kg*m^2]
JL  = JA+JB+4*JH+m*r^2;     % Tröghetsmoment för lasten [kg*m^2]

% Data för motorn
n0   =  1000;               % Tomgångsvarvtal [rpm]

% Data för växeln
u   =    3;                 % Utväxling [-]
eta =    0.95;              % Verkningsgrad i växel [-]


%% Momentkurvor

% Definition av varvtalsvektorer för momentkurvor
n1     = n0:10000;                   % Motorvarvtal [rpm]
n2     = 0:10000;                    % Utgående varvtal från koppling [rpm]
n3     = n2/u;                      % Varvtal på bakaxel [rpm]

% Beräkning av vinkelhastighetersvektorer
omega1 = n1*2*pi/60;                % Vinkelhastighet på motoraxel [1/s]
omega2 = n2*2*pi/60;                % Utgående vinkelhastighet från koppling [1/s]
omega3 = n3*2*pi/60;                % Vinkelhastighet på bakaxel [1/s]
v      = omega3*r;                  % Hastighet hos gokart [m/s]

% Beräkning av momentvektorer
MM     = motormoment(omega1);       % Motormoment [Nm]
MK     = kopplingsmoment(omega1);   % Maximalt överförbart kopplingsmoment [Nm]
ML     = lastmoment(v);             % Lastmoment [Nm]
% 



% % Svar till resultablad
i  = find(MK>0, 1 );
nA =n1(i);
disp(['Kopplingen börjar slira vid n1=' num2str(nA) ' rpm'])
i  = find(MK>=(ML(1)/(eta*u)), 1 );
nB = n1(i);
disp(['Bakaxeln börjar rotera vid n1=' num2str(nB) ' rpm'])
i = find(MK>=MM, 1 );
nC = n1(i);
disp(['Kopplingen slutar slira vid n1=' num2str(nC) ' rpm'])
i = find(n2>=n0);
j = find(MM<=(ML(i)/(eta*u)), 1 );
nD = n1(j);
disp(['Systemet når kontiumerlig drift vid n1=' num2str(nD) ' rpm'])
vD = 2*pi/60*nD/u*r;
disp(['Teoretisk topphastighet hos fordonet är ' num2str(vD*3.6) ' km/h'])

% Plot av momentkurvor
figure(1)
plot(n1,MM,n1,MK,n2,ML/(eta*u),'LineWidth',2); hold on
plot(1265,0,'*',2260,11.3385,'*', 'LineWidth',1)
grid on
legend('Motor','Koppling','Last', 'n1 - begynnande slirning', 'n2 - koppling')
xlabel('n [rpm]')
ylabel('M [Nm]')
ylim([0 12])
xlim([0 ceil(nD/1e3)*1e3])

%% Simulering av accelerationsförlopp
% % Definiera tidsintervall och begynnelsevillkor
tspan  = [0:0.001:90];                      % Start- och sluttid för simulering [s]
IC     = [n0*2*pi/60 0 0 0 0];              % Vid start går motorn på tomgång och lasten står still

% Lös problem med lämplig ode-lösare
opt    = odeset('RelTol',1e-9);             % Sätter toleranser på lösaren
[T,Y] = ode45('derivataacceleration',tspan,IC,opt);    % Anropar ode-lösare

% Dela upp tillståndsmatris på vektorer
Omega1 = Y(:,1);                            % Vinkelhastighet på motoraxel [rad/s]
Omega2 = Y(:,2);                            % Vinklehastighet på kopplingstrumma [rad/s]
Omega3 = Y(:,3);                            % Vinkelhastighet på bakaxel [rad/s]
V      = Y(:,4);                            % Hastighet på gokart [m/s]
S      = Y(:,5);                            % Färdsträcka för gokart [m]

%Beräkna varvtalsvektorer
N1 = Omega1*60/(2*pi);                      % Varvtal på motoraxel [rpm]
N2 = Omega2*60/(2*pi);                      % Varvtal på kopplingstrumma [rpm]
N3 = Omega3*60/(2*pi);                      % Varvtal på bakaxel [rpm]

figure(2)
plot(T,[N1 N2 N3],'LineWidth',2)
grid on
legend('n_1','n_2','n_3')
xlabel('t [s]')
ylabel('n [rpm]')

% Graf som visar gokartenshastighet sfa tid
figure(3)
plot(T,V*3.6,'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('v [km/h]')

% Graf som visar gokartens hastighet sfa körsträcka
figure(4)
plot(S,V*3.6,'LineWidth',2)
grid on
xlabel('s [m]')
ylabel('v [km/h]')

% Svar till resultatblad
disp(' ')
i = find(N1>=nA, 1 );
tA = T(i);
disp(['Kopplingen börjar slira vid n1=' num2str(ceil(N1(i))) ' rpm och t=' num2str(T(i)) ' s'])
i = find(N2>0, 1);
disp(['Lasten börjar rotera vid n1=' num2str(ceil(N1(i))) ' rpm och t=' num2str(T(i)) ' s'])
i = find(N2<=N1, 1, 'last' );
tC = T(i);
disp(['kopplingen slutar slira vid n1=' num2str(ceil(N1(i))) ' rpm och t=' num2str(T(i)) ' s'])
i = find(N1>=0.95*nD, 1);
disp(['Motorn når 95% av kontinuerligtdriftvarvatal (' num2str(ceil(N1(i))) ' rpm) vid t=' num2str(T(i)) ' s'])
disp(['Kopplingen slirar under totalt ' num2str(tC-tA) ' s'])
disp(['Fordonet når 95% hav teoretisk topphastighet efter ' ceil(num2str(S(i))) ' m och ' num2str(T(i)) ' s'])
i = find(V>(50*1e3/3600),1);
disp(['Fordonet når 50 km/h efter ' ceil(num2str(ceil(S(i)))) ' m och ' num2str(T(i)) ' s'])



%% EOF