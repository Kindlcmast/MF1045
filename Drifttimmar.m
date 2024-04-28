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
JV2 = 0.013;                    % Tröghetsmoment för växel hjul [kg*m^2]
JA  = 0.070;                % Tröghetsmoment för bakaxel [kg*m^2]
JB  = 0.015;                    % Tröghetsmoment för bromsskiva [kg*m^2]
JH  = 0.03;                    % Tröghetsmoment per hjul [kg*m^2]
JL  = JA+JB+4*JH+m*r^2;     % Tröghetsmoment för lasten [kg*m^2]

% Data för motorn
n0   =  1000;               % Tomgångsvarvtal [rpm]

% Data för växeln
u   =    3;                 % Utväxling [-]
eta =    0.95;              % Verkningsgrad i växel [-]

% Beräkning av vinkelhastighetersvektorer
n1     = n0:10000;                  % Motorvarvtal [rpm]
n2     = 0:10000;                   % Utgående varvtal från koppling [rpm]
omega1 = n1*2*pi/60;                % Vinkelhastighet på motoraxel [1/s]
omega2 = n2*2*pi/60;                % Utgående vinkelhastighet från koppling [1/s]

% Beräkning av momentvektorer
MM     = motormoment(omega1);       % Motormoment [Nm]
MK     = kopplingsmoment(omega1);   % Maximalt överförbart kopplingsmoment [Nm]
ii = find(MK>0,1);
i = find(MK>=MM, 1 ); %kollar upp när kopplingen binder
Mdriv=[MK(1:i-1), MM(i:end)]; %lägger ihop kopplinges momentkurva med motorns momentkurva

tspan  = [0:0.001:90];                      % Start- och sluttid för simulering [s]
IC     = [n0*2*pi/60 0 0 0 0];              % Vid start går motorn på tomgång och lasten står still

% Lös problem med lämplig ode-lösare
opt    = odeset('RelTol',1e-9);             % Sätter toleranser på lösaren
[T,Y] = ode45('derivataacceleration',tspan,IC,opt);    % Anropar ode-lösare

% Dela upp tillståndsmatris på vektorer
Omega1 = Y(:,1);                            % Vinkelhastighet på motoraxel [rad/s]
Omega2 = Y(:,2);                            % Vinklehastighet på kopplingstrumma [rad/s]

%Beräkna varvtalsvektorer
N1 = Omega1.*60/(2*pi);                      % Varvtal på motoraxel [rpm]
N2 = Omega2.*60/(2*pi);                      % Varvtal på kopplingstrumma [rpm]

Mdriv=spline(omega1,Mdriv,Omega1);           % Momentet från koppling till växel


%data för Lager, kedja och lagerinbyggnad

%lager
d1=11.75;                           % avstånd mellan innersta lager och drev [mm]
d2=11.75;                           % avstånd mellan drev och yttersta lager [mm]
s=30000;                            % körsträcka [km]
p=3;                                % Kullager
text=["61806-2RS1","W 61706-2RZ"];  % namnen för de två undersökta lagrena  
C=[4.49,0.65];                      % dynamiskt bärighetstal[kN] SKF 61806-2RS1
C0=[2.9,0.53];                      % statiskt bärighetstal [kN] SKF 61806-2RS1
f0=[14,8.9];                        % f0 Skf data SKF 61806-2RS1
f0FaC0faktorlista=[0.17,0.69,2.08,3.46,5.19];      % vektor med värden från SKF
elista= [0.23,0.3,0.4,0.45,0.5];                   % vektor med värden från SKF
Y1lista=[2.8,2.1,1.6,1.4,1.26];                    % vektor med värden från SKF
Y2lista=[3.7,2.8,2.15,1.85,1.7];                   % vektor med värden från SKF

ns1=2.3;                                        % säkerhetsfaktor från kriterie A B C
ns2=1;                                          % säkerhetsfaktor från kriterie D E
ns=ns1*ns2;                                     % total säkerhetsfaktor

%kedja
k=0.1;                              % faktor k, för spännign av kedja []
delning=12.7;                       % kedjans delning [mm]
mkpm=0.68; %Massa per meter för kedjan
z1=19;                              % antal tänder på drevet [st]
z2=u*z1;                            % antal tänder på hjulet [st]               
alfa=0;                             % vinkel på kedja [deg]
CC=354.75; %centrumavstånd mm
CC=361;
fy=1.5;
fz=1;
fi=1;
fa=1.10;
D0=230.54;                          % delningsdiameter  hjul

fF=[0.5,0.28];                      % tabellvärde för respektive kedja
d0=[77.16,66.97];                   % drevets olika delningsdiametrar för  respektive kedja [mm] 
Pt=[3100 2900 2600 2400 2300 2000 1800 1700 1600 1300 1200 1000]; % lista med tillåtet lagertryck för kedjan vid olika hastigheter
V=[0.1 0.4 0.8 1.5 2 3 4 5 6 8 10 12];                            % hastigheter för respektive lagertryck

% bomförband
Inre_diameter= 46; %mm
Yttre_diameter= 50; %mm
Bomantal= 8; %st
Navlangd= 7.2;%mm

%Kilförband bakaxel
hb=7;%höjd kil mm
Lb=15;%aktiv längd kil mm 
db=30;%diameter på axel mm 
bb=8;%bredd kil mm 
ptill=200; %tillåtet yttryck kilspårsförband [Mpa]

%skruvförband växelhjul
n=6;%antalet skruvar
ds=70*1e-3; %skruvcirkelns diameter m
my=0.4;%statiskt friktionstal mellan växelhjul och nav

%Kilförband koppling
hk=5;%höjd kil mm
Lk=29.5;%aktiv längd kil mm 
dk=20;%diameter på axel mm 
bk=5;%bredd kil mm 
ptill=200; %tillåtet yttryck kilspårsförband [Mpa]

% Beräkning av lagerkrafter från friläggning
vm=((delning*(10^-3)*z1*Omega2)./(2*pi)); %medelhastighet vid respektive kopplingsvarvtal
Fr2= Omega2.*Mdriv./(((1-k)*vm));   % belastning Fr2 vid varje kopplingsvarvtal [N]
Fr1=k.*Fr2;                         % F1 [N]
Fr=Fr1+Fr2;                         % total radiell last från drevet [N]
Fa= Fr*sind(alfa);                  % axiell last över varvtalsområdet [N]
Rr1=((Fr.*d1)/(d1+d2))-Fr;          % radiella reaktionskraften för det inre lagret
Rr2=-(Fr.*d1)/(d1+d2);              % radiella reaktionskraften för det yttre lagret
if Fa>0                             % om kraften riktas utåt kommer axiella kraften gå genom det yttre lagret
    Ra2=Fa;                         % axiella reaktionskraften för det inre lagret
    Ra1=0;                          % axiella reaktionskraften för det yttre lagret
elseif Fa<0                         % om kraften riktas innåt kommer axiella kraften gå genom det inre lagret
    Ra1=Fa;                         % axiella reaktionskraften för det yttre lagret
    Ra2=0;                          % axiella reaktionskraften för det inre lagret
else                                % om den axiella kraften är noll kommer ingen axiell reaktinskraft att uppstå
    Ra1=zeros(size(Fa));            % axiella reaktionskraften för det inre lagret
    Ra2=Ra1;                        % axiella reaktionskraften för det yttre lagret
end

%plottar Lagerkrafterna över motorvarvtalet
figure(1)
plot(N1,abs(Rr1),N1,abs(Rr2), N1, Ra2,N1, Ra1, 'LineWidth',2); hold on; grid on;
title('Lagerkraftdiagram')
xlabel('Motorvarvtal [rpm]')
ylabel('Lagerkraft [N]')
legend('Radiell lagerkraft 1', 'Radiell lagerkraft 2','Axiel lagerkraft lager 2','Axiel lagerkraft lager 1')


%Lagerlivslängd för SKF 61806-2RS1 kontra W 61706-2RZ
%beräknar enligt parade enradiga spårkulager i o-anordning SKF sida 257 Tabell 10

%Dimentionerar lager efter största belastning
if (0==nnz(Ra1)) FA=abs(Ra2); else FA=abs(Ra1);  end
if max(abs(Rr2))<max(abs(Rr1)) FR= abs(Rr1); else FR=abs(Rr2); end

%Beräknar dynamisk lagerlast
l=1;
while l<=2
f0FaC0faktor= f0(l).*FA.*(10^(-3))./C0(l);         % Beräknar faktorn för att kunna interpolera fram rätt koeficienter
e=spline(f0FaC0faktorlista, elista, f0FaC0faktor); % tar fram korresponderande e
if not(0==nnz(FA))
if FA./FR<=e
    Y1=spline(elista,Y1lista,e);
    P(1:90001,l)=((FR)+(Y1.*FA))*(10^(-3));                 %Ekvivalent dynamisk lagerbelastning [kN]
else
    Y2=spline(elista,Y2lista,e); 
    P(1:90001,l)=((0.75.*FR)+(Y2.*FA))*(10^(-3));           %Ekvivalent dynamisk lagerbelastning [kN]
end
else
    P(1:90001,l)=FR.*(10^(-3));                             %Ekvivalent dynamisk lagerbelastning [kN]
end
Pmax(l,:)=max(P(:,end));                       % maximala dynamiska lagerlasten under drift [N]
j(l,:)=find(Pmax(end)==P(:,end),1);            % vilken plats  i listan som ger högsta dynamiska lagerlasten

%% 

%Drifttimmar vid drift vid maximal lagerkraft i 3000mil
Nm=s/(r*u*pi*2000);                % antal rotationer av motoraxeln i miljoner varv
disp(['antal miljoner varv om drift skulle ske vid maximal dynamisk lagerlast i 3000mil' num2str(Nm(end))])

%% 

%Lagerlivslängdsberäkning
L10(:,l)=(C(l)./(Pmax(l)))^p;       % nominell livslängd  vid tillförlitligheten 90%
delta=Nm-L10(end);                  % differens mellan lagerlivslängd och antal rotationer under kravstäld sträcka
disp(text(l))

%avgör hurvida lagret är tillräckligt
if delta<0
    disp("kommer inte att klara 3000 mil vid maximal lagerlast")
    disp(['lagret saknar: ',num2str(abs(delta)),' miljoner varv i livslängd'])
else
    disp("kommer att klara 3000 mil vid maximal lagerlast")
    disp(['marginal: ',num2str(abs(delta)),' miljoner varv'])
end
disp(['L10 =' num2str(L10(end))])
l=l+1;
end
%% 

%plottar Lagerkrafterna över motorvarvtalet
figure(2)
%plot(N1,P(1:end,1), N1,P(1:end,2),'LineWidth',2); hold on; grid on;
%plot(N1(j(1)),Pmax(1),'*',N1(j(2)),Pmax(2),'*')
plot(N1,P(1:end,2),'LineWidth',2); hold on; grid on;
plot(N1(j(2)),Pmax(2),'*')
title('dynamisk lagerlast över varvtalsområdet')
xlabel('motorvarvtal [rpm]')
ylabel('lagerlast [kN]')
%legend('dynamisk Lagerlast 61806-2RS1','dynamisk Lagerlast W 61706-2RZ' ,'maximal dynamisk lagerlast 61806-2RS1', 'maximal dynamisk lagerlast W 61706-2RZ')
legend('dynamisk Lagerlast W 61706-2RZ' , 'maximal dynamisk lagerlast W 61706-2RZ')

%___Kedja____
jj=1;
while jj<=2
nP=V.*(19098./d0(jj));                          % varvtal för olika tillåtna lagertryck
PL=(Mdriv./(d0(jj)./2*10^-3))./fF(jj);          % uppkommet lagertryck
Ptn=spline(nP,Pt,N2);                           % splinar tillåtet lagertryck för varvtalsområdet
delta(1:90001,jj)=Ptn-PL;                       % tar fram deltan mellan tillåtet lagertryck och uppkommet lagertryck för varje varvtal
delta0(jj,:)=round(spline((delta(n1(ii):end,jj)),transpose(N2(n1(ii):end)),0)); %tar fram det varvtal där deltan blir noll
jj=jj+1;
end
P=max(Mdriv.*Omega2)/1000; %maximal effekt som går ut i kedjan, [kw]
k=fy*fz*fi*fa+0.5; 
Pd=P*k; %diagrameffekt
disp(['diagrameffekt kedja: ', num2str(Pd),' [kw]'])
L=2*CC+1.57*(d0(1)+D0)+((d0(1)+D0)^2)/(4*CC); %Längd kedja [mm]
X=round(L/delning); %antal kedjelänka
disp(['antal kedjelänkar för kedja med delning ' num2str(delning) 'mm : ' num2str(X)])
kedjemassa = mkpm*L/1000; %Beräknar kedjas massa
disp(['vikt kedja: ' num2str(kedjemassa) ' kg'])
%% 

%plottar de två kurvorna för utvärdering
figure(3) 
plot(N2,delta(1:end,1),'LineWidth',1.5)
hold on 
plot(delta0(1),0,'o','LineWidth',1)
plot(N2,delta(1:end,2),'LineWidth',1.5)
plot(delta0(2),0,'o','LineWidth',1)
legend('08B-1',num2str(delta0(1)),'06B-1',num2str(delta0(2)))
xlabel('motorvarvtal [rpm]')
ylabel('delta tillåtet - uppkommet lagertryck [N/cm^2]')
title('Kedjors lagertryck under drift')
ylim([-2000, 2000])
xlim([n1(ii),6000])
grid on
%% 


%bomförband
M= max(Mdriv).*1000; %Nmm
yttrycktill=30; %MPa
h=(Yttre_diameter-Inre_diameter)*0.5; %mm
Resulterande_yttryck=M./(0.75*Bomantal*h*Navlangd*Inre_diameter*0.5); %MPa
% Mamimalt_overforbart_moment=Maximaltyttryck*Inre_diameter*0.5*0.75*Bomantal*h*Navlangd/1000;%Nm
nbom=yttrycktill/max(Resulterande_yttryck);
disp(['resulterande säkeretsfaktor för bomförband koppling-drev ' num2str(nbom)])
%% 


%Kilspår bakaxel
M=Mdriv.*3*eta;%moment för bakaxeln
p=(4*1000*M)/(hb*Lb*db); %yttryck MPa
Nphjul=ptill/max(p);
disp(['resulterande säkeretsfaktor för kilföbarband växelhjulsnav-bakaxel ' num2str(Nphjul)])
%% 

%skruvförband
Ni=max((2*M)/(my*n*ds)); %erfoderlig kraft i skruv [N]
disp(['erforderlig axiell kraft i varje skruv ' num2str(Ni) ' [N]'])
%% 

%Kilspår motoraxel
M=Mdriv;%moment för motoraxel
p=(4*1000*M)/(hk*Lk*dk); %yttryck MPa 
Npkoppling=ptill/(max(p));
disp(['resulterande säkeretsfaktor för kilförband koppling-motoraxel ' num2str(Npkoppling)])

%% 


































