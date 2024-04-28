clear all
clc
close all
b=23; %yttre radie
a=37/2; %inre radie
ts=190; %tillåten skjuvspänning
ns=2.3; %säkerhetsfaktor
A=pi*(-(a^2)+(b^2));
Fr1=34.467; %F1 [N]
Fr2=344.6730;% F2 [N]
Fr=Fr1+Fr2;
M=11.3;
sigma=sqrt(((M/((pi/(2*b))*((10^(-3))*((b^4)-(a^4)))))^2)+((Fr/A)^2))*ns;
sigmatill=ts/ns;
n=sigmatill/sigma;



