clc; clear; close all;
%% Treball de disseny d'un compressor axial per a l'assignatura de sistemes
%propulsios d'aeronaus.
%Introduccio de dades:
g=9.81; %m/s2
TAUc=300000; %J/kg
G=21.38; %kg/s
Pat=1; %kg/cm2
Pat=Pat*g*10000;
Tat=288; %K
gamma=1.4;
Rgas=286.8; %J/kgK
Rreac=0.5;
M=0.8;
Cp=1004; %J/kgK
lambda=0.7;
sigmAl=(39/4)*10^6; %kg/m2
rho=2800; %kg/m2
%Inicialitzacio dels valor comuns de solidesa i coeficient de flux
sigma=[0.4 0.6 0.8 1 1.2];
flux=[0.4 0.5 0.6 0.7 0.8];
N=size(sigma,2);
for i=1:N %Per a les N sigma
    for j=1:N %Per als N coef de flux
        %Calcul betas
        betA(i,j)=atan(1/2*(1.55/(1+1.55*sigma(i))+1/flux(j)));
        betB(i,j)=atan(1/flux(j)-tan(betA(i,j)));
        betM(i,j)=atan(1/(2*flux(j)));
        %Calcul CL, CLi i CD
        a=tan(betM(i,j))*0.018;
        c=tan(betM(i,j))*(0.021*0.02/2.5*sigma(i))-2*sigma(i)*(tan(betA(i,j))-tan(betB(i,j)))*cos(betM(i,j));
        p=[a 1 c];
        root=roots(p);
        CL(i,j)=root(2);
        CD(i,j)=0.021+0.02/2.5*sigma(i)+0.018*CL(i,j)^2;
        CLi(i,j)=2*sigma(i)*(tan(betA(i,j))-tan(betB(i,j)))*cos(betM(i,j));
        %Rendiment per esglao
        ETAesg(i,j)=1-CD(i,j)/CLi(i,j)*(2*flux(j)+1/(2*flux(j)));
        %Velocitat axial. Aquesta velocitat haura d'estar entre 150 i 180
        %m/s
        Vz(i,j)=sqrt((M^2*gamma*Rgas*Tat)/(1/(cos(betA(i,j)))^2+(M^2*gamma*Rgas)/(2*Cp*(cos(betB(i,j)))^2)));
        %Velocitat tangencial. El valor maxim d'aquesta velocitat es 320
        %m/s
        U(i,j)=Vz(i,j)/flux(j);
        %Treball per esglao
        TAUesg(i,j)=U(i,j)*Vz(i,j)*(tan(betA(i,j))-tan(betB(i,j)));
        rire(i,j)=(U(i,j)^2-sigmAl*g/(2*lambda*rho))/(U(i,j)^2+sigmAl*g/(2*lambda*rho)); %OJO unidades
        
        Wa(i,j)=Vz(i,j)/cos(betA(i,j));
        Wb(i,j)=Vz(i,j)/cos(betB(i,j));
        Wm(i,j)=(Wa(i,j)+Wb(i,j))/2;
        
        Va(i,j)=sqrt((Vz(i,j)*sin(betA(i,j))/cos(betA(i,j))-U(i,j))^2+Vz(i,j)^2);
        Vb(i,j)=sqrt((Vz(i,j)*sin(betB(i,j))/cos(betB(i,j))-U(i,j))^2+Vz(i,j)^2); %REVISAR
        
        Ta(i,j)=Tat-Va(i,j)^2/(2*Cp);
        Pa(i,j)=Pat/(1+Wa(i,j)^2/(2*Rgas*Ta(i,j)));
        rhoa(i,j)=Pa(i,j)/(Rgas*Ta(i,j));
        re(i,j)=sqrt(G/(pi*(1-rire(i,j)^2)*Vz(i,j)*rhoa(i,j)));
        ri(i,j)=rire(i,j)*re(i,j);
        h(i,j)=re(i,j)-ri(i,j);
        rm(i,j)=(re(i,j)+ri(i,j))/2;
        RPM(i,j)=60*U(i,j)/(2*pi*rm(i,j));        
    end
end

%% Seleccio de parametres
%print -depsc2 myplot.eps
n = [ 7 8 9 10];
TAUstage = TAUc./n;
for ii = 1 : 4
    TAUstagePlot(:,ii) = TAUstage(ii)*ones(N,1);
end

% Seleccionar flux de tall amb el necessari (mirant el grafic anterior)
FLUX(1) = 0;
FLUX(2) = interpolarPunt( 3.75e4, 3.788e4, 3.701e4, 0.6, 0.5 ); %[Y, Xa, Xb, Yb, Ya]
FLUX(3) = interpolarPunt( 3.333e4, 3.348e4, 3.224e4, 0.7, 0.6 );
FLUX(4) = interpolarPunt( 3e4, 3.051e4, 2.939e4, 0.7, 0.6 );
FLUX(5) = 0;

%Trobar rendiments per cada flux
for ii = 1 : 5
    FLUXplot(:,ii) = FLUX(ii)*ones(N,1);
end
% (mirant el grafic anterior)
ETA(1) = 0;
ETA(2) = interpolarPunt( FLUX(2), 0.5, 0.6, 0.9026, 0.9003 );
ETA(3) = interpolarPunt( FLUX(3), 0.6, 0.7, 0.9022, 0.9041 );
ETA(4) = interpolarPunt( FLUX(4), 0.6, 0.7, 0.9013, 0.9036 );
ETA(5) = 0;

Cd(1) = 0;
Cd(2) = interpolarPunt( FLUX(2), 0.5, 0.6, 0.03545, 0.03396 );
Cd(3) = interpolarPunt( FLUX(3), 0.6, 0.7, 0.04168, 0.0401 );
Cd(4) = interpolarPunt( FLUX(4), 0.6, 0.7, 0.04618, 0.04427 );

CLi_n(1) = 0;
CLi_n(2) = interpolarPunt( FLUX(2), 0.5, 0.6, 0.7404, 0.6815 );
CLi_n(3) = interpolarPunt( FLUX(3), 0.6, 0.7, 0.9009, 0.8505 );
CLi_n(4) = interpolarPunt( FLUX(4), 0.6, 0.7, 0.9892, 0.9339 );

CL_n(1) = 0;
CL_n(2) = interpolarPunt( FLUX(2), 0.5, 0.6, 0.7322, 0.6732 );
CL_n(3) = interpolarPunt( FLUX(3), 0.6, 0.7, 0.8906, 0.8398 );
CL_n(4) = interpolarPunt( FLUX(4), 0.6, 0.7, 0.9769, 0.9211 );

ETA23(1) = 1 - n(2)*Cd(2)/CLi_n(2)*(2*FLUX(2)+1/(2*FLUX(2))); % 8 etapas
ETA23(2) = 1 - n(3)*Cd(3)/CLi_n(3)*(2*FLUX(3)+1/(2*FLUX(3))); % 9 etapas
ETA23(3) = 1 - n(4)*Cd(4)/CLi_n(4)*(2*FLUX(4)+1/(2*FLUX(4))); % 10 etapas


%% Numero d'aleps en primer esglao i longitud total
% ------------- 8 etapas
et = 8;
sigma_e(1) = round((Cd(2)-0.021-0.018*CL_n(2)^2)*2.5/0.02,1);
index_i = find(sigma == sigma_e(1));
index_j = find(flux == round(FLUX(2),1));
h_e(1) = h(index_i,index_j);
r_e(1) = rm(index_i,index_j);
N(1) = alabes(h_e(1),r_e(1),sigma_e(1));

% ------------- 9 etapas
sigma_e(2) = round((Cd(3)-0.021-0.018*CL_n(3)^2)*2.5/0.02,1);
index_i = find(sigma == sigma_e(2));
index_j = find(flux == round(FLUX(3),1));
h_e(2) = h(index_i,index_j);
r_e(2) = rm(index_i,index_j);
N(2) = alabes(h_e(2),r_e(2),sigma_e(2));
% ------------- 10 etapas
sigma_e(3) = round((Cd(4)-0.021-0.018*CL_n(4)^2)*2.5/0.02,1);
index_i = find(sigma == sigma_e(3));
index_j = find(flux == round(FLUX(4),1));
h_e(3) = h(index_i,index_j);
r_e(3) = rm(index_i,index_j);
N(3) = alabes(h_e(3),r_e(3),sigma_e(3));

N_blades = fix(N);
