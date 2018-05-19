%%Treball de disseny d'un compressor axial per a l'assignatura de sistemes
%propulsios d'aeronaus.
%Introducció de dades:
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
%Inicialització dels valor comuns de solidesa i coeficient de flux
sigma=[0.4 0.6 0.8 1 1.2];
flux=[0.4 0.5 0.6 0.7 0.8];
N=size(sigma,2);
for i=1:N %Per a les N sigma
    for j=1:N %Per als N coef de flux
        %Cálcul betas
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
        %Rendiment per esglaó
        ETAesg(i,j)=1-CD(i,j)/CLi(i,j)*(2*flux(j)+1/flux(j));
        %Velocitat axial. Aquesta velocitat haurá d'estar entre 150 i 180
        %m/s
        Vz(i,j)=sqrt((M^2*gamma*Rgas*Tat)/(1/(cos(betA(i,j)))^2+(M^2*gamma*Rgas)/(2*Cp*(cos(betB(i,j)))^2)));
        %Velocitat tangencial. El valor máxim d'aquesta velocitat es 320
        %m/s
        U(i,j)=Vz(i,j)/flux(j);
        %Treball per esglaó
        TAUesg(i,j)=U(i,j)*Vz(i,j)*(tan(betA(i,j))-tan(betB(i,j)));
        rire(i,j)=(U(i,j)^2-sigmAl*g/(2*lambda*rho))/(U(i,j)^2+sigmAl*g/(2*lambda*rho)); %OJO unidades
        %A partir de aqui no funciona
        Wa(i,j)=Vz(i,j)/cos(betA(i,j));
        Va(i,j)=sqrt((Vz(i,j)*sin(betA(i,j))/cos(betA(i,j))-U(i,j))^2+Vz(i,j)^2);
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
%Print surfaces: 
figure;
surf(sigma,flux,betA);
figure;
surf(sigma,flux,betB);
figure;
surf(sigma,flux,CL);
figure;
surf(sigma,flux,CD);
figure;
surf(sigma,flux,ETAesg);
figure;
surf(sigma,flux,Vz);
figure;
surf(sigma,flux,U);
figure;
surf(sigma,flux,TAUesg);
figure;
surf(sigma,flux,rire);
figure;
surf(sigma,flux,re);
figure;
surf(sigma,flux,ri);
figure;
surf(sigma,flux,rm);
figure;
surf(sigma,flux,h);
figure;
surf(sigma,flux,RPM);


