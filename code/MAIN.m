%%Treball de disseny d'un compressor axial per a l'assignatura de sistemes
%propulsios d'aeronaus.
%Introducció de dades:
TAUc=300000; %J/kg
G=21.38; %kg/s
Patm=1; %kg/cm2
Tatm=288; %K
gamma=1.4;
Rgas=286.8; %J/kgK
Rreac=0.5;
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
        %Rendiment per esglaó:
        ETAesg(i,j)=1-CD(i,j)/CLi(i,j)*(2*flux(j)+1/flux(j));
    end
end




