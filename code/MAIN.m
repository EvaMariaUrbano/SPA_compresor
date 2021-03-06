clc; clear; close all;
addpath(genpath('./plots'))
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
%export values
matrix2latex(betA, 'matrix2Latex/betA.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(betB, 'matrix2Latex/betB.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(betM, 'matrix2Latex/betM.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(CL, 'matrix2Latex/CL.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(CD, 'matrix2Latex/CD.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(CLi, 'matrix2Latex/CLi.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(Vz, 'matrix2Latex/Vz.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(U, 'matrix2Latex/U.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(TAUesg, 'matrix2Latex/TAUesg.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(rire, 'matrix2Latex/rire.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(Wa, 'matrix2Latex/Wa.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(Wb, 'matrix2Latex/Wb.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(Wm, 'matrix2Latex/Wm.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(Va, 'matrix2Latex/Va.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(Vb, 'matrix2Latex/Vb.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(Ta, 'matrix2Latex/Ta.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(Pa, 'matrix2Latex/Pa.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(re, 'matrix2Latex/re.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(ri, 'matrix2Latex/ri.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(h, 'matrix2Latex/h.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(rm, 'matrix2Latex/rm.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
matrix2latex(RPM, 'matrix2Latex/RPM.tex', 'alignment', 'c', 'format', '%-6.2f', 'size', 'large')
%% Seleccio de parametres
%print -depsc2 myplot.eps
n = [ 7 8 9 10];
TAUstage = TAUc./n;
for ii = 1 : 4
    TAUstagePlot(:,ii) = TAUstage(ii)*ones(N,1);
end
%Treure gr�fics per l'informe
    %Dibuixar normal
figure();
hold on
plot(flux, TAUesg'); 
plot(flux, TAUstagePlot);
hold off
savefig('./figures/TAUSesg.fig')
close
    %Dibuixar bonic
plt = Plot('./figures/TAUS.fig');
plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = 'Flux, \Psi';
plt.YLabel = 'Treball escal�, \tau';
plt.LineWidth(:) = 2;
plt.LineStyle(6:9) = {':',':',':',':'}; %:,-,--,-.
plt.Legend = {'S/C = 0.4','S/C = 0.6','S/C = 0.8','S/C = 1','S/C = 1.2', ...
    '\tau_{N7}','\tau_{N8}','\tau_{N9}','\tau_{N10}',};
%plt.Legend(6:9) = {'n = 7','n = 8','n = 9','n = 10'};
plt.Title = '\tau en funci� de S/C i \Psi';

print -depsc2 figures/parametres/TAUSesg.eps %guardadr foto per l'informe

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
figure()
hold on
plot(flux, ETAesg');
plot(FLUXplot(:,2:end-1), linspace(0.88,0.91,N) );
hold off
savefig('./figures/FLUXs.fig')
close
    %Dibuixar bonic
plt = Plot('./figures/FLUXs.fig');
plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = 'Flux, \Psi';
plt.YLabel = 'Rendiment, \eta_{esg}';
plt.LineWidth(:) = 2;
plt.LineStyle(6:9) = {':',':',':',':'}; %:,-,--,-.
plt.Legend = {'S/C = 0.4','S/C = 0.6','S/C = 0.8','S/C = 1','S/C = 1.2', ...
   '\Psi_{N8}','\Psi_{N9}','\Psi_{N10}',};
plt.Title = '\eta en funci� de S/C i \Psi';

print -depsc2 figures/parametres/FLUXs.eps %guardadr foto per l'informe
% (mirant el grafic anterior)
ETA(1) = 0;
ETA(2) = interpolarPunt( FLUX(2), 0.5, 0.6, 0.9026, 0.9003 );
ETA(3) = interpolarPunt( FLUX(3), 0.6, 0.7, 0.9022, 0.9041 );
ETA(4) = interpolarPunt( FLUX(4), 0.6, 0.7, 0.9013, 0.9036 );
ETA(5) = 0;

%Rendiment total del compressor
    %interpolar CD
% figure()
% hold on
% plot(flux, CD');
% plot(FLUXplot(:,2:end-1), linspace(0.025,0.055,N) );
% hold off
Cd(1) = 0;
Cd(2) = interpolarPunt( FLUX(2), 0.5, 0.6, 0.03545, 0.03396 );
Cd(3) = interpolarPunt( FLUX(3), 0.6, 0.7, 0.04168, 0.0401 );
Cd(4) = interpolarPunt( FLUX(4), 0.6, 0.7, 0.04618, 0.04427 );

%interpolar CLi
% figure()
% hold on
% plot(flux, CLi');
% plot(FLUXplot(:,2:end-1), linspace(0.4,1,N) );
% hold off

CLi_n(1) = 0;
CLi_n(2) = interpolarPunt( FLUX(2), 0.5, 0.6, 0.7404, 0.6815 );
CLi_n(3) = interpolarPunt( FLUX(3), 0.6, 0.7, 0.9009, 0.8505 );
CLi_n(4) = interpolarPunt( FLUX(4), 0.6, 0.7, 0.9892, 0.9339 );

%interpolar CL
% figure()
% hold on
% plot(flux, CL');
% plot(FLUXplot(:,2:end-1), linspace(0.4,1.1,N) );
% hold off

CL_n(1) = 0;
CL_n(2) = interpolarPunt( FLUX(2), 0.5, 0.6, 0.7322, 0.6732 );
CL_n(3) = interpolarPunt( FLUX(3), 0.6, 0.7, 0.8906, 0.8398 );
CL_n(4) = interpolarPunt( FLUX(4), 0.6, 0.7, 0.9769, 0.9211 );

ETA23(1) = 1 - n(2)*Cd(2)/CLi_n(2)*(2*FLUX(2)+1/(2*FLUX(2))); % 8 etapas
ETA23(2) = 1 - n(3)*Cd(3)/CLi_n(3)*(2*FLUX(3)+1/(2*FLUX(3))); % 9 etapas
ETA23(3) = 1 - n(4)*Cd(4)/CLi_n(4)*(2*FLUX(4)+1/(2*FLUX(4))); % 10 etapas

VZ = interpolarPunt( FLUX(2), flux(2), flux(3), Vz(2,2), Vz(2,3) ); %[X, Xa, Xb, Yb, Ya]
Uval = interpolarPunt( FLUX(2), flux(2), flux(3), U(2,2), U(2,3) ); %[X, Xa, Xb, Yb, Ya]
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

%% Print surfaces:
figure;
surf(sigma,flux,betA);
view(40,22)
savefig('./figures/betA.fig')
close
% plt = Plot('figure', 'true');
    %Dibuixar bonic
plt = Plot('./figures/betA.fig', 'true');
% plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = '1/Solidesa, S/C';
plt.YLabel = 'Flux, \Psi';
plt.ZLabel = '\beta_{a}';
plt.LineWidth(:) = 2;
plt.Title = '\beta_{a} en funci� de S/C i \Psi';

print -depsc2 figures/parametres/betA.eps %guardar foto per l'informe
%---------------
figure;
surf(sigma,flux,betB);
view(-156,45)
savefig('./figures/betB.fig')
close
% plt = Plot('figure', 'true');
    %Dibuixar bonic
plt = Plot('./figures/betB.fig', 'true');
% plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = '1/Solidesa, S/C';
plt.YLabel = 'Flux, \Psi';
plt.ZLabel = '\beta_{b}';
plt.LineWidth(:) = 2;
plt.Title = '\beta_{b} en funci� de S/C i \Psi';

print -depsc2 figures/parametres/betB.eps %guardar foto per l'informe
%---------------------------------------
figure;
surf(sigma,flux,CL);
savefig('./figures/CL.fig')
close
% plt = Plot('figure', 'true');
    %Dibuixar bonic
plt = Plot('./figures/CL.fig', 'true');
% plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = '1/Solidesa, S/C';
plt.YLabel = 'Flux, \Psi';
plt.ZLabel = 'C_{L}';
plt.LineWidth(:) = 2;
plt.Title = 'C_{L} en funci� de S/C i \Psi';

print -depsc2 figures/parametres/CL.eps %guardar foto per l'informe
%---------------------------------------

figure;
surf(sigma,flux,CD);
savefig('./figures/CD.fig')
close
% plt = Plot('figure', 'true');
    %Dibuixar bonic
plt = Plot('./figures/CD.fig', 'true');
% plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = '1/Solidesa, S/C';
plt.YLabel = 'Flux, \Psi';
plt.ZLabel = 'C_{D}';
plt.LineWidth(:) = 2;
plt.Title = 'C_{D} en funci� de S/C i \Psi';

print -depsc2 figures/parametres/CD.eps %guardar foto per l'informe
%---------------------------------------

figure;
surf(sigma,flux,ETAesg);
savefig('./figures/ETAesg.fig')
close
% plt = Plot('figure', 'true');
    %Dibuixar bonic
plt = Plot('./figures/ETAesg.fig', 'true');
% plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = '1/Solidesa, S/C';
plt.YLabel = 'Flux, \Psi';
plt.ZLabel = '\eta_{esgla�}';
plt.LineWidth(:) = 2;
plt.Title = '\eta_{esgla�} en funci� de S/C i \Psi';

print -depsc2 figures/parametres/ETAesg.eps %guardar foto per l'informe
%---------------------------------------

figure;
surf(sigma,flux,Vz);
savefig('./figures/Vz.fig')
close
% plt = Plot('figure', 'true');
    %Dibuixar bonic
plt = Plot('./figures/Vz.fig', 'true');
% plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = '1/Solidesa, S/C';
plt.YLabel = 'Flux, \Psi';
plt.ZLabel = 'V_{z}';
plt.LineWidth(:) = 2;
plt.Title = 'V_{z} en funci� de S/C i \Psi';

print -depsc2 figures/parametres/Vz.eps %guardar foto per l'informe
%---------------------------------------

figure;
surf(sigma,flux,U);
view(-33,33)
savefig('./figures/U.fig')
close
% plt = Plot('figure', 'true');
    %Dibuixar bonic
plt = Plot('./figures/U.fig', 'true');
% plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = '1/Solidesa, S/C';
plt.YLabel = 'Flux, \Psi';
plt.ZLabel = 'U';
plt.LineWidth(:) = 2;
plt.Title = 'U en funci� de S/C i \Psi';

print -depsc2 figures/parametres/U.eps %guardar foto per l'informe
%---------------------------------------
figure;
surf(sigma,flux,TAUesg);
view(31,30)
savefig('./figures/TAUesg.fig')
close
% plt = Plot('figure', 'true');
    %Dibuixar bonic
plt = Plot('./figures/TAUesg.fig', 'true');
% plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = '1/Solidesa, S/C';
plt.YLabel = 'Flux, \Psi';
plt.ZLabel = '\tau_{esgla�}';
plt.LineWidth(:) = 2;
plt.Title = '\tau_{esgla�} en funci� de S/C i \Psi';

print -depsc2 figures/parametres/TAUesg.eps %guardar foto per l'informe
%---------------------------------------

figure;
surf(sigma,flux,rire);
view(146,31)
savefig('./figures/rire.fig')
close
% plt = Plot('figure', 'true');
    %Dibuixar bonic
plt = Plot('./figures/rire.fig', 'true');
% plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = '1/Solidesa, S/C';
plt.YLabel = 'Flux, \Psi';
plt.ZLabel = 'r_{i}/r_{e}';
plt.LineWidth(:) = 2;
plt.Title = 'r_{i}/r_{e} en funci� de S/C i \Psi';

print -depsc2 figures/parametres/rire.eps %guardar foto per l'informe
%---------------------------------------

figure;
surf(sigma,flux,re);
view(30.74,34.64)
savefig('./figures/re.fig')
close
% plt = Plot('figure', 'true');
    %Dibuixar bonic
plt = Plot('./figures/re.fig', 'true');
% plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = '1/Solidesa, S/C';
plt.YLabel = 'Flux, \Psi';
plt.ZLabel = 'r_{e}';
plt.LineWidth(:) = 2;
plt.Title = 'r_{e} en funci� de S/C i \Psi';

print -depsc2 figures/parametres/re.eps %guardar foto per l'informe
%---------------------------------------

figure;
surf(sigma,flux,ri);
view(-31,32)
savefig('./figures/ri.fig')
close
% plt = Plot('figure', 'true');
    %Dibuixar bonic
plt = Plot('./figures/ri.fig', 'true');
% plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = '1/Solidesa, S/C';
plt.YLabel = 'Flux, \Psi';
plt.ZLabel = 'r_{i}';
plt.LineWidth(:) = 2;
plt.Title = 'r_{i} en funci� de S/C i \Psi';

print -depsc2 figures/parametres/ri.eps %guardar foto per l'informe
%---------------------------------------


figure;
surf(sigma,flux,rm);
view(-25.9,46.8)
savefig('./figures/rm.fig')
close
% plt = Plot('figure', 'true');
    %Dibuixar bonic
plt = Plot('./figures/rm.fig', 'true');
% plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = '1/Solidesa, S/C';
plt.YLabel = 'Flux, \Psi';
plt.ZLabel = 'r_{m}';
plt.LineWidth(:) = 2;
plt.Title = 'r_{m} en funci� de S/C i \Psi';

print -depsc2 figures/parametres/rm.eps %guardar foto per l'informe
%---------------------------------------

figure;
surf(sigma,flux,h);
savefig('./figures/h.fig')
close
% plt = Plot('figure', 'true');
    %Dibuixar bonic
plt = Plot('./figures/h.fig', 'true');
% plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = '1/Solidesa, S/C';
plt.YLabel = 'Flux, \Psi';
plt.ZLabel = 'h';
plt.LineWidth(:) = 2;
plt.Title = 'h en funci� de S/C i \Psi';

print -depsc2 figures/parametres/h.eps %guardar foto per l'informe
%---------------------------------------

figure;
surf(sigma,flux,RPM);
savefig('./figures/RPM.fig')
close
% plt = Plot('figure', 'true');
    %Dibuixar bonic
plt = Plot('./figures/RPM.fig', 'true');
% plt.XLim(2) = plt.XLim(2) + .15;
plt.BoxDim = [6 4];
plt.XLabel = '1/Solidesa, S/C';
plt.YLabel = 'Flux, \Psi';
plt.ZLabel = 'N(rpm)';
plt.LineWidth(:) = 2;
plt.Title = 'N en funci� de S/C i \Psi';

print -depsc2 figures/parametres/RPM.eps %guardar foto per l'informe
%---------------------------------------



