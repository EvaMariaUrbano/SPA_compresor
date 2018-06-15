% para el caso de ocho etapas
clc; clear all;
%% Constantes

g=9.81; %m/s2
gamma=1.4;
Rgas=286.8; %J/kgK
Cp=1004; %J/kgK

%% Datos

Rreac=0.5;
Pat=1*g*10000;
Tat=288; %K
m = 21.38;  %kg/s
M=0.8;
lambda=0.7;
sigmAl=(39/4)*10^6; %kg/m2
rho_al=2800; %kg/m2

%% Datos calculados

betA = 0.9511; %rad
betB = 0.5393; %rad
betM = 0.7854;
alphaA = betB;

rm  = 0.221;
rire = 0.5890;
Wesg = 37500;  %J/kg
w = 1.3828e+03; %rad/s
Vz = 153.5685; %m/s

%% Calculos

cp = 1-(cos(betA)/cos(betB))^2;
pb_pa = 1 + 0.5*gamma*cp*M^2;
% salto de presiones en una etapa
pi_et = pb_pa^2;

et = 8;
h = zeros(et,2); %1a columna -> seccion a. 2a columna -> seccion b.

for i=1:et
    if i == 1
        % Entramos con M pero luego seguimos el loop con Va que es cte
        Ta = Tat/(1+((gamma-1)/2)*M^2);
        Pa = Pat/((1+((gamma-1)/2)*M^2)^(gamma/(gamma-1)));
        a = sqrt(gamma*Rgas*Ta);
        Va = M*a;
        ua = Va*cos(alphaA);
    else
        % Va se mantiene de la etapa anterior
        Ta = Tat - Va^2/(2*Cp);
        % Tat viene de la etapa anterior Tat_i = Tbt_(i-1)
        Pa = Pbt*(Ta/Tat)^(gamma/(gamma-1));
    end
    rhoa = Pa/(Rgas*Ta);
    % forumla de gasto
    re=sqrt(m/(pi*(1-rire^2)*Vz*rhoa));
    ri = rire*re;
    h(i,1) = re - ri;

    %----------------------------------------- 
    Tbt = Tat + Wesg/Cp;
    Pbt = Pat*pb_pa;
    % Con triangulo obtenemos Vb
    % se supone que la geometria se repite y que rm es cte.
    ub=ua;
    vbR = ub*tan(betB);
    vb = w*rm - vbR;
    Vb = sqrt(ub^2+vb^2);
    Tb = Tbt - Vb^2/(2*Cp);
    Pb = Pbt*(Tb/Tbt)^(gamma/(gamma-1));
    rhob = Pb/(Rgas*Tb);
    % forumla de gasto
    re=sqrt(m/(pi*(1-rire^2)*Vz*rhob));
    ri = rire*re;
    h(i,2) = re - ri;

    % para siguiente iteracion
    Tat = Tbt;
    Pat = Pbt*pb_pa;
end

C_h = 1/2.5;
C_h = linspace(C_h,1.25*C_h,et);
L = 0;
for i=1:et
    h_r = (h(i,1)+h(i,2))/2;
    if i == et
        h_e = (h(i,2)+h(i,1))/2;
    else
        h_e = (h(i,2)+h(i+1,1))/2;
    end
    L_stage = longitud(C_h(i),h_r,h_e,et,betM);
    L = L + L_stage;   
end
h_vect = zeros(2*et,1);
for i=1:et
    h_vect(i*2-1) = h(i,1);
    h_vect(i*2) = h(i,2);
end
L_vect = linspace(0,L,size(h_vect,1));
Zero_line = linspace(0,0,size(h_vect,1));
vert_line = zeros(5,et);

% aux = linspace(0,L,et);
% for i = 1 : et
%     vert_line(:,i) = linspace(-h(i,1),h(i,1),5);
% end
% for j = 1:5
%     L_vect_u(j,:) = aux;
% end
plot(L_vect,h_vect,L_vect,-h_vect,L_vect,Zero_line)
    