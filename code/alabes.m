function [N] = alabes(h,r_m,S_C)
%Computes number of blades
% from class:
c_h = 1/2.5;
C = c_h*h;
S = S_C*C;
N = 2*pi*r_m/S;
end

