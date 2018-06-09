function [N] = alabes(h,r_m,sigma)
%Computes number of blades
% from teacher:
c_h = 1/2.5;
C = c_h*h;
S = (1/sigma)*C;
N = 2*pi*r_m/S;
end

