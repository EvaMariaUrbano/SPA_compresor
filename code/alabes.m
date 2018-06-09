function [N] = alabes(h,r_m,sigma)
%Computes number of blades
C = c_h*h;
S = (1/sigma)*C;
N = 2*pi*r_m/S;
end

