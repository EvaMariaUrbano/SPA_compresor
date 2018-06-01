function [ FLUX ] = interpolarPunt( TAU, TAUa, TAUb, fluxB, fluxA )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
FLUX = fluxA + (TAU - TAUa)*(fluxB-fluxA)/(TAUb-TAUa);

end

