function L = longitud(c_h,hr,he,et,betM)
%computes length of stage
L = 0;
Cr = c_h*hr;
Ce = c_h*he;

if et == 1
    % IGV
    L = L+1.25*Cr;
end

L = L + 0.4*Cr+Cr*cos(betM)+0.25*Cr+Ce*cos(betM);
end

