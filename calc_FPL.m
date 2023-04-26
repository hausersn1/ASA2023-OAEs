function [S] = calc_FPL(f, complexP, C)

% C contains calibration parameters
% f (Hz) and complexP contains stimulus data

% interpolate first
fi = C.freq.*1000; %in Hz
VtoFPL = interp1(fi,C.Pfor,f,'linear','extrap');

S.P_fpl = complexP .* VtoFPL'; 
S.f = f;

return