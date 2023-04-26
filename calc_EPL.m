function [S] = calc_EPL(f, complexP, C)
% calculate rthe inital outgoing OAE wave at the drum
% C contains calibration parameters
% D contains OAE data

if ~(numel(f) == numel(C.Rs)) %need to interpolate first
    fi = 1000.*linspace(0,C.SamplingRate/2,numel(C.Rs)); %in kHz
    Rec = interp1(fi,C.Rec,f,'linear','extrap');
    Rs = interp1(fi,C.Rs,f,'linear','extrap');
else
    Rec = C.Rec;
    Rs = C.Rs;
end
% calculate the delay
S.tec_fres = 1./(2*C.fres);
delay = exp(-i*2*pi*f*S.tec_fres);

S.P_epl = (complexP.*(1-Rs.*Rec)')./(delay.*(1+Rs))'; %EPL

S.f = f;

return



