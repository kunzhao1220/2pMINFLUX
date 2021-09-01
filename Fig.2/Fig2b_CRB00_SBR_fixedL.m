close all

NN=100; % photon number
K=4;
wavelength=[647,800];
fwhm_option=wavelength*0.51/1.4;
L=50;
SBR = 10.^(-0.75:0.05:1);
CRB = [SBR;SBR];

line_width=1;
label=[{'1p'},{'2p'}];

for photon=1:2
    index=2*photon;
    fwhm=fwhm_option(photon);
    
    CRB(photon,:) = L .* sqrt((1./SBR + 1).*(3/K./SBR + 1)/NN/2) ./(index * (1 - (log(2)*L.^2 /fwhm^2)) ) ;
%     CRB(photon,:) = sqrt((1./SBR + 1).*(3/K./SBR + 1)) ;
    
    plot(SBR,CRB(photon,:)), hold on
%     semilogx(SBR,CRB(photon,:)), hold on
end

axis([0 10 0 10])
title('CRB vs SBR_L_0')
line([1 1],[0 10],'linestyle','--')
box off
ylabel('CRB (nm)'), xlabel('SBR')
legend(label)