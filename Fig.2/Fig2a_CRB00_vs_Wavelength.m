close all
clear all

wavelength=[400:10:700;700:10:1000];
fwhm=wavelength*0.51/1.4;

NN=100; % photon number
L=50;
SBRall=[2,3,4];
K=4;
CRB=fwhm;
% cc={{'--'},{'-'}};

%% 1p
for photon=1
figure
jj=1;
    for ii=1:length(SBRall)
        SBR=SBRall(ii);
        index=2*photon;
        if SBR < Inf
            CRB(photon,:) = L .* sqrt((1./SBR + 1).*(3/K./SBR + 1)./NN/2)  ...
                ./(index * (1 - (log(2)*L^2 ./fwhm(photon,:).^2)) ) ;
        else
            CRB(photon,:) = L .* sqrt(1./NN/2) ./(index * (1 - (log(2)*L^2 ./fwhm(photon,:).^2)) ) ;
        end
        plot(wavelength(photon,:),CRB(photon,:),'--'),hold on
        name=sprintf('SBR=%.0f',SBR);
        label(jj)=cellstr(name);
        jj=jj+1;
    end
    legend(label,'Location','northeast')
    title(sprintf('CRB vs Wavelength, %dp',photon))
    xlabel('Wavelength (nm)')
    ylabel('CRB (nm)')
    box off
end

%% 2p
for photon=2
figure
jj=1;
    for ii=1:length(SBRall)
        SBR=SBRall(ii);
        index=2*photon;
        if SBR < Inf
            CRB(photon,:) = L .* sqrt((1./SBR + 1).*(3/K./SBR + 1)./NN/2)  ...
                ./(index * (1 - (log(2)*L^2 ./fwhm(photon,:).^2)) ) ;
        else
            CRB(photon,:) = L .* sqrt(1./NN/2) ./(index * (1 - (log(2)*L^2 ./fwhm(photon,:).^2)) ) ;
        end
        plot(wavelength(photon,:),CRB(photon,:)),hold on
        name=sprintf('SBR=%.0f',SBR);
        label(jj)=cellstr(name);
        jj=jj+1;
    end
    legend(label,'Location','northeast')
    title(sprintf('CRB vs Wavelength, %dp',photon))
    xlabel('Wavelength (nm)')
    ylabel('CRB (nm)')
    box off
end
figure(1); axis([400 700  1.8 2.8])
figure(2); axis([700 1000 1 2])