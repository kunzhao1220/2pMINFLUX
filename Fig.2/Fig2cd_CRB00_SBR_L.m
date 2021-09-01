clear all
close all
clc
NN=100;% photon number
K=4;
wavelength=[647,800];
fwhm_option=wavelength*0.51/1.4;
Lmax=100;
L=0:0.1:Lmax;
L0=50;
SBRL0=[3,2,4];  % !=0
numSBR0=numel(SBRL0);

line_width=1;
label={'test'};

for photon=1:2
    index=2*photon;
    fwhm=fwhm_option(photon);
    
    for ii=1:numSBR0
        if SBRL0(ii) < Inf
            SBR=(L.^index/(L0^index)).*exp(index*(log(2)/fwhm^2).*(L0^2-L.^2))*SBRL0(ii);
            CRB=L .* sqrt((1./SBR + 1).*(3/K./SBR + 1)/NN/2)  ...
                ./(index * (1 - (log(2)*L.^2 /fwhm^2)) ) ;
        else
            CRB = L .* sqrt(1/NN/2) ./(index * (1 - (log(2)*L.^2 /fwhm^2)) ) ;
        end
        if ii==1    c1='k';c2='r'; end
        if ii==2    c1='k--';c2='r--'; end
        if ii==3    c1='k-.';c2='r-.'; end
        if ii==4    c1='k:';c2='r:'; end
        if photon==2
            figure(1), plot(L,CRB,c2,'LineWidth',line_width), hold on
            figure(2), semilogy(L,SBR,c2,'LineWidth',line_width); hold on
        else
            figure(1), plot(L,CRB,c1,'LineWidth',line_width), hold on
            figure(2), semilogy(L,SBR,c1,'LineWidth',line_width); hold on
        end
        
        minCRB(photon,ii)=min(CRB);
        minCRB_index(photon,ii)=find(CRB==min(CRB));
        L_min(photon,ii)=L(minCRB_index(photon,ii));
        SBR_min(photon,ii)=SBR(minCRB_index(photon,ii));
        CRB_L50(photon,ii)=CRB(L==L0);
        
        name=sprintf('%dp, SBR_L_0=%.1f',photon, SBRL0(ii));
        label(ii+photon*numSBR0-numSBR0) = cellstr(name);
    end
end

figure(1)
title(['CRB versus L & SBR, L_0=',num2str(L0),'nm'])
axis([0 Lmax 1 5])
legend(label,'Location','north')
hold on, scatter(L_min(1,:),minCRB(1,:),'MarkerFaceColor','k');
hold on, scatter(L_min(2,:),minCRB(2,:),'MarkerFaceColor','r');
ylabel('CRB (nm)'), xlabel('L (nm)')
line([50 50],[0 10])
box off

figure(2)
title(['SBR vs L, L_0=',num2str(L0),'nm'])
axis([0 Lmax 0.5 10])
legend(label,'Location','northwest')
ylabel('SBR'), xlabel('L (nm)')
hold on, scatter(L_min(1,:),SBR_min(1,:),'MarkerFaceColor','k');
hold on, scatter(L_min(2,:),SBR_min(2,:),'MarkerFaceColor','r');
line([50 50],[0.01 100])
box off

Ratio_LminL50 = 1 - minCRB./CRB_L50
RatioMin_L50  = CRB_L50(1,:)./CRB_L50(2,:)
SBR_min
L_min