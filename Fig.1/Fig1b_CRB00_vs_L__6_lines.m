close all
clear all
clc
figure
x=0;
y=x;
L0=50;
K=4;
%wavelength=[647,800];
fwhm_option=[360,420];%wavelength/2/1.4;

NN=[100,400,1600];% photon number
L=40:0.1:200;
SBR_L0=3;
label={'test'};

jj=1;
for ii=1:length(NN)
    
    for photon=1:2
        index=2*photon;
        fwhm=fwhm_option(photon);
        
        if SBR_L0 < Inf
             SBR(ii,:,photon)=(L.^index/(L0^index)).*exp(index*(log(2)/fwhm^2).*(L0^2-L.^2))*SBR_L0;
            CRB(ii,:,photon) = L .* sqrt((1./SBR(ii,:,photon) + 1).*(3/K./SBR(ii,:,photon) + 1)./NN(ii)/2)  ...
                ./(index * (1 - (log(2)*L.^2 /fwhm^2)) ) ;
            else
            CRB(ii,:,photon) = L .* sqrt(1./NN(ii)/2) ./(index * (1 - (log(2)*L.^2 /fwhm^2)) ) ;
        end
       
        if ii==1    c1='k';c2='r'; end
        if ii==2    c1='k--';c2='r--'; end
        if ii==3    c1='k:';c2='r:'; end
        figure(1)
        if photon==2    plot(L,CRB(:,:,photon),c1),hold on
        else            plot(L,CRB(:,:,photon),c2),hold on
        end
        name1=sprintf('N=%.0f, %.0fp',NN(ii),photon);
        label1(jj)=cellstr(name1);
        jj=jj+1;
        min_CRBresult(photon,ii)=min(CRB(ii,:,photon));
        min_CRBposition(photon,ii)=find(CRB(ii,:,photon)==min(CRB(ii,:,photon)));
        L_result(photon,ii)=L(min_CRBposition(photon,ii));
       
    end
    name2=sprintf('N=%.0f',NN(ii));
    label2(ii)=cellstr(name2);
    RCRB(ii,:)=CRB(ii,:,1)./CRB(ii,:,2);
    figure(2), plot(L,RCRB(ii,:)),hold on;
     
end
figure(1)
title(['CRB vs L, SBR_L_0=',num2str(SBR_L0)])%, SBR=',num2str(SBR)])
axis([min(L) max(L) 0 10])
legend(label1,'Location','north')

% figure(2)
% title(['RCRB vs L'])%, SBR=',num2str(SBR)])
% axis([min(L) max(L) 0 3])
% legend(label2,'Location','north')
