close all

x=-400:1:400;
index1=find(x==25);

fwhm=647*0.51/1.4; 
I1p = (4*exp(1)*log(2)*(x.^2/fwhm^2).*exp(-4*log(2)*x.^2/fwhm^2) ) ;

fwhm=800*0.51/1.4;
I2pEx = (4*exp(1)*log(2)*(x.^2/fwhm^2).*exp(-4*log(2)*x.^2/fwhm^2) ) ;
I2pEm=I2pEx.^2;

%% Figure RR1
subplot(2,2,1)
colororder({'r','k'})
yyaxis left
plot(x,I2pEx,'linewidth',1);
ylabel('2p intensity (I/Ib)')
yyaxis right
plot(x,I1p,'linewidth',1);
ylabel('1p intensity (I/Ia)')
hold on, scatter(x(index1),I2pEx(index1),'r','filled')
hold on, scatter(x(index1),I1p(index1),'k','filled')
box off
xlabel('r (nm)')
title('(a)    donut excitation 2p VS 1p')
legend('2p Ex','1p Ex','location','northeast')
xlim([-400 400])

subplot(2,2,2), plot(x,I2pEm,'r',x,I1p,'k','linewidth',1);
hold on, scatter(x(index1),I2pEx(index1)^2,'r','filled')
hold on, scatter(x(index1),I1p(index1),'k','filled')
box off
xlabel('r (nm)')
ylabel('Normalized intensity')
title('(b)    donut emission 2p VS 1p')
legend('2p Em','1p Em','location','northeast')
xlim([-400 400])

subplot(2,2,3)
yyaxis left
plot(x,I2pEx/I2pEx(index1),'linewidth',1);
ylabel('2p intensity (I/Ib)')
ylim([0 20])
yyaxis right
plot(x,I1p/I1p(index1),'linewidth',1);
ylabel('1p intensity (I/Ia)')
ylim([0 20])
hold on, scatter(x(index1),1,'k','filled')
hold on, scatter(x(index1),1,'r','filled')
box off
xlabel('r (nm)')
title('(c)    MINFLUX excitation 2p VS 1p')
legend('2p Ex','1p Ex','location','northeast')
xlim([-400 400])

subplot(2,2,4), plot(x,I2pEm/I2pEx(index1)^2,'r--',x,I1p/I1p(index1),'k--','linewidth',1);
hold on, scatter(x(index1),1,'k','filled')
hold on, scatter(x(index1),1,'r','filled')
box off
xlabel('r (nm)')
ylabel('Normalized intensity')
title('(d)    MINFLUX emission 2p VS 1p')
legend('2p Em','1p Em','location','northeast')
xlim([-400 400])
