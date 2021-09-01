close all
clear x
syms x
N=50;
% sig=diff(normcdf([-1 1]));
fwhm=[647,800]*0.51/1.4;
sig=0.95;
[l,lci] = poissfit(N,1-sig);
r0 = 25; %nm
% CI_r0_p1 = sqrt(lci)/sqrt(N)*r0
% CI_r0_p2 = nthroot(lci,4)/nthroot(N,4)*r0
Dnt1p = @(x)( (x^2)/fwhm(1)^2 .*exp(-4*log(2)*(x^2)/fwhm(1)^2) ) ;
Dnt2p = @(x)( (x^2)/fwhm(2)^2 .*exp(-4*log(2)*(x^2)/fwhm(2)^2) )^2;
for ii=1:2
    CI_r0_p1(ii) = vpasolve(Dnt1p(x)/Dnt1p(r0)==lci(ii)/N,x,r0);
    CI_r0_p2(ii) = vpasolve(Dnt2p(x)/Dnt2p(r0)==lci(ii)/N,x,r0);
end
CI_r0_p1=double(CI_r0_p1);
CI_r0_p2=double(CI_r0_p2);


xmax=75;
x_3points=[xmax,xmax/2,0];
x= -xmax:(1/2):xmax;
bg0=0;
bg=1/(1/bg0-1); % solve bg0=bg/(1+bg)

fwhm=647*0.51/1.4; I1p =   (4*exp(1)*log(2)*(x.^2/fwhm^2).*exp(-4*log(2)*x.^2/fwhm^2) + bg) /(1+bg) ;
fwhm=800*0.51/1.4; I2p = ( (4*exp(1)*log(2)*(x.^2/fwhm^2).*exp(-4*log(2)*x.^2/fwhm^2) ).^2 + bg) /(1+bg);
I1p=I1p/I1p(x==r0)*N;
I2p=I2p/I2p(x==r0)*N;
index1=(x==round(CI_r0_p1(1)) | (x==round(CI_r0_p1(2))));
index2=(x==round(CI_r0_p2(1)) | (x==round(CI_r0_p2(2))));

plot(x,I1p,'r',x,I2p,'b','linewidth',1);
legend('single-photon (647nm)','two-photon (800nm)')
hold on, scatter([-CI_r0_p1,CI_r0_p1],[0 0 0 0],'r.')
hold on, scatter([-CI_r0_p2,CI_r0_p2],[0 0 0 0],'b.')
hold on, scatter([0 0 0 0],[-lci;lci],'.')
box off
ylim([0 N*1.75])
xlabel('r (nm)')
ylabel('lambda')
% title(['fluorescence difference'])

% x = [-20:.01:20];
% norm = N-20+500*normpdf(x,0,10);
% hold on, plot(x,norm)

x1p=[-CI_r0_p1,CI_r0_p1];
ymax=[lci(1),lci(2),lci(1),lci(2)];
for ii=1:4
     line([x1p(ii),x1p(ii)],[0 ymax(ii)],'color','red')
end
x2p=[-CI_r0_p2,CI_r0_p2];
for ii=1:4
     line([x2p(ii),x2p(ii)],[0 ymax(ii)],'color','blue')
end