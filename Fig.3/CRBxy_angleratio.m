clear all
close all
clc
%% Definition
L=50;
angle=[0,pi/9,2*pi/9,pi/3];
% TCPlabel=[{'Triangle+Origin'}];
wavelength=[647,800];
NN=100;
SBR_L0=3;
L0=50;
for mm=1:length(angle)
    name=sprintf('angle = %.0f¡ã',angle(mm)*180/pi);
    label(mm)=cellstr(name);
end

r=0:0.1:L/2;
n_r=numel(r);
xx=zeros(length(angle),n_r);
yy=zeros(length(angle),n_r);
%% Donut position
DD=L0/2;
x0=0; y0=0;
times=2/3*pi*1;    x1=DD*cos(times); y1=DD*sin(times);
times=2/3*pi*2;    x2=DD*cos(times); y2=DD*sin(times);
times=2/3*pi*3;    x3=DD*cos(times); y3=DD*sin(times);
%times=2/3*pi*3;    x4=DD*cos(times); y4=DD*sin(times);
R1=[x1,y1]; R2=[x2,y2]; R3=[x3,y3];
%R4=[x4,y4];% for LMS


%% Calculation
for p2=0:1
    if p2==0
        wavelength=647;
        fwhm=wavelength/2/1.4;
        Dntxy = @(x,y) 4*exp(1)*log(2)*(x^2+y^2)/fwhm^2 .*exp(-4*log(2)*(x^2+y^2)/fwhm^2);
    else
        wavelength=800;
        fwhm=wavelength/2/1.4;
        Dntxy = @(x,y)( 4*exp(1)*log(2)*(x^2+y^2)/fwhm^2 .*exp(-4*log(2)*(x^2+y^2)/fwhm^2) )^2 ;
    end
    
    %% SBR
    I0 = @(x,y) Dntxy(x-x0,y-y0);
    I1 = @(x,y) Dntxy(x-x1,y-y1);
    I2 = @(x,y) Dntxy(x-x2,y-y2);
    I3 = @(x,y) Dntxy(x-x3,y-y3);
    %I4 = @(x,y) Dntxy(x-x4,y-y4);
    %I_sum=@(x,y) I0(x,y)+I1(x,y)+I2(x,y)+I3(x,y)+I4(x,y);%Square+Ori
    I_sum=@(x,y) I0(x,y)+I1(x,y)+I2(x,y)+I3(x,y); %Triangle+Ori
    %sbr=@(x,y,l,index) (I_sum(x,y)/(I0(0,0)+I1(0,0)+I2(0,0)+I3(0,0)+I4(0,0)))*(l.^index/(L0^index)).*exp(index*(log(2)/fwhm^2).*(L0^2-l.^2))*SBR_L0; %Square+Ori
    sbr=@(x,y,l,index)(I_sum(x,y)/(I0(0,0)+I1(0,0)+I2(0,0)+I3(0,0)))*(l.^index/(L0^index)).*exp(index*(log(2)/fwhm^2).*(L0^2-l.^2))*SBR_L0; %Triangle+Ori
    
    %% Calculation
    
    for mm=1:length(angle)
        [xx(mm,:),yy(mm,:)]=pol2cart(angle(mm),r);
        for jj=1:n_r
            if p2==0
                index=2*(p2+1);
                SBR(mm,jj,p2+1)=sbr(xx(mm,jj),yy(mm,jj),L,index);
                % for jj=1:n_r
                %         index=2*p12(mm);
                %         SBR(mm,p2+1)=(L(mm).^index/(L0^index)).*exp(index*(log(2)/fwhm^2).*(L0^2-L(mm).^2))*SBR_L0;
                CRBcut1D_1p(mm,jj) = CRB_1P_T4(xx(mm,jj),yy(mm,jj),fwhm,L,NN,SBR(mm,jj,p2+1),4);
                %CRBcut1D_1p(mm,:) = CRB_1P_S5_SBR(xx(mm,:),yy(mm,:),fwhm,L0,NN,SBR_L0,5);
                %  end
            else
                index=2*(p2+1);
                SBR(mm,jj,p2+1)=sbr(xx(mm,jj),yy(mm,jj),L,index);
                %     for jj=1:n_r
                %         index=2*p12(mm);
                %         SBR(mm,p2+1)=(L(mm).^index/(L0^index)).*exp(index*(log(2)/fwhm^2).*(L0^2-L(mm).^2))*SBR_L0;
                CRBcut1D_2p(mm,jj) = CRB_2P_T4(xx(mm,jj),yy(mm,jj),fwhm,L,NN,SBR(mm,jj,p2+1),4);
                %CRBcut1D_2p(mm,:) = CRB_2P_S5_SBR(xx(mm,:),yy(mm,:),fwhm,L0,NN,SBR_L0,5);
            end
        end
    end
end
Ratio=CRBcut1D_1p./CRBcut1D_2p;

%% plotting
plot(r,Ratio(:,:));
title(['CRB Ratio, N=',num2str(NN)]);%,', SBR=',num2str(SBR)])
legend(label,'Location','northeast')
%hold on ,scatter([-L/2,L/2],[zeros(1,2*length(L))],'ro')
ylim([0 max(max(Ratio))]), xlim([0 max(r)])
ylabel('CRB Ratio'), xlabel('y (nm)'), box off
