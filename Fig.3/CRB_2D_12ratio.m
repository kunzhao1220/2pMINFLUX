clc
clear all
close all
clear pp
%% Definition
L0=50; %nm
NN=100;% photon number
SBR_L0=3;
K=4;  
%% Donut position
DD=L0/2;
x0=0; y0=0;
times=2/3*pi*1;    x1=DD*cos(times); y1=DD*sin(times);
times=2/3*pi*2;    x2=DD*cos(times); y2=DD*sin(times);
times=2/3*pi*3;    x3=DD*cos(times); y3=0;
R1=[x1,y1]; R2=[x2,y2]; R3=[x3,y3]; % for LMS
% [x0 x1 x2 x3;y0 y1 y2 y3]'
%% Simulation
 for p2=0:1 % 2p or not
    if p2==0
        wavelength=647;
        fwhm=wavelength/2/1.4;
        Dntxy = @(x,y) 4*exp(1)*log(2)*(x^2+y^2)/fwhm^2 .*exp(-4*log(2)*(x^2+y^2)/fwhm^2);
    else
        wavelength=800;
        fwhm=wavelength/2/1.4;
        Dntxy = @(x,y)( 4*exp(1)*log(2)*(x^2+y^2)/fwhm^2 .*exp(-4*log(2)*(x^2+y^2)/fwhm^2) )^2 ;
    end
    
    I0 = @(x,y) Dntxy(x-x0,y-y0);
    I1 = @(x,y) Dntxy(x-x1,y-y1);
    I2 = @(x,y) Dntxy(x-x2,y-y2);
    I3 = @(x,y) Dntxy(x-x3,y-y3);
    I_sum=@(x,y) I0(x,y)+I1(x,y)+I2(x,y)+I3(x,y);
    sbr=@(x,y,l,index) (I_sum(x,y)/(I0(0,0)+I1(0,0)+I2(0,0)+I3(0,0)))*(l.^index/(L0^index)).*exp(index*(log(2)/fwhm^2).*(L0^2-l.^2))*SBR_L0;
    
    
  

    %% Calculation
    ddd=0.5;
    [xx,yy]=meshgrid(-25:ddd:25);
    [ix,iy]=size(xx);
    CRB2d=zeros(ix,iy,2);
    
    for ii=1:ix
        for jj=1:iy
            %if xx(ii,jj)^2+yy(ii,jj)^2 < R^2
                if p2==1
                    index=2*(p2+1);
                SBR(ii,jj,p2+1)=sbr(xx(ii,jj),yy(ii,jj),L0,index);
                CRB2d(ii,jj,p2+1) = CRB_2P_T4(xx(ii,jj),yy(ii,jj),fwhm,L0,NN,SBR(ii,jj,p2+1),K);
               % CRB2d(ii,jj,p2+1) = CRB_2P_S5_SBR(xx(ii,jj),yy(ii,jj),fwhm,L,NN,SBR,K);
                else
                    index=2*(p2+1);
                 SBR(ii,jj,p2+1)=sbr(xx(ii,jj),yy(ii,jj),L0,index);
                 CRB2d(ii,jj,p2+1) = CRB_1P_T4(xx(ii,jj),yy(ii,jj),fwhm,L0,NN,SBR(ii,jj,p2+1),K);
               % CRB2d(ii,jj,p2+1) = CRB_1P_S5_SBR(xx(ii,jj),yy(ii,jj),fwhm,L,NN,SBR,K);
                end
            %end
        end
    end
    max_num(p2+1)=max(max(CRB2d(:,:,p2+1)));
    min_num(p2+1)=min(min(CRB2d(:,:,p2+1)));
        
    N2(:,:,p2+1) = imresize(CRB2d(:,:,p2+1),int8(500/ix),'nearest');
    
end
ratio=N2(:,:,1)./N2(:,:,2);

%% Plotting
CRB_2D_12ratio_plotting