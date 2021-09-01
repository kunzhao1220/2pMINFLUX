close all
clc


%% Definition
L0=50; %nm
N12=[400,100];  % photon number for 1p and 2p
SBR_L0=3;
K=4;
L=50;
x = sym('x', [1,2]);
CRB=3.3;
%% Donut position
DD=L/2;
x0=0; y0=0;
times=2/3*pi*1;    x1=DD*cos(times); y1=DD*sin(times);
times=2/3*pi*2;    x2=DD*cos(times); y2=DD*sin(times);
times=2/3*pi*3;    x3=DD*cos(times); y3=0;
R1=[x1,y1]; R2=[x2,y2]; R3=[x3,y3]; % for LMS
% [x0 x1 x2 x3;y0 y1 y2 y3]'
%% Simulation
for photon=1:2  % 2p or not
    NN=N12(photon);
    if photon==1
        wavelength=647;
        fwhm=wavelength/2/1.4;
        index=2*photon;
        Dntxy = @(x,y) 4*exp(1)*log(2)*(x^2+y^2)/fwhm^2 .*exp(-4*log(2)*(x^2+y^2)/fwhm^2);
    else
        wavelength=800;
        fwhm=wavelength/2/1.4;
        index=2*photon;
        Dntxy = @(x,y)( 4*exp(1)*log(2)*(x^2+y^2)/fwhm^2 .*exp(-4*log(2)*(x^2+y^2)/fwhm^2) )^2 ;
    end
    
    I0 = @(x,y) Dntxy(x-x0,y-y0);
    I1 = @(x,y) Dntxy(x-x1,y-y1);
    I2 = @(x,y) Dntxy(x-x2,y-y2);
    I3 = @(x,y) Dntxy(x-x3,y-y3);
    I_sum=@(x,y) I0(x,y)+I1(x,y)+I2(x,y)+I3(x,y);
    
    %% SBR Definition
   SBR=@(x,y,l,index) (I_sum(x,y)/(I0(0,0)+I1(0,0)+I2(0,0)+I3(0,0)))*(l.^index/(L0^index)).*exp(index*(log(2)/fwhm^2).*(L0^2-l.^2))*SBR_L0;
    
    %% MLE Definition
    P0 = @(x,y,l,index) (SBR(x,y,l,index)/(SBR(x,y,l,index)+1))*I0(x,y)/I_sum(x,y)+1/(K*(SBR(x,y,l,index)+1));
    P1 = @(x,y,l,index) (SBR(x,y,l,index)/(SBR(x,y,l,index)+1))*I1(x,y)/I_sum(x,y)+1/(K*(SBR(x,y,l,index)+1));
    P2 = @(x,y,l,index) (SBR(x,y,l,index)/(SBR(x,y,l,index)+1))*I2(x,y)/I_sum(x,y)+1/(K*(SBR(x,y,l,index)+1));
    P3 = @(x,y,l,index) (SBR(x,y,l,index)/(SBR(x,y,l,index)+1))*I3(x,y)/I_sum(x,y)+1/(K*(SBR(x,y,l,index)+1));
    lnp0x=@(x,y,l,index) diff(P0(x,y,l,index),x)/P0(x,y,l,index);  lnp0y=@(x,y,l,index) diff(P0(x,y,l,index),y)/P0(x,y,l,index);
    lnp1x=@(x,y,l,index) diff(P1(x,y,l,index),x)/P1(x,y,l,index);  lnp1y=@(x,y,l,index) diff(P1(x,y,l,index),y)/P1(x,y,l,index);
    lnp2x=@(x,y,l,index) diff(P2(x,y,l,index),x)/P2(x,y,l,index);  lnp2y=@(x,y,l,index) diff(P2(x,y,l,index),y)/P2(x,y,l,index);
    lnp3x=@(x,y,l,index) diff(P3(x,y,l,index),x)/P3(x,y,l,index);  lnp3y=@(x,y,l,index) diff(P3(x,y,l,index),y)/P3(x,y,l,index);
    
    
    
    %% Simulation
    
    xx=6/sqrt(2);   % minimum distace: xx*sqrt2
    ro=[0 0 0 2*xx -2*xx xx xx -xx -xx;...
        0 -2*xx 2*xx 0 0 xx -xx xx -xx];
    alfa = 25.8 * 3.1415926 / 180.0; % rotation angle
    for jj=1:length(ro(1,:))
        [TH,Rth] = cart2pol(ro(1,jj),ro(2,jj));
        ro(1,jj)=Rth*cos(TH+alfa);
        ro(2,jj)=Rth*sin(TH+alfa);
    end
    
    
    times=100; % repeat
    res=zeros(times*length(ro(1,:)),2);  %r estimate
    opt=optimset('Display','off');
    
    res_avg_r0=[0,0];
    prec_xy_r0=[0,0];
    
    for m=1:length(ro(1,:))
        tic
%         xx=ro(1,m);
%         yy=ro(2,m);
%         r0=ro(:,m)';
        
%         p0 = P0(xx,yy);  p1 = P1(xx,yy);    p2 = P2(xx,yy);    p3 = P3(xx,yy);
        
        parfor jj=(1+(m-1)*times):(m*times)
            Rxy=normrnd(0,CRB/photon);
            THxy=rand*2*pi;
            xx=Rxy*cos(THxy);
            yy=Rxy*sin(THxy);
            p0=P0(xx,yy,L,index); p1=P1(xx,yy,L,index); p2=P2(xx,yy,L,index); p3=P3(xx,yy,L,index);
            
            R = mnrnd(NN,[p0,p1,p2,p3]);          % Mulitnomial
            n0=R(1); n1=R(2); n2=R(3); n3=R(4);
            
            func = matlabFunction([n0*lnp0x(x(1),x(2),L,index) + n1*lnp1x(x(1),x(2),L,index) + n2*lnp2x(x(1),x(2),L,index) + n3*lnp3x(x(1),x(2),L,index), ...
                n0*lnp0y(x(1),x(2),L,index) + n1*lnp1y(x(1),x(2),L,index) + n2*lnp2y(x(1),x(2),L,index) + n3*lnp3y(x(1),x(2),L,index)], 'Vars',{[x(1), x(2)]});
            
            r_LMS = @(x,y,z) (x*R1 + y*R2 + z*R3) / NN  / (-1 + L^2*log(2)/fwhm^2) / photon; %original estimation of solution utilizing least mean squared to put in fsolve
            solution = fsolve(func, r_LMS(n1,n2,n3),opt); 
            % solution = fsolve(func, r0,opt);   % solve
            res(jj,:)=solution-[xx,yy]+ro(:,m)';
        end
        res_avg_r0(m,:) = mean(res((1+(m-1)*times):(m*times),:)); % average
        prec_xy_r0(m,:) = std(res((1+(m-1)*times):(m*times),:));  % std (precision)
        toc
    end
    prec_xy_all = mean(prec_xy_r0);
    
    [N,Xedges,Yedges] = histcounts2(res(:,1),res(:,2),[-12:(0.5):12],[-12:(0.5):12]);
    N2 = imresize(N,5,'nearest');
    % muticolor %  nn=size(N2);M=zeros([nn,3]);M(:,:,1)=N2;M(:,:,3)=N2;M2=M/max(M(:)); figure,imshow(M2,[])
    
    figure,imagesc(N2),colormap hot, colorbar;
    line([185 235],[230 230],'linewidth',5,'color','w')  % scale bar 5nm  % 0.1nm/pixel
    axis off
    axis equal
    hold on; scatter(120.5-10*ro(1,:),120.5+10*ro(2,:),'ro','MarkerFaceColor','b')  % ground truth
    
    if photon==1
        title(   ['MINFLUX MLE, N=',num2str(N12(1))])
        res1p=res; % data storage
    else
        title(['2p-MINFLUX MLE, N=',num2str(N12(2))])
        res2p=res; % data storage
    end
end