clear all
close all
clc


%% Definition
L0=50; %nm
N12=[100,100];  % photon number for 1p and 2p
SBR_L0=3;
L=50;
K=2;
z = sym('z');
CRBL100nm=3.3/sqrt(2);
%% Donut position
DD=L0/2;
z1=DD; z2=-DD;
R1=z1; R2=z2;  % for LMS
% [x0 x1 x2 x3;y0 y1 y2 y3]'
%% Simulation
for photon=1:2  % 1p or 2p
    NN=N12(photon);
    if photon==1
        wavelength=647;
        fwhm=wavelength/2/1.4;
        index=2*photon;
        Dntz = @(z) z^2;
    else
        wavelength=800;
        fwhm=wavelength/2/1.4;
        index=2*photon;
        Dntz = @(z) z^4;
    end
    
    I1 = @(z) Dntz(z-z1);
    I2 = @(z) Dntz(z-z2);
    I_sum=@(z) I1(z)+I2(z);
    
    %% SBR Definition
    SBR=@(z,l,index) (I_sum(z)/(I1(0)+I2(0)))*(l.^index/(L0^index)).*exp(index*(log(2)/fwhm^2).*(L0^2-l.^2))*SBR_L0;
    
    
    %% MLE Definition
    P1 = @(z,l,index) (SBR(z,l,index)/(SBR(z,l,index)+1))*I1(z)/I_sum(z)+1/(K*(SBR(z,l,index)+1));
    P2 = @(z,l,index) (SBR(z,l,index)/(SBR(z,l,index)+1))*I2(z)/I_sum(z)+1/(K*(SBR(z,l,index)+1));
    lnp1z=@(z,l,index) diff(P1(z,l,index),z)/P1(z,l,index);  
    lnp2z=@(z,l,index) diff(P2(z,l,index),z)/P2(z,l,index);  
    
    
    
    %% Simulation
    
    ro=[-5 0 5];
    
    times=100; % repeat
    res=zeros(times,length(ro));  %r estimate
    opt=optimset('Display','off');
    
    res_avg_r0=[0];
    prec_z_r0=[0];
    
    for m=1:length(ro)
        tic
%         xx=ro(1,m);
%         yy=ro(2,m);
%         r0=ro(:,m)';
        
%         p0 = P0(xx,yy);  p1 = P1(xx,yy);    p2 = P2(xx,yy);    p3 = P3(xx,yy);
        
        parfor jj=1:times
            zz=normrnd(0,CRBL100nm/photon);
           p1=P1(zz,L,index); p2=P2(zz,L,index); 
            
            R = mnrnd(NN,[p1,p2]);          % Mulitnomial
            n1=R(1); n2=R(2);
            
            func = matlabFunction(n1*lnp1z(z,L,index) + n2*lnp2z(z,L,index), 'Vars',{z});
            
            r_LMS = @(x,y) (x*R1 + y*R2) / NN  / (-1 + L^2*log(2)/fwhm^2) / photon;
            solution = fsolve(func, r_LMS(n1,n2),opt); 
            % solution = fsolve(func, r0,opt);   % solve
            res(jj,m)=solution-zz+ro(m);
        end
       res_avg(m,:) = mean(res(:,m));
       prec_z(m,:) = std(res(:,m));  % std (precision)
        toc
    end
    prec_z_all = mean(prec_z_r0);
    
    resx=normrnd(0,1/photon,[times*length(ro),1]);
    scale=5;
    step=0.5;
    [N,Xedges,Yedges] = histcounts2(res(:),resx,[-10:step:10],[-5:step:5]);
    N2 = imresize(N,scale,'nearest');
    
    
    % muticolor %  nn=size(N2);M=zeros([nn,3]);M(:,:,1)=N2;M(:,:,3)=N2;M2=M/max(M(:)); figure,imshow(M2,[])
    
    figure,imagesc(N2),colormap hot, colorbar;
    line([185 235],[230 230],'linewidth',5,'color','w')  % scale bar 5nm  % 0.1nm/pixel
    axis off
    axis equal
    hold on; scatter(scale*(scale/step)*ones(1,3),100-(scale/step)*ro(:),'ro','MarkerFaceColor','b')
    
    if photon==1
        title(   ['MINFLUX z-axis MLE, N=',num2str(N12(1))])
        res1p=res; % data storage
    else
        title(['2p-MINFLUX z-axis MLE, N=',num2str(N12(2))])
        res2p=res; % data storage
    end
end