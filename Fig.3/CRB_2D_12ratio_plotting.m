close all
%% Plotting
%colorbar data
type=find(max(max_num));
k=[min(min_num) max(max_num)];
%draw circle
factor=round(500/ix);
DD=factor*L0/ddd;
r=DD/2;
delta=ceil(factor*ix/2);
% x_donut=delta+r*[cos(1*pi/2),cos(2*pi/2),cos(3*pi/2),cos(2*pi)];
% y_donut=delta+r*[sin(1*pi/2),sin(2*pi/2),sin(3*pi/2),sin(2*pi)];

 x_donut=delta+r*[cos(2*pi/3),cos(4*pi/3),cos(6*pi/3)];
 y_donut=delta+r*[sin(2*pi/3),sin(4*pi/3),0];

for pp=1:3
    if pp<3
        figure, imagesc(N2(:,:,pp)),colormap jet,caxis(k), colorbar; %subplot(1,2,pp)
        if pp==1   title(['1p-MINFLUX CRB'])
        else       title(['2p-MINFLUX CRB'])
        end
    else
        figure, imagesc(ratio),colormap jet, colorbar
        title(['CRB Ratio'])
    end
    line([400 450],[470 470],'linewidth',5,'color','w')  % scale bar 5nm
    hold on
    scatter(x_donut,y_donut,'ro','MarkerFaceColor','b');
    hold on
    alpha=0:pi/40:2*pi;
    xxx=r*cos(alpha)+delta;
    yyy=r*sin(alpha)+delta;
    plot(xxx,yyy,'--w','LineWidth',2);
    axis equal;
    axis off
    % text(150,550,['N=',num2str(NN),', L=',num2str(L),', SBR=',num2str(roundn(SBR,-1))])
end
