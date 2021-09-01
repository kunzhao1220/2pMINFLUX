for photon=1
    if photon==1
        res=res1p;
        tt=['MINFLUX MLE, N=',num2str(N12(1))];
    else
        res=res2p;
        tt=['2p-MINFLUX MLE, N=',num2str(N12(2))];
    end
    
    [N,Xedges,Yedges] = histcounts2(res(:,1),res(:,2),[-12:(0.5):12],[-12:(0.5):12]);
    N2 = imresize(N,5,'nearest');
    % nn=size(N2);M=zeros([nn,3]);M(:,:,1)=N2;M(:,:,3)=N2;M2=M/max(M(:)); figure,imshow(M2,[])
    
    figure,imagesc(N2),colormap hot, colorbar;
    line([185 235],[230 230],'linewidth',5,'color','w')  % scale bar 5nm  % 0.1nm/pixel
    axis off
    axis equal
    hold on; scatter(120.5-10*ro(1,:),120.5+10*ro(2,:),'ro','MarkerFaceColor','b')
    title(tt)
    
end