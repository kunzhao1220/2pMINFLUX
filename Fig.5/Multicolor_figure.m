index=reshape(1:times*size(ro,2),times,size(ro,2));
color1=int8([3,5,8,9]);
color2=int8([1,2,4,6,7]);
res1=res2p(index(:,color1),:);
res2=res2p(index(:,color2),:);
[N1,Xedges,Yedges] = histcounts2(res1(:,1),res1(:,2),[-12:(0.5):12],[-12:(0.5):12]);
[N2,Xedges,Yedges] = histcounts2(res2(:,1),res2(:,2),[-12:(0.5):12],[-12:(0.5):12]);

nn=size(N1);M=zeros([nn,3]);
M(:,:,1)=N2/max(N1(:));
M(:,:,2)=N1/max(N1(:));
M2 = imresize(M,5,'nearest');
figure,imshow(M2)
line([185 235],[230 230],'linewidth',5,'color','w')  % scale bar 5nm  % 0.1nm/pixel
axis off
axis equal
hold on; scatter(120.5-10*ro(1,:),120.5+10*ro(2,:),'ro','MarkerFaceColor','b')
title(['2p-mc-MINFLUX MLE, N=',num2str(N12(2))]);