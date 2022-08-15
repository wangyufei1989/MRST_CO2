function G=orth2radial(G)
% for i=1:G.cartDims(3)
% G.faces.areas((i-1)*(G.cartDims(1)+1)+1:i*(G.cartDims(1)+1))=2.*pi.*G.faces.centroids((i-1)*(G.cartDims(1)+1)+1:i*(G.cartDims(1)+1),1);
% 
% end
if G.cartDims(1)>20
    k=10;
    
else
    k=2;
end
xtot=G.faces.centroids(G.cartDims(1)+1,1);
dx0=2*xtot/G.cartDims(1)/(G.cartDims(1)/k);%increase linearly from dx0 to (G.CartDims(1)/k-1)*dx0
dxend=(G.cartDims(1)/k-1)*dx0;
dX=linspace( dx0,dxend,G.cartDims(1));
for i=1:G.cartDims(1)
x(i)=sum(dX(1:i))-0.5*dX(i);
end
G.cells.centroids(:,1)=repmat(x',G.cartDims(3),1);

fcx(1)=0;
for i=2:G.cartDims(1)+1
fcx(i)=sum(dX(1:i-1));
end
G.faces.centroids(1:(G.cartDims(1)+1)*G.cartDims(3),1)=repmat(fcx',G.cartDims(3),1);
G.nodes.coords(:,1)=repmat(fcx',(G.cartDims(3)+1)*(G.cartDims(2)+1),1);
G.faces.centroids((G.cartDims(1)+1).*G.cartDims(2).*G.cartDims(3)+(G.cartDims(2)+1).*G.cartDims(1).*G.cartDims(3)+1:G.faces.num)=...
    repmat(G.cells.centroids(1:G.cartDims(1),1),G.cartDims(3)+1,1);


dZ=(G.faces.centroids(G.cells.faces(6:6:end,1),3)-G.faces.centroids(G.cells.faces(5:6:end-1,1),3));
for i=1:G.cartDims(3)
    dz((i-1)*(G.cartDims(1)+1)+1:i*(G.cartDims(1)+1),1)=dZ(i*G.cartDims(1)); 
end
G.faces.areas(1:(G.cartDims(1)+1)*G.cartDims(2)*G.cartDims(3))=2.*pi.*G.faces.centroids(1:(G.cartDims(1)+1)*G.cartDims(2)*G.cartDims(3),1)...
    .*dz;
Dr=G.faces.centroids(2:G.cartDims(1)+1,1)-G.faces.centroids(1:G.cartDims(1),1);
for i=1:G.cartDims(3)+1
  dr((i-1)*(G.cartDims(1))+1:i*(G.cartDims(1)),1)=Dr; 
end

G.faces.areas((G.cartDims(1)+1).*G.cartDims(2).*G.cartDims(3)+(G.cartDims(2)+1).*G.cartDims(1).*G.cartDims(3)+1:G.faces.num)=...
    2.*pi.*G.faces.centroids((G.cartDims(1)+1).*G.cartDims(2).*G.cartDims(3)+(G.cartDims(2)+1).*G.cartDims(1).*G.cartDims(3)+1:G.faces.num,1).*dr;
A=G.faces.areas((G.cartDims(1)+1).*G.cartDims(2).*G.cartDims(3)+(G.cartDims(2)+1).*G.cartDims(1).*G.cartDims(3)+1:(G.cartDims(1)+1).*G.cartDims(2).*G.cartDims(3)+(G.cartDims(2)+1).*G.cartDims(1).*G.cartDims(3)+G.cartDims(1));
for i=1:G.cartDims(3)
  a((i-1)*(G.cartDims(1))+1:i*(G.cartDims(1)),1)=A; 
end
G.cells.volumes=a.*(G.faces.centroids(G.cells.faces(6:6:end,1),3)-G.faces.centroids(G.cells.faces(5:6:end-1,1),3));
G.faces.normals(1:(G.cartDims(1)+1)*G.cartDims(3),1)=G.faces.areas(1:(G.cartDims(1)+1)*G.cartDims(3));
fz=(G.cartDims(1)+1)*G.cartDims(3)*G.cartDims(2)+(G.cartDims(2)+1)*G.cartDims(3)*G.cartDims(1)+1:...
    (G.cartDims(1)+1)*G.cartDims(3)*G.cartDims(2)+(G.cartDims(2)+1)*G.cartDims(3)*G.cartDims(1)+(G.cartDims(3)+1)*G.cartDims(2)*G.cartDims(1);
G.faces.normals(fz,3)=G.faces.areas(fz);
end