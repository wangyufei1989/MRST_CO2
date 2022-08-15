
function G=tensor_face_cell(G,nr,ny,nz)
%[xy,xz,yx,yz,zx,zy]=1,2,3,4,5,6
index_x=1:(nr+1)*nz*ny;
face_in=(G.faces.neighbors(:,1).*G.faces.neighbors(:,2))~=0;
face_in_x=index_x(face_in(index_x));

[neighborship, ~] = getNeighbourship(G, 'Topological', true);
face_cell_x=neighborship(face_in_x,:);
[cellNo, cf, ~] = getCellNoFaces(G);%cellNo records the location of faces for each cell. cf records the face marks of each cell.
f_xy1=cf((face_cell_x(:,1)-1).*6+3);
f_xy2=cf((face_cell_x(:,2)-1).*6+3);
f_xy3=cf((face_cell_x(:,1)-1).*6+4);
f_xy4=cf((face_cell_x(:,2)-1).*6+4);
c_xy1=neighborship(f_xy1,1);
c_xy2=neighborship(f_xy2,1);
c_xy3=neighborship(f_xy3,2);
c_xy4=neighborship(f_xy4,2);
face_in_in=c_xy1.*c_xy2.*c_xy3.*c_xy4~=0;
face_in_in_x=face_in_x(face_in_in);
G.faces4.c1=c_xy1(face_in_in);
G.faces4.c2=c_xy2(face_in_in);
G.faces4.c3=c_xy3(face_in_in);
G.faces4.c4=c_xy4(face_in_in);
G.faces4.x_dyz=(G.cells.centroids(c_xy3(face_in_in),2)-G.cells.centroids(c_xy1(face_in_in),2))./2;
G.faces4.describe=ones(size(G.faces4.c4));
G.faces4.n=face_in_in_x';
%%
index_x=1:(nr+1)*nz*ny;
face_in=(G.faces.neighbors(:,1).*G.faces.neighbors(:,2))~=0;
face_in_x=index_x(face_in(1:(nr+1)*nz*ny));

[neighborship, ~] = getNeighbourship(G, 'Topological', true);
face_cell_x=neighborship(face_in_x,:);
[cellNo, cf, ~] = getCellNoFaces(G);%cellNo records the location of faces for each cell. cf records the face marks of each cell.
f_xy1=cf((face_cell_x(:,1)-1).*6+5);
f_xy2=cf((face_cell_x(:,2)-1).*6+5);
f_xy3=cf((face_cell_x(:,1)-1).*6+6);
f_xy4=cf((face_cell_x(:,2)-1).*6+6);
c_xy1=neighborship(f_xy1,1);
c_xy2=neighborship(f_xy2,1);
c_xy3=neighborship(f_xy3,2);
c_xy4=neighborship(f_xy4,2);
face_in_in=c_xy1.*c_xy2.*c_xy3.*c_xy4~=0;
face_in_in_x=face_in_x(face_in_in);
G.faces4.c1=[G.faces4.c1;c_xy1(face_in_in)];
G.faces4.c2=[G.faces4.c2;c_xy2(face_in_in)];
G.faces4.c3=[G.faces4.c3;c_xy3(face_in_in)];
G.faces4.c4=[G.faces4.c4;c_xy4(face_in_in)];
G.faces4.x_dyz=[G.faces4.x_dyz;(G.cells.centroids(c_xy3(face_in_in),3)-G.cells.centroids(c_xy1(face_in_in),3))./2];
G.faces4.describe=[G.faces4.describe;ones(size(c_xy4(face_in_in)))*2];
G.faces4.n=[G.faces4.n;face_in_in_x'];
%%
index_x=(nr+1)*nz*ny+1: (nr+1)*nz*ny+(nr)*nz*(ny+1);
face_in=(G.faces.neighbors(:,1).*G.faces.neighbors(:,2))~=0;
face_in_x=index_x(face_in(index_x));

[neighborship, ~] = getNeighbourship(G, 'Topological', true);
face_cell_x=neighborship(face_in_x,:);
[cellNo, cf, ~] = getCellNoFaces(G);%cellNo records the location of faces for each cell. cf records the face marks of each cell.
f_xy1=cf((face_cell_x(:,1)-1).*6+1);
f_xy2=cf((face_cell_x(:,2)-1).*6+1);
f_xy3=cf((face_cell_x(:,1)-1).*6+2);
f_xy4=cf((face_cell_x(:,2)-1).*6+2);
c_xy1=neighborship(f_xy1,1);
c_xy2=neighborship(f_xy2,1);
c_xy3=neighborship(f_xy3,2);
c_xy4=neighborship(f_xy4,2);
face_in_in=c_xy1.*c_xy2.*c_xy3.*c_xy4~=0;
face_in_in_x=face_in_x(face_in_in);
G.faces4.c1=[G.faces4.c1;c_xy1(face_in_in)];
G.faces4.c2=[G.faces4.c2;c_xy2(face_in_in)];
G.faces4.c3=[G.faces4.c3;c_xy3(face_in_in)];
G.faces4.c4=[G.faces4.c4;c_xy4(face_in_in)];
G.faces4.x_dyz=[G.faces4.x_dyz;(G.cells.centroids(c_xy3(face_in_in),1)-G.cells.centroids(c_xy1(face_in_in),1))./2];
G.faces4.describe=[G.faces4.describe;ones(size(c_xy4(face_in_in)))*3];
G.faces4.n=[G.faces4.n;face_in_in_x'];
%%

%%
index_x=(nr+1)*nz*ny+1: (nr+1)*nz*ny+(nr)*nz*(ny+1);
face_in=(G.faces.neighbors(:,1).*G.faces.neighbors(:,2))~=0;
face_in_x=index_x(face_in(index_x));

[neighborship, ~] = getNeighbourship(G, 'Topological', true);
face_cell_x=neighborship(face_in_x,:);
[cellNo, cf, ~] = getCellNoFaces(G);%cellNo records the location of faces for each cell. cf records the face marks of each cell.
f_xy1=cf((face_cell_x(:,1)-1).*6+5);
f_xy2=cf((face_cell_x(:,2)-1).*6+5);
f_xy3=cf((face_cell_x(:,1)-1).*6+6);
f_xy4=cf((face_cell_x(:,2)-1).*6+6);
c_xy1=neighborship(f_xy1,1);
c_xy2=neighborship(f_xy2,1);
c_xy3=neighborship(f_xy3,2);
c_xy4=neighborship(f_xy4,2);
face_in_in=c_xy1.*c_xy2.*c_xy3.*c_xy4~=0;
face_in_in_x=face_in_x(face_in_in);
G.faces4.c1=[G.faces4.c1;c_xy1(face_in_in)];
G.faces4.c2=[G.faces4.c2;c_xy2(face_in_in)];
G.faces4.c3=[G.faces4.c3;c_xy3(face_in_in)];
G.faces4.c4=[G.faces4.c4;c_xy4(face_in_in)];
G.faces4.x_dyz=[G.faces4.x_dyz;(G.cells.centroids(c_xy3(face_in_in),3)-G.cells.centroids(c_xy1(face_in_in),3))./2];
G.faces4.describe=[G.faces4.describe;ones(size(c_xy4(face_in_in)))*4];
G.faces4.n=[G.faces4.n;face_in_in_x'];
%%
%%
index_x= (nr+1)*nz*ny+(nr)*nz*(ny+1)+1:(nr+1)*nz*ny+(nr)*nz*(ny+1)+(nr)*(nz+1)*(ny);
face_in=(G.faces.neighbors(:,1).*G.faces.neighbors(:,2))~=0;
face_in_x=index_x(face_in(index_x));

[neighborship, ~] = getNeighbourship(G, 'Topological', true);
face_cell_x=neighborship(face_in_x,:);
[cellNo, cf, ~] = getCellNoFaces(G);%cellNo records the location of faces for each cell. cf records the face marks of each cell.
f_xy1=cf((face_cell_x(:,1)-1).*6+1);
f_xy2=cf((face_cell_x(:,2)-1).*6+1);
f_xy3=cf((face_cell_x(:,1)-1).*6+2);
f_xy4=cf((face_cell_x(:,2)-1).*6+2);
c_xy1=neighborship(f_xy1,1);
c_xy2=neighborship(f_xy2,1);
c_xy3=neighborship(f_xy3,2);
c_xy4=neighborship(f_xy4,2);
face_in_in=c_xy1.*c_xy2.*c_xy3.*c_xy4~=0;
face_in_in_x=face_in_x(face_in_in);
G.faces4.c1=[G.faces4.c1;c_xy1(face_in_in)];
G.faces4.c2=[G.faces4.c2;c_xy2(face_in_in)];
G.faces4.c3=[G.faces4.c3;c_xy3(face_in_in)];
G.faces4.c4=[G.faces4.c4;c_xy4(face_in_in)];
G.faces4.x_dyz=[G.faces4.x_dyz;(G.cells.centroids(c_xy3(face_in_in),1)-G.cells.centroids(c_xy1(face_in_in),1))./2];
G.faces4.describe=[G.faces4.describe;ones(size(c_xy4(face_in_in)))*5];
G.faces4.n=[G.faces4.n;face_in_in_x'];
%%



%%
index_x= (nr+1)*nz*ny+(nr)*nz*(ny+1)+1:(nr+1)*nz*ny+(nr)*nz*(ny+1)+(nr)*(nz+1)*(ny);
face_in=(G.faces.neighbors(:,1).*G.faces.neighbors(:,2))~=0;
face_in_x=index_x(face_in(index_x));

[neighborship, ~] = getNeighbourship(G, 'Topological', true);
face_cell_x=neighborship(face_in_x,:);
[cellNo, cf, ~] = getCellNoFaces(G);%cellNo records the location of faces for each cell. cf records the face marks of each cell.
f_xy1=cf((face_cell_x(:,1)-1).*6+3);
f_xy2=cf((face_cell_x(:,2)-1).*6+3);
f_xy3=cf((face_cell_x(:,1)-1).*6+4);
f_xy4=cf((face_cell_x(:,2)-1).*6+4);
c_xy1=neighborship(f_xy1,1);
c_xy2=neighborship(f_xy2,1);
c_xy3=neighborship(f_xy3,2);
c_xy4=neighborship(f_xy4,2);
face_in_in=c_xy1.*c_xy2.*c_xy3.*c_xy4~=0;
face_in_in_x=face_in_x(face_in_in);
G.faces4.c1=[G.faces4.c1;c_xy1(face_in_in)];
G.faces4.c2=[G.faces4.c2;c_xy2(face_in_in)];
G.faces4.c3=[G.faces4.c3;c_xy3(face_in_in)];
G.faces4.c4=[G.faces4.c4;c_xy4(face_in_in)];
G.faces4.x_dyz=[G.faces4.x_dyz;(G.cells.centroids(c_xy3(face_in_in),2)-G.cells.centroids(c_xy1(face_in_in),2))./2];
G.faces4.describe=[G.faces4.describe;ones(size(c_xy4(face_in_in)))*6];
G.faces4.n=[G.faces4.n;face_in_in_x'];

end