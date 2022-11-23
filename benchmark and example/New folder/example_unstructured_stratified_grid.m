%% install the mrst_co2
 clear;
%add the path of the folder of mrst_co2
addpath('../../')
run startup.m
%% Horizons for internal and external geology
% We begin by building the horizons that define the top and bottom
% structure of the sector, as well as one internal erosions.
% Define areal mesh
[xmax,ymax, n]  = deal(1000*meter, 1000*meter, 30);
[x, y] = meshgrid(linspace(0,xmax,n+1), linspace(0,ymax,n+1));
[x, y] = deal(x',y');
% Basic dome structure
dome = 1-exp(sqrt((x - xmax/2).^2 + (y - ymax/2).^2)*1e-3);
% Periodic perturbation
[xn,yn] = deal(pi*x/xmax,pi*y/ymax);
perturb = sin(5*xn) + .5*sin(4*xn+6*yn) + cos(.25*xn*yn./pi^2) + cos(3*yn);
perturb = perturb/3.5;
% Random small-scale perturbation
rng(0);
[h, hr] = deal(8,1);
zt = 50 + h*perturb + rand(size(x))*hr - 20*dome;
zb = zt + 200;
zmb = min(zb -54 + 0.08*x - 0.020*y + hr*rand(size(x)), zb);
zmt = max(zb -135 -0.01*x - 0.025*y + hr*rand(size(x)), zt);
horizons = {struct('x', x, 'y', y, 'z', zt+1500), ...
    struct('x', x, 'y', y, 'z', zmt+1500), ...
    struct('x', x, 'y', y, 'z', zmb+1500), ...
    struct('x', x, 'y', y, 'z', zb+1500)};

surf(x,y,zt+1500-.2, 'EdgeC','r','FaceC',[.8 .8 .8]),  hold on
mesh(x,y,zmt+1500-.2,'EdgeC','g','FaceC',[.7 .7 .7]), 
mesh(x,y,zmb+1500-.2,'EdgeC','k','FaceC',[.6 .6 .6]), 
mesh(x,y,zb+1500-.2, 'EdgeC','b','FaceC',[.5 .5 .5]); hold off
set(gca,'ZDir','reverse')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
xticks([0 500 1000])
yticks([0 500 1000])
pbaspect([1 1 .3])
view(-50,10); %axis off
%% Interpolate to build unfaulted corner-point grid
dims = [20, 20]; layers = [1 2 1];
grdecl = convertHorizonsToGrid(horizons, 'dims', dims, 'layers', layers);
G = processGRDECL(grdecl);
[~,~,k] = gridLogicalIndices(G);
k(k==2|k==3)=3;
figure;plotCellData(G,k,'EdgeAlpha',.2); view(3);
colormap(.5*jet + .5*ones(size(jet)));
view(-50,10)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
xticks([0 500 1000])
yticks([0 500 1000])
pbaspect([1 1 .3])
%axis off
%% Insert faults
[X,Y,Z]  = buildCornerPtNodes(grdecl);
i=27:40; Z(i,:,:) = Z(i,:,:) + .022*min(0,Y(i,:,:)-550);
j= 1:30; Z(:,j,:) = Z(:,j,:) + .021*min(0,X(:,j,:)-400);
%j=27:40; Z(:,j,:) = Z(:,j,:) + .023*min(0,X(:,j,:)-750);
grdecl.ZCORN = Z(:);
G = processGRDECL(grdecl);
G = computeGeometry(G);
[~,~,k] = gridLogicalIndices(G);
k(k==2|k==3)=3;
figure, plotCellData(G,k,'EdgeAlpha',.2); view(3);
plotFaces(G,find(G.faces.tag>0),'EdgeColor','r','FaceColor',[.8 .8 .8]);
colormap(.5*jet + .5*ones(size(jet)));
view(-50,10)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
xticks([0 500 1000])
yticks([0 500 1000])
pbaspect([1 1 .3])
%% remove very small grids
[G,gs] = extractSubgrid(G, G.cells.volumes>1e-1*max(G.cells.volumes));
figure, plotCellData(G,k(gs),'EdgeAlpha',.2); view(3);
colormap(.5*jet + .5*ones(size(jet)));
view(-50,10)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
xticks([0 500 1000])
yticks([0 500 1000])
pbaspect([1 1 .3])

%% Petrophysics
% Set up permeability based on K-indices and introduce anisotropy by
% setting K_z = .1*K_x
rng(357371);
[K,L] = logNormLayers(G.cartDims, [10 400  50]*milli*darcy);
K = K(G.cells.indexMap);
perm = [K, K, 0.1*K];
rock = makeRock(G, perm, 0.3);
% Plot horizontal permeability
figure
%K = convertTo(K,milli*darcy);
plotCellData(G, log10(K),'EdgeAlpha',.1)
view(-50,10)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
xticks([0 500 1000])
yticks([0 500 1000])
pbaspect([1 1 .3])
colorbar
%% define the fluid property; the retension curve and relative permeability;
%J-fucntion; Leverette function 
pJ=sqrt(mean(rock.poro)./harmmean(max(rock.perm,[],2))).*sqrt(max(rock.perm,[],2)./rock.poro);
% residual saturation (s_lr,s_gr); maximum relative permeability
% (k_lm,k_gm);
s_lr=0;s_gr=0;
k_lm=1;k_gm=1;
alpha_van=5;%[bar]
m_van=0.8;
%  van Genuchtten model or simple model can be used; in the simple model;
%  the relative permeability is simply proportional to the saturation.
van_genuchtten=1;
if van_genuchtten==0
    fluid = initFluidPPVGM('sr' , [ s_lr, s_gr] , ...
'kwm', [ k_lm, k_gm], 'pc_scale', pJ,'alpham',[alpha_van,m_van]);
else
fluid = initCoreyFluidPP_3p_benchmark('sr' , [ s_lr, s_gr] , ...
'kwm', [ 1, 1], 'pc_scale', pJ,'alpham',[alpha_van,m_van]);
end
% set the gravity; if the gravity is not necessary, then use 'gravity reset off';
gravity reset on
% set the molality of NaCl:  m_NaCl;
m_NaCl=0.5;
% set the temperature
Tk=273.15+60;
%% define the well location and rate
% the injected volume rate;
pv      = poreVolume(G, rock);
rate= 1*sum(pv)/(10*365*24*3600);
% injected phase; [0 1] means inject gas; [1 0] means inject brine;
phase=[0 1];
% well radius 
r_w=0.1;
% add the vertical injection well; w_cell denotes the index of pierced cells.
pW=50.001*barsa;
w_cell=find((abs(G.cells.centroids(:,1)-24.9)<2)&(abs(G.cells.centroids(:,2)-24.9)<2)&((G.cells.centroids(:,3))>1633)&((G.cells.centroids(:,3))<1714));
W = addWell([],G, rock,w_cell, 'Type', 'rate', ...
    'Val', rate,'name', 'i1', 'radius', r_w, 'Comp_i', phase,'pressure',pW);
w_cell=find((abs(G.cells.centroids(:,1)-975)<2)&(abs(G.cells.centroids(:,2)-25)<2)&((G.cells.centroids(:,3))>1620)&((G.cells.centroids(:,3))<1762));
W = addWell(W,G, rock,w_cell, 'Type', 'rate', ...
    'Val', rate,'name', 'i2', 'radius', r_w, 'Comp_i', phase,'pressure',pW);
% Plot the wells
figure
plotWell(G, W,'color','k')
[~,~,k] = gridLogicalIndices(G);
k(k==2|k==3)=3;
k([W(1).cells;W(2).cells])=2;
plotCellData(G, k,'EdgeAlpha',.2)
colormap(.5*jet + .5*ones(size(jet)));
view(-50,10)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
xticks([0 500 1000])
yticks([0 500 1000])
pbaspect([1 1 .3])
close all
hT = computeTrans(G, rock);% calculate the half transmisivity of each face
%% initialize the fluid system
%initial liquid pressure and gas pressure at the top (p_l0,p_g0); the
%pressures increase downard under gravity force;
p_l0=150*barsa;
p_g0=1.00*barsa;
density_ref=[1000 640 2700];
viscosity_ref=[1e-4,1e-5];
property=5;% 
structuregrid=0;
if structuregrid==1
x=initStatePP_equilibrium(G,W,[p_l0,p_g0],density_ref,viscosity_ref,Tk,fluid,m_NaCl,rock,property); 
else
 x=initStatePP_equilibrium1(G,W,[p_l0,p_g0],density_ref,viscosity_ref,Tk,fluid,m_NaCl,rock,property); 
end
% the longitudinal dispersivity (x.al), tranverse dispersivity (x.at) and
% diffusion coefficient (x.dm); first brine and last gas; 
x.al=[5 5]; x.at=[1 1]; x.dm=[1e-9 1e-9];
% if you want to use the initial data from observation, modify the 'x' at
% here; e.g., x.pressure(:,1)=[...], x.pressure(:,2)=[...],
%%  pressure boundary
% the boudnary pressure now only allow the pressure boudary.
% if bc=[], then all the boundaries are closed. here 'xmax' means the right
% most boundary; similarly,'zmin' means the bottom boudnary; we give the
% pressure at the top layer;
pB=150.0*barsa;
bc  = pside([], G, 'ymax', pB,x,fluid); 
bc.type='pressure';
% if you need more than one boundary 
%bc  = pside(bc, G, 'xmin', 150.0*barsa(),x,fluid); 
%bc.type='pressure';
%% simulation time and initial time step 
% assign the total simulation time 
T= 365*24*3600*1;
% initial time step. very small. 
dt=1e-1;
%% error monitor
% tolerance is based on the maximum allowed residual pressure; when the
% pressure change is smaller than dplim, we deem the result converges.
x.dplim=1e-3;
% we monitor the total mass of the co2 component (mc) and water component (mb);
mc=sum((x.rho(:,1).*x.s(:,1).*x.frac(:,2)).*rock.poro.*G.cells.volumes)+sum((x.rho(:,2).*x.s(:,2).*x.frac(:,4)).*rock.poro.*G.cells.volumes);
mb=sum((x.rho(:,1).*x.s(:,1).*x.frac(:,1)).*rock.poro.*G.cells.volumes)+sum((x.rho(:,2).*x.s(:,2).*x.frac(:,3)).*rock.poro.*G.cells.volumes);
co2=0;
water=0;
%record the error
ERR=[];
%record the time step;
DT=[];
%% store frequency 
dst=T/50;
%% main loop
jc=1;
t=0;
x.solve=1;
while t<T
    if x.solve==1
        t=t+dt;
    end

        [x,W,dt]  =   kr3p(x, G, fluid,rock ,dt,'wells', W,'bc',bc);
    if  abs(t-jc*dst)<3*dt && x.solve==1&&t<T
        dissolve=(x.poro-x.poro0)./dt;
        savei=num2str(jc);
        savex=strcat('x_stra_k',savei,'.mat');
        save(savex,'x')
        savet=strcat('t_stra_k',savei,'.mat');
        save(savet,'t')
        jc=jc+1;
    end
    DT=[DT;dt];
    x.meandt=mean(DT);
    dt=Updatedt_3pn(x,dt,DT)
    max( x.s(:,2))
end


%% figures
save('G.mat','G')
plotCellData(G, x.frac(:,2),'EdgeAlpha',.1)
view(-50,10)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
xticks([0 500 1000])
yticks([0 500 1000])
pbaspect([1 1 .3])

