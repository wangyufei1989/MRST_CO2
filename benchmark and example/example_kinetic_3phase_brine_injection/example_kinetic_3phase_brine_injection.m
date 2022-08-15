%% install the mrst_co2
 clear;
%add the path of the folder of mrst_co2
addpath('../../')
run startup.m
%% define the grid system

% grid discretization; the numbers of grids in the (x,y,z) direction are (nx,ny,nz);
nx=40;
ny=1;
nz=20;
dims= [nx ny nz];
% the domain sizes in the (x,y,z) directions are (distance, thickness,
% depth);
distance=0.02;
thickness=0.005;
depth=0.01;
domain=[distance thickness  depth]; 
% generate cartesian grid  system; the dimension can be three;
Grid= computeGeometry(cartGrid(dims,domain));
% load the data for natural log permeability that has zero mean and unit variance; 
load('perm_data.dat');
% transfer the natural log permeability to permeability with geometric mean of 1e-11;
Rock.perm=exp(log(1e-11)+sqrt(1).*perm_data(:,2));
% assign uniform porosity to all the grid;
Rock.poro = 0.2*ones( nx*ny*nz,1);
%% define the fluid property; the retension curve and relative permeability;

%J-fucntion; Leverette function 
pJ=sqrt(mean(Rock.poro)./harmmean(Rock.perm)).*sqrt(Rock.perm./Rock.poro);
% residual saturation (s_lr,s_gr); maximum relative permeability
% (k_lm,k_gm);
s_lr=0;s_gr=0;
k_lm=1;k_gm=1;
alpha_van=5;%[bar]
m_van=0.8;
%  van Genuchten model or simple model can be used; in the simple model;
%  the relative permeability is simply proportional to the saturation.
van_genuchtten=1;
if van_genuchtten==1
    Fluid = initFluidPPVGM('sr' , [ s_lr, s_gr] , ...
'kwm', [ k_lm, k_gm], 'pc_scale', pJ,'alpham',[alpha_van,m_van]);
else
Fluid = initCoreyFluidPP_3p_benchmark('sr' , [ s_lr, s_gr] , ...
'kwm', [ 1, 1], 'pc_scale', pJ,'alpham',[alpha_van,m_van]);
end
% set the gravity; if the gravity is not necessary, then use 'gravity reset off';
gravity reset off
% set the molality of NaCl:  m_NaCl;
m_NaCl=0.5;
% set the temperature
Tk=273.15+60;
%% define the well location and rate
% the injected volume rate;
darcyv=1e-5;
rate=depth*thickness*darcyv;
% injected phase; [0 1] means inject gas; [1 0] means inject brine;
phase=[1 0];
% well radius 
r_w=min(distance/10/nx,thickness/3);
% add the vertical injection well;  the (x,y) location of the well is
% (wx,wy), which means the x location of the well is the wx-th grid; the y
% location of the well is the wy-th grid; the well pierces from the wz1
% layer to the wz2 layer in the vertical direaction; 1<=wx<=nx; 1<=wy<=ny; 1<=wz1<=wz2<=nz;
wx=1;wy=1;wz1=1;wz2=nz;
% initial well pressure; normally, this value should be slightly larger or equal to the initial pressure; 
pW=150.01*barsa;
Well = verticalWell([],Grid, Rock,wx,wy, (wz1:wz2), 'Type', 'rate', ...
    'Val', rate,'name', 'i1', 'radius', r_w, 'Comp_i', phase,'pressure',pW);
%mass fractions of four components at gas pressure 145 bar; the mineral
% species has been excluded.
Well.frac=[0.914449132425116,0.0589177386616989,0,-4.77863038865764e-05];
% if you need more than one injection well, then use
%Well = verticalWell(Well,Grid, Rock,wx1,wy1, (wz1:wz2), 'Type', 'rate', ...
   % 'Val', rate,'name', 'i1', 'radius', r_w, 'Comp_i', phase,'pressure',pW);
%% initialize the fluid system

%initial liquid pressure and gas pressure at the top (p_l0,p_g0); the
%pressures increase downward under gravity force;
p_l0=150*barsa;
p_g0=1.00*barsa;
density_ref=[1000 640 2700];
%dispersion parameters of two phases: [dispersivity_longitudinal,
%dispersivity_transverse,molecular_diffusion];
D_liquid=[0 0 1e-9];D_gas=[0 0 1e-9];
State=initStatePP_kinetic_unstructure(Grid,Well,[p_l0,p_g0],density_ref,Tk,Fluid,m_NaCl,Rock,D_liquid,D_gas); 
% if you want to use the initial data from observation, modify the 'State' at
% here; e.g., State.pressure(:,1)=[...], State.pressure(:,2)=[...],

%%  pressure boundary
% the boundnary pressure now only allow the pressure boundary.
% if bc=[], then all the boundaries are closed. here 'xmax' means the right
% most boundary; similarly,'zmax' means the bottom boundary; we give the
% pressure at the right side;
pB=150.0*barsa;
bc  = pside([], Grid, 'xmax', pB,State,Fluid); 
% if you need more than one boundary 
%bc  = pside(bc, Grid, 'xmin', 150.0*barsa(),State,Fluid); 
%bc.type='pressure';
%% simulation time and initial time step 
% assign the total simulation time 
T= distance/darcyv;
% initial time step. very small. 
dt=1e-5;

%% error monitor
% tolerance is based on the maximum allowed residual pressure; when the
% pressure change is smaller than dplim, we deem the result converges.
State.dplim=1e-3;
% we monitor the total mass of the co2 component (mc) and water component (mb);
mc=sum((State.rho(:,1).*State.s(:,1).*State.frac(:,2)).*Rock.poro.*Grid.cells.volumes)+sum((State.rho(:,2).*State.s(:,2).*State.frac(:,4)).*Rock.poro.*Grid.cells.volumes);
mb=sum((State.rho(:,1).*State.s(:,1).*State.frac(:,1)).*Rock.poro.*Grid.cells.volumes)+sum((State.rho(:,2).*State.s(:,2).*State.frac(:,3)).*Rock.poro.*Grid.cells.volumes);
co2=0;
water=0;
%record the error
ERR=[];
%record the time step;
DT=[];
%% store frequency 
dst=T/50;
%%
jc=1;
t=0;
State.solve=1;
while t<T
    
    if State.solve==1
        t=t+dt;
    end
    [State,Well,dt]  =   kr3p(State, Grid, Fluid,Rock ,dt,'wells', Well,'bc',bc);
    co2=co2+Well.val*State.rhoref(1)*dt.*Well.frac(2)-State.outflux(2);
    water=water+Well.val*State.rhoref(1)*dt.*Well.frac(1)-State.outflux(1);
    err2=(sum((State.rho(:,1).*State.s(:,1).*State.frac(:,2)).*State.poro.*Grid.cells.volumes)...
        -co2 -mc)...
        /(sum((State.rho(:,1).*State.s(:,1).*State.frac(:,2)).*State.poro.*Grid.cells.volumes));
    err1=(sum((State.rho(:,1).*State.s(:,1).*State.frac(:,1)).*State.poro.*Grid.cells.volumes)...
        -water -mb)...
        /(sum((State.rho(:,1).*State.s(:,1).*State.frac(:,1)).*State.poro.*Grid.cells.volumes));
    
    ERR=[ERR;err1;err2];
    
    if  abs(t-jc*dst)<3*dt && State.solve==1&&t<T
        dissolve=(State.poro-State.poro0)./dt;
        savei=num2str(jc);
        savex=strcat('State_kinetic',savei,'.mat');
        save(savex,'State')
        savet=strcat('t_kinetic',savei,'.mat');
        save(savet,'t')
        jc=jc+1;
    end
    
    DT=[DT;dt];
    State.meandt=mean(DT);
    dt=Updatedt_kinetic(State,dt);
end


%% figures
% porosity change
figure
X=Grid.cells.centroids(1:nx);
Y=linspace(depth,0,nz);
[X,Y]=meshgrid(X,Y);
S=reshape((State.poro(:)-0.2)./0.2,nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
xlabel('r (m)')
ylabel('z (m)')
colorbar
pbaspect([2 1 1])
colormap hot

% permeability
figure
X=Grid.cells.centroids(1:nx);
Y=linspace(depth,0,nz);
[X,Y]=meshgrid(X,Y);
S=reshape(log10(Rock.perm(:)),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
xlabel('r (m)')
ylabel('z (m)')
colorbar
pbaspect([2 1 1])
colormap hot

% liquid pressure
figure
 X=Grid.cells.centroids(1:nx);
Y=linspace(depth,0,nz);
[X,Y]=meshgrid(X,Y);
S=reshape(State.pressure(:,1),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlabel('r (m)')
ylabel('z (m)')
colorbar
pbaspect([2 1 1])
colormap hot

% gas pressure
figure
 X=Grid.cells.centroids(1:nx);
Y=linspace(depth,0,nz);
[X,Y]=meshgrid(X,Y);
S=reshape(State.pressure(:,2),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlabel('r (m)')
ylabel('z (m)')
colorbar
pbaspect([2 1 1])
colormap hot

  
% pH
figure
 X=Grid.cells.centroids(1:nx);
Y=linspace(depth,0,nz);
[X,Y]=meshgrid(X,Y);
S=reshape(-log10(State.species(:,4)./State.species(:,1)./0.001),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
xlabel('r (m)')
ylabel('z (m)')
colorbar
pbaspect([2 1 1])
colormap hot
  
 
   % mass fraction of co2 component
 figure
 X=Grid.cells.centroids(1:nx);
Y=linspace(depth,0,nz);
[X,Y]=meshgrid(X,Y);
S=reshape(State.frac(:,2),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlim([0 distance])
xlabel('r (m)')
ylabel('z (m)')
colorbar
pbaspect([2 1 1])
colormap hot

% mineral component
 figure
 X=Grid.cells.centroids(1:nx);
Y=linspace(depth,0,nz);
[X,Y]=meshgrid(X,Y);
S=reshape(State.frac(:,3),nx,nz);
contourf(X,Y,S', 'LineStyle', 'none')
xlabel('x (m)')
ylabel('z (m)')
colorbar
pbaspect([2 1 1])
colormap hot

 % velocity in x direction
 figure
X=Grid.cells.centroids(1:nx,1);
X=[0;X];
Y=linspace(depth,0,nz);
[X,Y]=meshgrid(X,Y);
S=reshape(State.flux(1:(nx+1)*nz,1),nx+1,nz);
contourf(X,Y,S')
pbaspect([2 1 1])
xlabel('r (m)')
ylabel('z (m)')
colorbar
% velocity in y direction
 figure
X=Grid.cells.centroids(1:nx);
Y=linspace(depth,0,nz+1);
[X,Y]=meshgrid(X,Y);
S=reshape(State.flux(((nx+1)*nz+nx*nz*2+1):end,1),nx,nz+1);
contourf(X,Y,S')
pbaspect([2 1 1])
xlabel('x (m)')
ylabel('z (m)')
colorbar



