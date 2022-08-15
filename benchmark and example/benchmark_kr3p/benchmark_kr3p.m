%% install the mrst_co2
 clear;
%add the path of the folder of mrst_co2
addpath('..\..\')
run startup.m
%% define the grid system

% grid discretization; the numbers of grids in the (x,y,z) direction are (nx,ny,nz);
nx=500;ny=1;nz=1;
dims= [nx ny nz];
% the domain sizes in the (x,y,z) directions are (distance, thickness, depth);
distance=10;thickness=1;depth=1;
domain=[distance thickness  depth]; 
% generate cartesian grid  system; the dimension can be three.
G= computeGeometry(cartGrid(dims,domain));
% transfer the cartesian grid system to radial grid system; currently, we
% only have vertically one- or two-dimensional radial choice, i.e., ny=1;
%G=orth2radial(G);

% assign the permeability and porosity 
rock.perm=ones(nx*ny*nz,1).*1e-11;
rock.poro = 0.2*ones( nx*ny*nz,1);
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
 fluid = initFluidPPVG('sr' , [ s_lr, s_gr,s_lr,s_gr] , ...
 'kwm', [ k_lm, k_gm], 'pc_scale', pJ,'alpham',[alpha_van,m_van]);
end
% set the gravity; if the gravity is not necessary, then use 'gravity reset off';
gravity reset on
% set the molality of NaCl:  m_NaCl;
m_NaCl=0.5;
% set the temperature
Tk=273.15+60;
%% define the well location and rate
% the injected volume rate;
rate=1.3808e-4*733/640;
% injected phase; [0 1] means inject gas; [1 0] means inject brine;
phase=[0 1];
% well radius 
r_w=0.01; 
% add the vertical injection well;  the (x,y) location of the well is
% (wx,wy), which means the x location of the well is the wx-th grid; the y
% location of the well is the wy-th grid; the well pierces from the wz1
% grid to the wz2 grid in the vertical direaction; 1<=wx<=nx; 1<=wy<=ny; 1<=wz1<=wz2<=nz;
wx=1;wy=1;wz1=1;wz2=nz;
% initial gas of well pressure; normally, this value should be slightly larger or equal to the initial pressure; 
pW=50.001*barsa;
W = verticalWell([],G, rock,wx,wy, (wz1:wz2), 'Type', 'rate', ...
    'Val', rate,'name', 'i1', 'radius', r_w, 'Comp_i', phase,'pressure',pW);

% if you need more than one injection well, then use
%W = verticalWell(W,G, rock,wx1,wy1, (wz1:wz2), 'Type', 'rate', ...
   % 'Val', rate,'name', 'i1', 'radius', r_w, 'Comp_i', phase,'pressure',pW);
%% initialize the fluid system

%initial liquid pressure and gas pressure at the top (p_l0,p_g0); the
%pressures increase downard under gravity force;
p_l0=150*barsa;
p_g0=1.00*barsa;
%dispersion parameter of liquid phase: [dispersivity_longitudinal,
%dispersivity_transverse,molecular_diffusion];
D_liquid=[0 0 1e-9];D_gas=[0 0 1e-9];
structuregrid=0;
if structuregrid==1
x=initStatePP_kinetic(G,W,[p_l0,p_g0],[1000,640,2700],Tk,fluid,0.5,rock,D_liquid,D_gas); 
else
    x=initStatePP_kinetic_unstructure(G,W,[p_l0,p_g0],[1000,640,2700],Tk,fluid,0.5,rock,D_liquid,D_gas); 
end

% The temporal porosiry
x.poro=rock.poro;
x.poro0=rock.poro;


% if you want to use the initial data from observation, modify the 'x' at
% here; e.g., x.pressure(:,1)=[...], x.pressure(:,2)=[...],

%%  pressure boundary
% the boudnary pressure now only allow the pressure boudary.
% if bc=[], then all the boundaries are closed. here 'xmax' means the right
% most boundary; similarly,'zmax' means the bottom boundary; we give the
% pressure at the right side;
pB=150.0*barsa;
bc  = pside([], G, 'xmax', pB,x,fluid); 

% if you need more than one boundary 
%bc  = pside(bc, G, 'xmin', 150.0*barsa(),x,fluid); 
%bc.type='pressure';
%% simulation time and initial time step 
% assign the total simulation time 
T= 9839;
% initial time step. very small. 
dt=1e-6;

%% error monitor
% tolerance is based on the maximum allowed residual pressure; when the
% pressure change is smaller than dplim, we deem the result converges.
x.dplim=1e-3;
% we monitor the total mass of the co2 component (mc) and water component (mb);
mc=sum((x.rho(:,1).*x.s(:,1).*x.frac(:,2)).*rock.poro.*G.cells.volumes)+sum((x.rho(:,2).*x.s(:,2)).*rock.poro.*G.cells.volumes);
mb=sum((x.rho(:,1).*x.s(:,1).*x.frac(:,1)).*rock.poro.*G.cells.volumes);
co2=0;
water=0;
%record the error
ERR=[];
%record the time step;
DT=[];
%% store frequency 
dst=T/100;
%%
jc=1;
t=0;

x.solve=1;

while t<T
    if x.solve==1
        t=t+dt;
    end
    W.val=rate*t^(-0.5);
    
    
    [x,W,dt]  = kr3p(x, G, fluid,rock,dt,'wells', W,'bc',bc);
    
    co2=co2+W.val*640*dt.*W.frac(4)-x.outflux(2);
    water=water+W.val*640*dt.*W.frac(3)-x.outflux(1);
    err2=(sum((x.rho(:,1).*x.s(:,1).*x.frac(:,2)).*x.poro.*G.cells.volumes)+sum((x.rho(:,2).*x.s(:,2)).*x.poro.*G.cells.volumes)...
        -co2 -mc)...
        /(sum((x.rho(:,1).*x.s(:,1).*x.frac(:,2)).*x.poro.*G.cells.volumes)+sum((x.rho(:,2).*x.s(:,2)).*x.poro.*G.cells.volumes));
    err1=(sum((x.rho(:,1).*x.s(:,1).*x.frac(:,1)).*x.poro.*G.cells.volumes)...
        -water -mb)...
        /(sum((x.rho(:,1).*x.s(:,1).*x.frac(:,1)).*x.poro.*G.cells.volumes))
    
    ERR=[ERR;err1;err2];
    
    
    
    if  abs(t-jc*dst)<3*dt && x.solve==1&&t<T
        dissolve=(x.poro-x.poro0)./dt;
        savei=num2str(jc);
        saved=strcat('dissolveck',savei,'.mat');
        save(saved,'dissolve')
        savex=strcat('xck',savei,'.mat');
        save(savex,'x')
        savet=strcat('tck',savei,'.mat');
        save(savet,'t')
        savet=strcat('dtck',savei,'.mat');
        save(savet,'dt')
        jc=jc+1;
    end
    
    DT=[DT;dt];
    x.meandt=mean(DT);
    dt=Updatedt_kinetic(x,dt);
    max( x.s(:,2))
    
    
end
save('x4d4896.mat','x')
save('W.mat','W')



