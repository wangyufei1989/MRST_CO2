function dt=Updatedt_comp_soluble_diffusion(state,dt,DT,dzmin,D)
% Update dt

dpm=100000*20;
dsm=1e-3;
state.ds=abs(state.ds);
state.dp=abs(state.dp);
a=1;
dt=dt*min([1.5,0.5*abs(dpm/max(max(state.dp(:,1)),max(state.dp(:,1))))*a]);
dt=min([dt, state.CFL*2e-1,dzmin^2/2/(D).*0.2,2e7]);
 if size(DT,1)>100 &&all(DT(end-20:end)<5) 
     dt=rand(1,1)*state.meandt*2;
 end
%for compressible insoluble 0.00001;
%for incompressible 0.001
%for compressible soluble0.01
%dt=dt*min([1,dpm/max(state.dp)]);
end

