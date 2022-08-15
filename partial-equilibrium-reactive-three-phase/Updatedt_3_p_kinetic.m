function dt=Updatedt_3_p_kinetic(state,dt,DT)






    dpm=100000*50;
dsm=1e-2;

state.ds=abs(state.ds);
state.dp=abs(state.dp);

a=1;

dt=dt*min([1.5,0.5*abs(dpm/max(max(state.dp(:,1)),max(state.dp(:,1))))*a]);
%dt=dt*min([1.5,0.5*abs(dpm/max(max(state.dp(:,1),state.dp(:,2))))*a,0.5*abs(dsm/max(state.ds(:,1)))]);


dt=min([dt, state.CFL*4e-1,2e5]);

 if size(DT,1)>100 &&all(DT(end-20:end)<5) 
     dt=rand(1,1)*state.meandt*2;
 end
%for compressible insoluble 0.00001;
%for incompressible 0.001
%for compressible soluble0.01

%dt=dt*min([1,dpm/max(state.dp)]);


end

