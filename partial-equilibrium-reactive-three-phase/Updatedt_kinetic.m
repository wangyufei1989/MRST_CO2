function dt=Updatedt_kinetic(state,dt)






    dpm=100000*5;
dsm=1e-3;

state.ds=abs(state.ds);
state.dp=abs(state.dp);
a=1;

dt=dt*min([5,0.5*abs(dpm/max(max(state.dp(:,1))))*a,0.5*abs(dsm/max(state.ds(:,1)))]);




dt=min([dt, state.CFL*0.1,5e3]);

end

