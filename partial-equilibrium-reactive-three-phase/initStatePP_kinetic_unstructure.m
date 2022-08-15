function state = initStatePP_kinetic_unstructure(G, W, p0,rhoref,Tk,fluid,m_NaCl,rock,D_liquid,D_gas)
% initialize the system.
%
% SYNOPSIS:
%  state = initStatePP_kinetic1(G, W, p0,rhoref,muref,Tk,fluid,m_NaCl,rock,property)
%
% DESCRIPTION:
%   This function serves to initialize the system.


%
% REQUIRED PARAMETERS:
%   G  - grid system
%  W  -well information
%  p0  - initial pressure at the top of the aquifer
%   s0  - initial saturation at the top of the aquifer
% frac0  - intial mass fraction at the top of the aquifer
%rhoref, muref  - reference density and viscosity; not used if the density
%and viscosity change with pressure.
%fluid  - relative permeability and retention curves;
% m_NaCl,Tk  -  (constant) NaCl concentration and (constant) temperature 
%rock  - an object including the initial porosity and intrinsic permeability 
%property  - indicates what phases and reactions do we have.


% RETURNS:
%   state   - an object containing the initialized variables. 
%

%{

This file is part of mrst_co2 based on MRST.

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

   if numel(p0) == 2
       state.pressure = repmat(p0, [G.cells.num, 1]);
   end
    if numel(rhoref) == 3
      state.rho = repmat(rhoref, [G.cells.num, 1]);
      state.rhoref=rhoref;
   end
state.al=[D_liquid(1) D_gas(1)];
state.at=[D_liquid(2) D_gas(2)]; state.dm=[D_liquid(3) D_gas(3)];
s0=fluid.S(state);
nu=[0,0,-1,-1,1,1,0,0,0;0,1,0,0,0,0,-1,0,0;1,0,0,-1,0,-1,1,0,0;0,0,0,1,0,-1,0,1,0;-1,0,0,1,0,0,0,0,1];
omega=[1,0,0,0,-1,1,0,1,1;0,1,0,0,-1,1,1,1,0;0,0,1,0,1,0,0,0,0;0,0,0,1,2,-1,0,-2,-1];
omega_t=[1,0,0,0,-0.449498477967962,0.295245587295344,0,0.300204969254612,1.05920743179680;0,1,0,0,-1.09810868805829,0.721274398938001,1,0.733389991501275,0;0,0,1,0,2.49730275961874,0,0,0,0;0,0,0,1,0.0502919307350666,-0.0165167084582985,0,-0.0335882950890700,-0.0592544684854186];
state.nu=nu;
state.omega=omega;
state.omega_t=omega_t;
species=zeros(G.cells.num,9);
component=zeros(G.cells.num,4);
frac=zeros(G.cells.num,4);
salinity=zeros(G.cells.num,1);
species(:,1)=G.cells.volumes.*rock.poro.*rhoref(1);
species(:,2)=G.cells.volumes.*rock.poro.*rhoref(2);%we assume the CO2 has a large quantity
species(:,3)=G.cells.volumes.*(1-rock.poro).*rhoref(3);
state.Tk=Tk;
        state.m_NaCl=m_NaCl;
        state.frac=frac;
        state.frac0=frac;
        state.salinity=salinity;
         state.salinity0=salinity;
        state.species=species;
         state.species0=species;
        state.component=component;
        state.component0=component;
        state.poro=rock.poro;
       state.poro0=rock.poro;
state.s=s0;
state.s0=s0;
state=species_equilibrium(state,G,rock);
property=5;
 if (p0(:,2)>p0(:,1)) && property==5 % two phases; compressible soluble& reaction with rock
    
      [rho0,mu0]=rhomu_p_frac_kinetic(state);
       s0=fluid.S(state);
       rhoeff0=effective_rho(rho0,mu0,s0,fluid);
       g=gravity;
       state.pressure(:,1)=state.pressure(:,1)+rhoeff0(:,1).*(G.cells.centroids(:,3)-min(G.cells.centroids(:,3))).*g(3);
       state.pressure(:,2)=state.pressure(:,2)+rhoeff0(:,1).*(G.cells.centroids(:,3)-min(G.cells.centroids(:,3))).*g(3);   
     
  
      [rho0,mu0]=rhomu_p_frac_kinetic(state);
       [state.rho,state.mu]=rhomu_p_frac_kinetic(state);
        state=species_equilibrium(state,G,rock);
          s0=fluid.S(state);
       rhoeff0=effective_rho(rho0,mu0,s0,fluid);
         
  

      
      state=species_equilibrium(state,G,rock);

 
        state.s0=fluid.S(state);
        state.s=state.s0;
      [state.rho0,state.mu0]=rhomu_p_frac_kinetic(state);
      state.rho=state.rho0;
      state.mu=state.mu0;
      state.species0=state.species;
      state.component0=state.component;
      state.frac0=state.frac;
      state.salinity0=state.salinity;
  end
       
       
       
        if (p0(:,2)<=p0(:,1))  && property==5% one phase
  
      [rho0,mu0]=rhomu_p_frac_kinetic(state);
       s0=fluid.S(state);
       rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        g=gravity;
       state.pressure(:,1)=state.pressure(:,1)+rhoeff0(:,1).*(G.cells.centroids(:,3)-min(G.cells.centroids(:,3))).*g(3);
           state.pressure(:,2)=state.pressure(:,2)+rhoeff0(:,1).*(G.cells.centroids(:,3)-min(G.cells.centroids(:,3))).*g(3);   
     
          [rho0,mu0]=rhomu_p_frac_kinetic(state);
          [state.rho,state.mu]=rhomu_p_frac_kinetic(state);
               state=species_equilibrium(state,G,rock);
           s0=fluid.S(state);
       rhoeff0=effective_rho(rho0,mu0,s0,fluid);
       

for refine=1:50
         state=species_equilibrium(state,G,rock);
end

 
        state.s0=fluid.S(state);
        state.s=state.s0;
      [state.rho0,state.mu0]=rhomu_p_frac_kinetic(state);
      state.rho=state.rho0;
      state.mu=state.mu0;
      state.species0=state.species;
      state.component0=state.component;
      state.frac0=state.frac;
      state.salinity0=state.salinity;
    
  end
 


%state=species_equilibrium(state,G,rock);

   if ~isempty(W)
      state.wellSol = initWellSol(W, p0);
   end
end

function rhoeff = effective_rho(rho0,mu0,s0, fluid)
   
 
   kr        = fluid.relperm(s0);

   mob    = bsxfun(@rdivide, kr, mu0);
   totmob = sum(mob, 2);
   rhoeff = sum(bsxfun(@times, mob, rho0(:,1:2)), 2) ./ totmob;
end

function state=species_equilibrium(state,G,rock)


for i=1:G.cells.num
   % component_t=state.omega_t*state.species(i,:)';
species_t=reaction_equilibrium(state.pressure(i,1),state.pressure(i,2),state.nu,state.omega,state.species(i,:));
state.species(i,:)=species_t';
state.species(i,2)=state.rho(i,2).*state.s(i,2).*rock.poro(i).*G.cells.volumes(i);
state.species(i,3)=state.rho(i,3).*(1-rock.poro(i)).*G.cells.volumes(i);
state.species(i,1)=state.rho(i,1).*state.s(i,1).*rock.poro(i).*G.cells.volumes(i).*state.species(i,1)./(sum(state.species(i,4:end),2)+state.species(i,1)*state.m_NaCl*0.05844+state.species(i,1));
state.species(i,4:end)=state.rho(i,1).*state.s(i,1).*rock.poro(i).*G.cells.volumes(i).*state.species(i,4:end)./(sum(state.species(i,4:end),2)+state.species(i,1)*state.m_NaCl*0.05844+state.species(i,1));

state.component(i,:)=state.omega_t*species_t;
speciesw=[state.species(i,1);0;0;0];
state.frac(i,:)=(state.omega_t(:,4:end)*state.species(i,4:end)'+speciesw)'./(sum(state.species(i,:))-state.species(i,2)-state.species(i,3)+state.species(i,1)*state.m_NaCl*0.05844);
state.salinity(i)=(sum(state.species(i,4:end),2)+state.species(i,1)*state.m_NaCl*0.05844)./(sum(state.species(i,4:end),2)+state.species(i,1)*state.m_NaCl*0.05844+state.species(i,1));
end

     

end


