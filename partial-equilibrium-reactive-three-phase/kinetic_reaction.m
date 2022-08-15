function state=kinetic_reaction(state,V,DT)

% speciation calculation based on the kinetic reaction given in section called
%Partial-Equilibrium Reactive Multi-Component Three-Phase Flow Mode in the user'guides
%
% SYNOPSIS:
% state=kinetic_reaction(state,V,DT)
%
% DESCRIPTION:
%   This function serves to calculate the mass fraction of each species and components. 
%Description is already included in the user's guide. 

%
% REQUIRED PARAMETERS:
% state   -state varibles
% V         - cell volume
% DT      -time step 

% RETURNS:
%   state   - contains the mass fractions of all species and components at equilibrium state; 
%                and also their derivatives over gas pressure. 

% SEE ALSO:
%  kr3p

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

%% main file for kinetic reaction

state.species(:,3)=(1-state.poro0).*state.rhoref(3).*V;
species=(state.species);
component=state.component;
frac=state.frac;
salinity=(state.salinity);
state.dsp=zeros(size(V));
state.dcp=zeros(size(V));
for i=1:size(state.pressure,1)


      Sr=V(i).*state.poro0(i).*4.69e4/1e0;
   %p1=ones(size(state.pressure,1),1).*1e7;
    [ nt,comp]=kinetic_reaction_sub(state.pressure(i,1),state.pressure(i,2),state.nu,state.omega,state.omega_t,state.species(i,:),state.component(i,:),state.Tk,Sr,DT);
     state.species(i,:)=nt;
     state.component(i,:)=comp';
   
  
%   
% if state.pressure(i,2)<state.pressure(i,1)
%           state.species(i,7)=state.species(i,7)+state.species(i,2);
%          state.species(i,2)=0;
%  
% end
%   

    
      %dpp=0.1.*state.pressure(i,2);%if we increase the pressure, the CO2 would not be enough if there is not gaseous CO2
      dpp=-1e0;
      dp=state.pressure(i,2)+dpp;
component(i,2)=component(i,2)*10;
     [ nt,comp]=kinetic_reaction_sub(state.pressure(i,1),dp,state.nu,state.omega,state.omega_t,species(i,:),component(i,:),state.Tk,Sr,DT);
   
     species(i,:)=nt;

   component(i,:)=comp';
%         if state.pressure(i,2)<state.pressure(i,1)
%           species(i,7)=species(i,7)+state.species(i,2);
%          species(i,2)=0;
%        end
   
   
 speciesw=[state.species(i,1);0;0;0];
state.frac(i,:)=(state.omega_t(:,4:end)*state.species(i,4:end)'+speciesw)'./(sum(state.species(i,:))-state.species(i,2)-state.species(i,3)+state.species(i,1)*state.m_NaCl*0.05844);
 specieswt=[species(i,1);0;0;0];
frac(i,:)=(state.omega_t(:,4:end)*species(i,4:end)'+specieswt)'./(sum(species(i,:))-species(i,2)-species(i,3)+species(i,1)*state.m_NaCl*0.05844);

state.salinity(i)=(sum(state.species(i,4:end),2)+state.species(i,1)*state.m_NaCl*0.05844)./(sum(state.species(i,4:end),2)+state.species(i,1)*state.m_NaCl*0.05844+state.species(i,1));
salinity(i)=(sum(species(i,4:end),2)+species(i,1)*state.m_NaCl*0.05844)./(sum(species(i,4:end),2)+species(i,1)*state.m_NaCl*0.05844+species(i,1));

mc0=state.species(i,7)./(sum(state.species(i,4:end),2)+state.species(i,1)*state.m_NaCl*0.05844+state.species(i,1));
mc=species(i,7)./(sum(species(i,4:end),2)+species(i,1)*state.m_NaCl*0.05844+species(i,1));
state.dcp(i)=(mc-mc0)'./dpp;
state.dsp(i)=(salinity(i)-state.salinity(i))'./dpp;

state.dfrac(i,:)=(frac(i,:)-state.frac(i,:))./dpp;
%state.dfrac(i,:)=state.dfrac(i,:)*0.6;
  state.dspecies(i)=(species(i,2)-state.species(i,2))/dpp;
% if state.pressure(i,2)<state.pressure(i,1)
%           state.species(i,7)=state.species(i,7)+state.species(i,2);
%          state.species(i,2)=0;
%          
% end
% speciesw=[state.species(i,1);0;0;0];
% state.frac(i,:)=(state.omega_t(:,4:end)*state.species(i,4:end)'+speciesw)'./(sum(state.species(i,:))-state.species(i,2)-state.species(i,3)+state.species(i,1)*state.m_NaCl*0.05844);
% state.salinity(i)=(sum(state.species(i,4:end),2)+state.species(i,1)*state.m_NaCl*0.05844)./(sum(state.species(i,4:end),2)+state.species(i,1)*state.m_NaCl*0.05844+state.species(i,1));

end
%state.poro=1-state.species(:,3)./state.rhoref(3)./V;
end





