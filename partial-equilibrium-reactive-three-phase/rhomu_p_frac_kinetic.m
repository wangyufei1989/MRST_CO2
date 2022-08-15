function [ rho,mu]=rhomu_p_frac_kinetic(state)
% Calculate the  density and viscosity of brine and co2 phase at 
% given pressure, mass fraction, NaCl concentration and temperature.
%
% SYNOPSIS:
%   [ rho,mu]=rhomu_p_frac_kinetic(state)
%
% DESCRIPTION:
%   This function serves to calculate the  density and viscosity of brine and co2 phase at 
% given pressure, mass fraction, NaCl concentration and temperature. 
%The following model for brine denisty can be found in Garcia [2001;2003]. 
%the following model for brine viscosity can be found in Garcia [2003]. 
%the following model for CO2 density can be found in Victor Vilarrasa [2012] phd thesis 
%or Spycher et al. [2003]. 
%the following model for CO2 viscosity can be found in Victor Vilarrasa [2012] phd thesis. 
%Description is already included in the user's guide. 

%
% REQUIRED PARAMETERS:
%  p0    - [liquid pressure, gas pressure];
%  frac   -mass fractions of [water in brine,  gas in brine, water in gas,  co2 in gas];
% m_NaCl   - molality of NaCl;
%  Tk     - temperature in K.


% RETURNS:
%   rho -density of brine and gas,  
%
%   mu   -  viscosity of brine and gas.
%
% SEE ALSO:
%    D_RHOMU_kinetic

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
% first it is brine density
if size(state.pressure,2)==1
    p0=[state.pressure,state.pressure];
else
    p0=state.pressure;
end
benchmark=0;
Tk=state.Tk;
T_c=Tk-273.15;
rhor=zeros(size(p0,1),2);
rho=zeros(size(p0,1),3);
%p0(:,1)=max(p0(:,1),1e7);

m_NaCl=state.m_NaCl;
x=-9.9595 *exp( -4.539*10^(-3)*m_NaCl )+7.0845*exp( -1.638*10^(-4).*T_c)...
+ 3.9093*exp(2.551*10^(-10).*(p0(:,1)));

rhor(:,1)=-3.033405*10^3+1.0128163*10^4 *x -8.750567*10^3* x.^2+2.66310*10^3* x.^3;
if isempty(state.species) %for wells
    mc=repmat([0],[size(p0,1),1]);
else
   mc=state.species(:,7)./(sum(state.species(:,4:end),2)+state.species(:,1)*state.m_NaCl*0.05844+state.species(:,1));
end
Mc=0.044; %molar mass of co2
vphi=3.751*10^(-5)-9.585*10^(-8).* T_c+8.74*10^(-10).* T_c^2-5.044*10^(-13).* T_c.^3;
rho(:,1)=rhor(:,1).*(1+mc.*(1-rhor(:,1).*vphi./Mc));



% now it is the co2 density

P=p0(:,2);
[v,~]=VCO2(P,state.Tk);    
     rho(:,2)=0.044./v;
rho(:,3)=state.rhoref(3);




%viscosity of water


Mc=0.044; %molar mass of co2
m_C=mc./Mc; % molarity of co2

mu=zeros(size(p0,1),2);

mu_w0=0.6527; % [pa*s]viscosity of water at (T,p_b=10^5 pa);here 
mu(:,1)=((3.85971-1.32561*10^(-2)* Tk)* m_NaCl+(-5.37539+1.90621*10^(-2)* Tk)*...
m_NaCl.^(1/2)+(8.79552-3.17229*10^(-2)* Tk)* m_C+(-7.22796+2.64498*10^(-2) *Tk)* m_C.^2+1.69956*10^(-9)...
*(p0(:,1)-10^5 )+mu_w0).*10^(-3); 



%viscosity of co2
tc=Tk/304;%temperature scaled by critical temperature.
mu0=(tc)^0.5*(27.2246461-16.6346068/(tc) +4.66920556/(tc )^2).*10^(-6);
a10  = 0.248566120;
a11  = 0.004894942;
a20  = -0.373300660;
a21  = 1.22753488;
a30  = 0.363854523;	
a31 =- 0.774229021;
a40  = -0.0639070755;
a41  = 0.142507049;
rhoc=rho(:,2)./468;%scaled by critical density
mu(:,2)=mu0*exp(a10*rhoc+a11*rhoc./tc+a20*rhoc.^2+a21*rhoc.^2./tc...
    +a30*rhoc.^3+a31*rhoc.^3./tc+a40*rhoc.^4+a41*rhoc.^4./tc);

if benchmark==1
    rho=[ones(size(rho,1),1).*1099,ones(size(rho,1),1).*733,ones(size(rho,1),1).*2700];
    mu=[ones(size(rho,1),1).*0.5e-3,ones(size(rho,1),1).*0.061e-3];
end
 
end