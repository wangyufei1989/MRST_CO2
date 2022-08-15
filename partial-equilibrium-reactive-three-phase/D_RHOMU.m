function [d_rho,d_mu]=D_RHOMU(state)

% Calculate the density gradient and the viscosity gradient of brine and co2 phase at
% given pressure, mass fraction, NaCl concentration and temperature; i.e.,
% drho/dp, dmu/dp.
%
% SYNOPSIS:
%  [d_rho,d_mu]=D_RHOMU(state)
%
% DESCRIPTION:
% This function serves to calculate the  density gradient and viscosity gradient of brine and co2 phase at
% given pressure, mass fraction, NaCl concentration and temperature.
%The following model for brine denisty can be found in Garcia [2001;2003].
%the following model for brine viscosity can be found in Garcia [2003].
%the following model for CO2 density can be found in Victor Vilarrasa [2012] phd thesis
%or Spycher et al. [2003].
%the following model for CO2 viscosity can be found in Victor Vilarrasa [2012] phd thesis.
%Description is already included in the user's guide.
%
% REQUIRED PARAMETERS:
% state   -state varibles
% RETURNS:
%   d_rho   -density gradient of brine and gas over pressure- [drho_l/dp_l,drho_g/dp_g,drho_l/dp_g,drho_g/dp_l];
%   d_mu    -  viscosity gradient of brine and gas over pressure-[dmu_l/dp_l,dmu_g/dp_g,dmu_l/dp_g,dmu_g/dp_l];  .
%
% SEE ALSO:
% EOSPP, rhomu_p_frac

%{

This file is a part of mrst_co2 based on MRST.

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
%%
% the brine density derivative
if size(state.pressure,2)==1
    state.pressure=[state.pressure,state.pressure];
end
p0=state.pressure(:,1);
if isempty(state.frac) %for wells
    state.frac=repmat([1,0,0,1],[size(p0,1),1]);
end
if isfield(state,'m_NaCl')
    m_NaCl=state.m_NaCl;%
else
    m_NaCl=0;
end
Mc=0.044; %molar mass of co2
T_c=state.Tk-273.15;
T=state.Tk;
vphi=3.751*10^(-5)-9.585*10^(-8)* T_c+8.74*10^(-10)* T_c^2-5.044*10^(-13)* T_c^3;
x=-9.9595 *exp( -4.539*10^(-3)*m_NaCl )+7.0845*exp( -1.638*10^(-4)*T_c)...
    + 3.9093*exp(2.551*10^(-10).*p0);
rhor=-3.033405*10^3+1.0128163*10^4 *x -8.750567*10^3* x.^2+2.66310*10^3* x.^3;
dx_p=3.9093*exp(2.551*10^(-10).*p0).*2.551*10^(-10);
drhor_p=1.0128163*10^4 *dx_p-2*8.750567*10^3.* x.*dx_p+3*2.66310*10^3* x.^2.*dx_p;
drhob_p=drhor_p.*(1+state.frac(:,2).*(1-rhor.*vphi./Mc))+rhor.*(-state.frac(:,2).*drhor_p.*vphi./Mc)+rhor.*(1-rhor.*vphi./Mc).*state.dfrac(:,1);
drhob_pg=rhor.*(1-rhor.*vphi./Mc).*state.dfrac(:,3);
maarti=0;
if maarti==1
    drhob_p=1037.12.*exp(-3.4e-4.*T_c+4.5e-10.*(p0(:,1)-0.1e6)+0.19.*0.07).*4.5e-10;
    drhob_pg=zeros(size(drhob_pg));
end
% now it is the co2 density derivative
frac=state.frac;
p0=state.pressure(:,2);
nordbotten=0;%
[v,dv]=VCO2(p0,T);
if maarti==1
    [v,dv]=VCO2M(p0,T);
end
xgC=frac(:,4)./0.044./(frac(:,4)./0.044+frac(:,3)./0.018);
xgH=1-xgC;
simple=1;
if simple==1
    drhoc_p=-(0.044)./v.^2.*dv;
    drhoc_pl=zeros(size(drhoc_p));
else
    drhoc_p=-(0.044.*xgC+0.018.*xgH)./v.^2.*dv+0.026./v.*(-44*18./(18+26.*frac(:,3)).^2.*state.dfrac(:,4));
    drhoc_pl=0.026./v.*(-44*18./(18+26.*frac(:,3)).^2.*state.dfrac(:,2));
end
%derivative of viscosity of water
Mc=0.044; %molar mass of co2
m_C=frac (:,2)./Mc./frac(:,1); % molarity of co2
dMub_p=1.69956*10^(-9)*10^(-3).*ones(size(state.pressure,1),1)+((8.79552-3.17229*10^(-2)* T)+(-7.22796+2.64498*10^(-2) *T)* m_C.*2).*1e-3./0.044./(1-frac(:,2)).^2.*state.dfrac(:,1);
dMub_pg=((8.79552-3.17229*10^(-2)* T)+(-7.22796+2.64498*10^(-2) *T)* m_C.*2).*1e-3./0.044./(1-frac(:,2)).^2.*state.dfrac(:,3);
if maarti==1
    dMub_p=zeros(size(dMub_p));
    dMub_pg=zeros(size(dMub_pg));
end
%viscosity of co2
tc=T/304;%temperature scaled by critical temperature.
a10  = 0.248566120;
a11  = 0.004894942;
a20  = -0.373300660;
a21  = 1.22753488;
a30  = 0.363854523;
a31 = -0.774229021;
a40  = -0.0639070755;
a41  = 0.142507049;
rhoc=state.rho(:,2)./468;%scaled by critical density
drhoc=drhoc_p./468;
dMuc_p=state.mu(:,2).*(a10.*drhoc+a11*drhoc./tc+2.*a20.*rhoc.*drhoc+2.*a21*rhoc.*drhoc./tc...
    +3.*a30.*rhoc.^2.*drhoc+3.*a31*rhoc.^2.*drhoc./tc+4.*a40*rhoc.^3.*drhoc+4*a41*rhoc.^3.*drhoc./tc);
if simple==1
    dMuc_pl=zeros(size(dMuc_p));
else
    drhoc=drhoc_pl./468;
    dMuc_pl=state.mu(:,2).*(a10.*drhoc+a11*drhoc./tc+2.*a20.*rhoc.*drhoc+2.*a21*rhoc.*drhoc./tc...
        +3.*a30.*rhoc.^2.*drhoc+3.*a31*rhoc.^2.*drhoc./tc+4.*a40*rhoc.^3.*drhoc+4*a41*rhoc.^3.*drhoc./tc); 
end
d_mu=[dMub_p,dMuc_p,dMub_pg,dMuc_pl];
d_rho=[drhob_p,drhoc_p,drhob_pg,drhoc_pl];
if nordbotten==1
    d_mu=zeros(size(d_mu));
    d_rho=zeros(size(d_rho));
end
end