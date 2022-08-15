function state=EOSPP(state)

% Calculate the mutual solubility reaction based on Spycher and Pruess
% (2005), ''CO2-H2O mixtures in the geological sequestration of CO2. II. Partitioning in chloride
%brines at 12–100°C and up to 600 bar''
% SYNOPSIS:
%  state=EOSPP(state)
%
% DESCRIPTION:
% This function serves to calculate the mass composition of CO2-rich phase
% and brine phase baed on Spycher and Pruess (2005);
%Description is already included in the user's guide.
%
% REQUIRED PARAMETERS:
% state   -state varibles;
% RETURNS:
%   state   -state variable with updated mass fractions;
% SEE ALSO:
%   VCO2

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

R=83.144598;
m_NaCl=state.m_NaCl;
p0=1;
p=state.pressure(:,2)/100000;
p1=state.pressure(:,1)/100000;
[v,dv]=VCO2(state.pressure(:,2),state.Tk);
[phi1,phi2,dphi1,dphi2]=PHI(v,dv,state.pressure(:,2),state.Tk);%%be careful about the unit of pressure bar or pa 
tc=state.Tk-273.15;
tk=state.Tk;
a1=-2.209;b1=3.097*10^(-2);c1=-1.098*10^(-4);d1=2.048*10^(-7);v1=18.1;
k10=10^(a1+b1*tc+c1*tc.^2+d1*tc.^3);
A=k10./phi1./p.*exp((p1-p0).*v1./R./tk);
a1=1.189;b1=1.304*10^(-2);c1=-5.446*10^(-5);d1=0;v2=32.6;
k20=10^(a1+b1*tc+c1*tc.^2+d1*tc.^3);
B=phi2.*p./55.508./k20.*exp(-(p1-p0).*v2./R./tk);
lambda=-0.411370585+6.07632013e-4*tk+97.5347708./tk-0.0237622469.*(p1)./tk...
    +0.0170656236.*(p1)./(630-tk)+1.41335834e-5.*tk.*log(p1);
xi=3.36389723e-4-1.98298980e-5.*tk+2.12220830e-3.*p1./tk-5.24873303e-3.*p1./(630-tk);
gamma=exp(2.*lambda.*m_NaCl+xi.*m_NaCl.*m_NaCl);
xcH0=(1-B)./(1./A-B);
xbC0=B.*(1-xcH0);
mc0=55.508.*xbC0./(1-xbC0);
mc=mc0./gamma;
xbC=mc./(mc+55.508+2*m_NaCl);
xbS=2*m_NaCl./(55.508+2*m_NaCl+mc);
xcH=A.*(1-xbC-xbS);
state.frac(:,3)=18*xcH./(18*xcH+44*(1-xcH));
state.frac(isnan(state.frac(:,3)),3)=0;
state.frac(:,4)=1-state.frac(:,3);
state.frac(:,2)=44.01*xbC./(18.015*(1-xbC).*(1+0.05844*m_NaCl)+44.01*xbC);
state.frac(isnan(state.frac(:,2)),2)=0;
state.frac(:,1)=(1-state.frac(:,2))./(1+0.05844.*state.m_NaCl);

%% derivative over p_g
k1=k10.*exp((p1-p0).*v1./R./tk);
df1pp=phi1+dphi1.*p;
dA=(-k1.*df1pp)./(phi1.*p).^2;
k2=k20.*exp((p1-p0).*v2./R./tk);
df2pp=phi2+dphi2.*p;
dB=1/55.508.*(k2.*df2pp)./(k2).^2;
dAB=1/55.508.*((dphi2.*k1).*phi1.*k2-(dphi1.*k2).*phi2.*k1)./(phi1.*k2).^2;
dxgH0=((dA-dAB).*(1-A.*B)+dAB.*(A-A.*B))./(1-A.*B).^2;
dxlC0=dB.*(1-xcH0)-B.*dxgH0;
dmlC0=55.508./(1-xbC0).^2.*dxlC0;
dmlC=(dmlC0.*gamma)./gamma.^2;
dxlC=(55.508+2*m_NaCl)./(mc+55.508+2.*m_NaCl).^2.*dmlC;   %24
dxlS=(-2*m_NaCl)./(mc+55.508+2.*m_NaCl).^2.*dmlC;   %25
dxgH=dA.*(1-xbC-xbS)-A.*(dxlC+dxlS); %26
dXgH=(18.015*44.01)./(44.01-25.995.*xcH).^2.*dxgH; %28
dXlC=dxlC.*18.015.*44.01.*(1+0.05844.*m_NaCl)./(18.015.*(1-xbC).*(1+0.05844.*m_NaCl)+44.01.*xbC).^2;
%dXlC=dxlC.*44.01.*18.015.*(1+0.05844.*m_NaCl)./(18.015.*(1-xbC).*(1+0.05844.*m_NaCl)+44.01.*xbC).^2;
state.dfrac(:,3)=dXlC./100000;
state.dfrac(:,4)=dXgH./100000;
%% derivative over pL
dk1p1=k1.*v1./R./tk;
dk2p1=k2.*v2./R./tk;
dgamma=gamma.*((-0.0237622469./tk...
    +0.0170656236./(630-tk)+1.41335834e-5.*tk./(p1)).*2.*m_NaCl+m_NaCl.^2.*(2.12220830e-3./tk-5.24873303e-3./(630-tk)));
dA=dk1p1.*phi1.*p./(phi1.*p).^2;
dB=1/55.508.*(-phi2.*p.*(dk2p1))./(k2).^2;
dAB=1/55.508.*((dk1p1.*phi2).*phi1.*k2-(dk2p1.*phi1).*phi2.*k1)./(phi1.*k2).^2;
dxgH0=((dA-dAB).*(1-A.*B)+dAB.*(A-A.*B))./(1-A.*B).^2; %34
dxlC0=dB.*(1-xcH0)-B.*dxgH0;  %35
dmlC0=55.508./(1-xbC0).^2.*dxlC0;
dmlC=(dmlC0.*gamma-dgamma.*mc0)./gamma.^2;%37
dxlC=(55.508+2*m_NaCl)./(mc+55.508+2.*m_NaCl).^2.*dmlC;   %
dxlS=(-2*m_NaCl)./(mc+55.508+2.*m_NaCl).^2.*dmlC;   %
dxgH=dA.*(1-xbC-xbS)-A.*(dxlC+dxlS); %
dXgH=(18.015*44.01)./(44.01-25.995.*xcH).^2.*dxgH; %
dXlC=dxlC.*44.01*18.015*(1+0.05844.*m_NaCl)./(18.015.*(1-xbC).*(1+0.05844.*m_NaCl)+44.01.*xbC).^2;
state.dfrac(:,1)=dXlC/100000;
state.dfrac(:,2)=dXgH/100000;
state.dfrac(isnan(state.dfrac))=0;
%%
% state.frac(:,3)=state.frac(:,2);state.frac(:,4)=1-state.frac(:,3);
%state.dfrac(:,2)=state.dfrac(:,1);state.dfrac(:,4)=state.dfrac(:,3);
simple=0;
if simple==1
    state.frac(:,3)=0;
     state.frac(:,4)=1;
     state.dfrac(:,2)=0;
     state.dfrac(:,4)=0;
end
end



function [phi1,phi2,dphi1,dphi2]=PHI(v,dv,p,Tk)
P=p./100000;
v=v.*1e6;
dv=dv.*1e5.*1e6;
R=83.144598;
bh2o=18.18; % parameter to calculate the fugacity [cm^3/mol]
bco2=27.8; % parameter to calculate the fugacity [cm^3/mol]
ach=7.89e7; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
aco2=7.54e7-4.13e4*Tk; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
%aco2=7.04e7;
am=aco2; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
bm=bco2; % parameter to calculate the fugacity [cm^3/mol]
logphi1=log(v./(v-bm))+bh2o./(v-bm)+am.*bh2o./R./Tk.^1.5./bm.^2.*(log((v+bm)./v)...
-bm./(v+bm))-2*ach./R./Tk.^1.5./bm*log((v+bm)./v)-log(P.*v./R./Tk);
phi1=exp(logphi1);
logphi2=log(v./(v-bm))+bco2./(v-bm)+am*bco2/R/Tk^1.5/bm^2.*(log((v+bm)./v)...
-bm./(v+bm))-2*aco2/R./Tk^1.5/bm.*log((v+bm)./v)-log(P.*v/R/Tk);
phi2=exp(logphi2);
dphi1=phi1.*(dv./v-dv./(v-bm)-bh2o./(v-bm).^2.*dv+am.*bh2o./R./Tk.^1.5./bm.^2.*(dv./(v+bm)-dv./v+bm./(v+bm).^2.*dv)-...
    2.*ach./R./Tk.^1.5./bm.*(dv./(v+bm)-dv./v)-(v+P.*dv)./P./v);
dphi2=phi2.*(dv./v-dv./(v-bm)-bco2./(v-bm).^2.*dv+am.*bco2./R./Tk.^1.5./bm.^2.*(dv./(v+bm)-dv./v+bm./(v+bm).^2.*dv)-...
    2.*aco2./R./Tk.^1.5./bm.*(dv./(v+bm)-dv./v)-(v+P.*dv)./P./v);
end
