function  [v,dv]=VCO2(p0,Tk)

% Calculate the volume  of CO2 and its gradient over pressure, at given pressure and
% temperature. Based on sovlving the Equation (B8) in ''Spycher et al. 2003 CO2-H2O mixtures
%in the geological sequestration of CO2. I. Assessment and calculation of
%mutual solubilities from 12 to 100°C and up to 600 bar'' with the method given in
%''Nickalls 1993 A new approach to solving the cubic: Cardan's solution revealed''
%
% SYNOPSIS:
%    [v,dv]=VCO2(p0,Tk)
%
% DESCRIPTION:
%   This function serves to calculate the CO2 volume V(P,T) and its
%   gradient over P, i.e., dV(P,T)/dP. The Redlich-Kwong equation [Redlich and Kwong,1949]
%   is solved based on the method given in Spycher et al. [2003]. The detail of this calculation can
%   be found in the appendix (the section called: density and viscosity) in user's guide. 
%
% REQUIRED PARAMETERS:
%  p0    - gas pressure
%  Tk     - temperature in K.
%
% RETURNS:
%   v -the volume  of CO2 , i.e., V(P,T),  
%   dv   -  the gradient of CO2 volume over P, i.e., dV(P,T)/dP.
%
% SEE ALSO:
%   EOSPP, D_RHOMU, rhomu_p_frac

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




R=83.144598;
P=p0./100000;
v=zeros(size(P));
V=zeros(size(P,1),3);
dv=zeros(size(P));
bco2=27.8; % parameter to calculate the fugacity [cm^3/mol]
aco2=7.54e7-4.13e4*Tk; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
%aco2=7.04e7;% this is from maarti, which is wrong
am=aco2; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
bm=bco2; % parameter to calculate the fugacity [cm^3/mol]
pv=[ones(size(P)), -R.*Tk./P, -(R.*Tk.*bm./P-am./P./Tk.^0.5+bm.^2), -am.*bm./P./Tk.^0.5];
E_A=pv(:,1);E_B=pv(:,2);E_C=pv(:,3);E_D=pv(:,4);
V_N=-E_B./(3.*E_A);
y_N=E_A.*V_N.^3+E_B.*V_N.^2+E_C.*V_N+E_D;
h=2.*E_A.*((E_B.^2-3.*E_A.*E_C)./(9.*E_A.^2)).^(3/2);
DIS=y_N.^2-h.^2;
r1=DIS>0;
r2=DIS==0;
r3= DIS<0;
delta=zeros(size(P));
delta(~r2)=sqrt((E_B(~r2).^2-3.*E_A(~r2).*E_C(~r2))./(9.*E_A(~r2).^2));
delta(r2)=nthroot((y_N(r2)./2./E_A(r2)),3);
theta=acos(-y_N(r3)./h(r3))./3;
V(r3,1)=V_N(r3)+2.*delta(r3).*cos(theta);
V(r3,2)=V_N(r3)+2.*delta(r3).*cos(theta+2*pi/3);
V(r3,3)=V_N(r3)+2.*delta(r3).*cos(theta+4*pi/3);
V(r2,1)=V_N(r2)+delta(r2);
V(r2,2)=V_N(r2)-2.*delta(r2);
V(r1,1)=V_N(r1)+nthroot((1./2./E_A(r1).*(-y_N(r1)+sqrt(y_N(r1).^2-h(r1).^2))),3)+nthroot((1./2./E_A(r1).*(-y_N(r1)-sqrt(y_N(r1).^2-h(r1).^2))),3);
[vg,im]=max(V,[],2);
[vl,in]=min(V,[],2);
w1=P.*(vg-vl);
w2=R.*Tk.*log((vg-bm)./(vl-bm))+am./Tk.^0.5./bm.*log((vg+bm).*vl./(vl+bm)./vg);
imin=w2<w1;
v(imin)=vl(imin)*1e-6;
v(~imin)=vg(~imin)*1e-6;
v(r1)=V(r1,1).*1e-6;
ii=w2(r3)<w1(r3);
iv=zeros(size(ii));
in3=in(r3);
im3=im(r3);
iv(ii)=in3(ii);
iv(~ii)=im3(~ii);

DthetaB=-1/3./sqrt(1-(y_N(r3)./h(r3)).^2).*(-1./h(r3).*(2./9.*E_B(r3).^2./E_A(r3).^2-E_C(r3)./3./E_A(r3))...
    +y_N(r3)./h(r3).^2.*2./3.*E_B(r3)./E_A(r3).*sqrt((E_B(r3).^2-3.*E_A(r3).*E_C(r3))./9./E_A(r3).^2));
DthetaC=-1/3./sqrt(1-(y_N(r3)./h(r3)).^2).*(-1./h(r3).*(-E_B(r3)./3./E_A(r3))...
    -y_N(r3)./h(r3).^2.*sqrt((E_B(r3).^2-3.*E_A(r3).*E_C(r3))./9./E_A(r3).^2));
DthetaD=-1/3./sqrt(1-(y_N(r3)./h(r3)).^2).*(-1./h(r3));
DV1B=-1./(3.*E_A(r3))+2.*E_B(r3)./3./abs(E_A(r3))./sqrt(E_B(r3).^2-3.*E_A(r3).*E_C(r3)).*cos(theta+(iv-1).*2/3.*pi)-2.*delta(r3).*sin(theta+(iv-1).*2/3.*pi).*DthetaB;
DV1C=-abs(E_A(r3))./E_A(r3)./sqrt(E_B(r3).^2-3.*E_A(r3).*E_C(r3)).*cos(theta+(iv-1).*2/3.*pi)-2.*delta(r3).*sin(theta+(iv-1).*2/3.*pi).*DthetaC;
DV1D=-2.*delta(r3).*sin(theta+(iv-1).*2/3.*pi).*DthetaD;
DBP=R.*Tk./P(r3).^2;
DCP=-1.*(-bm.*R.*Tk./P(r3).^2+am./sqrt(Tk)./P(r3).^2);
DDP=am.*bm./sqrt(Tk)./P(r3).^2;
dv(r3)=(DV1B.*DBP+DV1C.*DCP+DV1D.*DDP).*1e-5.*1e-6;

ii=w2(r2)<w1(r2);
iv=zeros(size(ii));
in2=in(r2);
im2=im(r2);
iv(ii)=in2(ii);
iv(~ii)=im2(~ii);
iv(iv==2)=-2;
DV1B=-1./3./E_A(r2)+iv.*nthroot((2.*E_A(r2)),-3)./3.*nthroot(y_N(r2).^(-2),3).*(2.*E_B(r2).^2./9./E_A(r2).^2-E_C(r2)./3./E_A(r2));
DV1C=iv.*nthroot((2.*E_A(r2)).^(-1),3)./3.*nthroot(y_N(r2).^(-2),3).*(-E_B(r2)./3./E_A(r2));
DV1D=iv.*nthroot((2.*E_A(r2)).^(-1),3)./3.*nthroot(y_N(r2).^(-2),3);
DBP=R.*Tk./P(r2).^2;
DCP=-1.*(-bm.*R.*Tk./P(r2).^2+am./sqrt(Tk)./P(r2).^2);
DDP=am.*bm./sqrt(Tk)./P(r2).^2;
dv(r2)=(DV1B.*DBP+DV1C.*DCP+DV1D.*DDP).*1e-5.*1e-6;
db=0.5.*(y_N(r1).^2-h(r1).^2).^(-0.5).*(2.*y_N(r1).*(2.*E_B(r1).^2./9./E_A(r1).^2-E_C(r1)./3./E_A(r1))-2.*h(r1).*2.*E_A(r1).*3./2.*((E_B(r1).^2-3.*E_A(r1).*E_C(r1))./9./E_A(r1).^2).^0.5.*2.*E_B(r1)./9./E_A(r1).^2);%(partial sqrt(y_N^2-h^2))
DV1B=-1./3./E_A(r1)+1./3.*nthroot((1./2./E_A(r1).*(-y_N(r1)+sqrt(y_N(r1).^2-h(r1).^2))).^(-2),3)./2./E_A(r1).*(-2.*E_B(r1).^2./9./E_A(r1).^2+E_C(r1)./3./E_A(r1)+db)+...
    1./3.*nthroot((1./2./E_A(r1).*(-y_N(r1)-sqrt(y_N(r1).^2-h(r1).^2))).^(-2),3)./2./E_A(r1).*(-2.*E_B(r1).^2./9./E_A(r1).^2+E_C(r1)./3./E_A(r1)-db);
dc=0.5.*(y_N(r1).^2-h(r1).^2).^(-0.5).*(2.*y_N(r1).*(-1.*E_B(r1)./3./E_A(r1))-2.*h(r1).*2.*E_A(r1).*3./2.*((E_B(r1).^2-3.*E_A(r1).*E_C(r1))./9./E_A(r1).^2).^0.5.*(-1)./3./E_A(r1));%(partial sqrt(y_N^2-h^2))
DV1C=1./3.*nthroot((1./2./E_A(r1).*(-y_N(r1)+sqrt(y_N(r1).^2-h(r1).^2))).^(-2),3)./2./E_A(r1).*(E_B(r1)./3./E_A(r1)+dc)+...
             1./3.*nthroot((1./2./E_A(r1).*(-y_N(r1)-sqrt(y_N(r1).^2-h(r1).^2))).^(-2),3)./2./E_A(r1).*(E_B(r1)./3./E_A(r1)-dc);
dd=0.5.*(y_N(r1).^2-h(r1).^2).^(-0.5).*(2.*y_N(r1));
DV1D=1./3.*nthroot((1./2./E_A(r1).*(-y_N(r1)+sqrt(y_N(r1).^2-h(r1).^2))).^(-2),3)./2./E_A(r1).*(-1+dd)+...
    1./3.*nthroot((1./2./E_A(r1).*(-y_N(r1)-sqrt(y_N(r1).^2-h(r1).^2))).^(-2),3)./2./E_A(r1).*(-1-dd);
DBP=R.*Tk./P(r1).^2;
DCP=-1.*(-bm.*R.*Tk./P(r1).^2+am./sqrt(Tk)./P(r1).^2);
DDP=am.*bm./sqrt(Tk)./P(r1).^2;
dv(r1)=(DV1B.*DBP+DV1C.*DCP+DV1D.*DDP).*1e-5.*1e-6;
end