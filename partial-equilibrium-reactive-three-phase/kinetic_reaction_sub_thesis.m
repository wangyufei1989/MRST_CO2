 function [n,h]=kinetic_reaction_sub_thesis(Pl,Pg,nu,G,G_t,n,H,Tk,Sr,dt)
 
 % speciation calculation based on the partial equilibrium reaction given in section called
%Partial Equilibrium Reactive Multi-Component Three-Phase Flow Mode in the user'guides
%
% SYNOPSIS:
%[n,h]=kinetic_reaction_sub(Pl,Pg,nu,G,G_t,n,H,Tk,Sr,dt)
%
% DESCRIPTION:
%   This function serves to calculate the mass fraction of each species and components. 
%Description is already included in the user's guide. 

%
% REQUIRED PARAMETERS:
% Pl       -liquid pressure
%Pg      - gas pressure
% nu   -stoichiometric matrix
% G        - kernel matrix based on molar abundance
% G_t      -kernel matrix based on mass
% n         -mass of the species
% H         - mass of the components
% Tk       -temperature in K
%Sr, dt   - useless for equilibrium reactions

% RETURNS:
%   n, h   -  mass of species and components

% SEE ALSO:
%  kinetic_reaction

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
 
 
%% born coefficient and effective electrostatic radius

nu_k=nu(1,:);
nu_e=nu(2:5,:);
G_e=G;
G_et=G_t;
G_e(:,[3,5])=0;
G_k=G-G_e;
G_et(:,[3,5])=0;
G_kt=G_t-G_et;
%%


Pl=Pl.*1e-5;
Pg=Pg.*1e-5;
%% Right hand side


n=n';
h=H';

n(5)=(h(3)-G_t(3,3)*n(3))/G_t(3,5);
if n(5)<0
    n(5)=1e-23;
end

he=h-G_t(:,5)*n(5);
he(3)=[];
h(1)=n(3);
h(2)=n(5);
h(3:5)=he;
nk=1;
dh=1;
MM=[0.0180150000000000,0,0,0,0,0,0,0,0;0,0.0440100000000000,0,0,0,0,0,0,0;0,0,0.100086900000000,0,0,0,0,0,0;0,0,0,0.00100780000000000,0,0,0,0,0;0,0,0,0,0.0400780000000000,0,0,0,0;0,0,0,0,0,0.0610170000000000,0,0,0;0,0,0,0,0,0,0.0440100000000000,0,0;0,0,0,0,0,0,0,0.0600090000000000,0;0,0,0,0,0,0,0,0,0.0170080000000000];

M=[55.5092978073827,0,0,0,0,0,0,0,0;0,22.7221086116792,0,0,0,0,0,0,0;0,0,9.99131754505335,0,0,0,0,0,0;0,0,0,992.260369120857,0,0,0,0,0;0,0,0,0,24.9513448774889,0,0,0,0;0,0,0,0,0,16.3888752314929,0,0,0;0,0,0,0,0,0,22.7221086116792,0,0;0,0,0,0,0,0,0,16.6641670416104,0;0,0,0,0,0,0,0,0,58.7958607714017];
n=M*n;

n=kinetic_reaction_sub_sub(Pl,Pg,h,n,G_e,G_et,nu_e,Tk,MM);

n0=n;
t=0;
v=VCO2S(Pg,Tk);
    [A,a,loga]=aloga(n,v,Pg,Tk);
 [r_k,dr_k]= rkdrk(nu_k,a,A,Tk,Pg,n,Sr);
 

% dt=1e-3
phase_3=2;
if phase_3==1% equilibrium three phase
if  dt>1e100
     %dt>1e100*min(abs(n(3:6)./r_k))
% dtt=min(0.1.*min(abs(n(3:6)./r_k)),dt-t);
 n2=1;n1=0;
 while t<dt&&max(abs(n2-n1)./n2)>1e-19
     err=1;
  nkk=1;
  n1=n;
while err>1e-12&&nkk<7
nkk=nkk+1;
%%



   
    [L,R]=LRH(dtt,n,n0,G_e,G_k,G_et,G_kt,nu_e,nu_k,r_k,dr_k,A,MM);
    dh=-L\R;

  
    h=dh+h;
    n(3)=n(3)+dh(1).*M(3,3); %
    n(5)=n(5)+dh(2).*M(5,5);
    n=kinetic_reaction_sub_sub(Pl,Pg,h,n,G_e,G_et,nu_e,Tk,MM);
   
err=max(abs(dh./h));
    [A,a,loga]=aloga(n,v,Pg,Tk);
   [r_k,dr_k]= rkdrk(nu_k,a,A,Tk,Pg,n,Sr);
end
 n2=n;
nk=nk+1;
   t=t+dtt;
    [A,a,loga]=aloga(n,v,Pg,Tk);
   [r_k,dr_k]= rkdrk(nu_k,a,A,Tk,Pg,n,Sr);
       %dtt=min(1e3.*min(abs(n(3:6)./(2.*r_k+dr_k(3:6).*dt))),dt-t);
  dtt=1e3;
 
      % dtt(dtt<1e-2)=1e-2;
      DTT(nk-1)=dtt;
     dtt=min(dtt,dt-t);
    
       n0=n;
 %test_co2(nk-1)=n(7);
  %test_calcite(nk-1)=n(3);
  %plot(1:length( test_calcite), test_calcite,'-')
  %plot(1:length(  test_co2),  test_co2)
  %plot(1:length(  DTT),  DTT)
%h=G_t*MM*n;
%% REACTION RATE 

 end
 
else
   % ntt=n;
   n=reaction_kinetic_equilibrium(Pl,Pg,nu,G_t,n,H,Tk); 
    
end
end


%%
if phase_3==2% true kinetic
if dt<1e100*min(abs(n(3:6)./r_k))
%dtt=min(0.1.*min(abs(n(3:6)./r_k)),dt-t);
%dtt=min(1e11.*min(abs(n(3:6)./(2.*r_k+dr_k(3:6).*dt))),dt-t);
 dtt=min(1e11.*min(abs(n(3:6)./(2.*r_k))),dt-t);
 n2=1;n1=0;
 while t<dt
     err=1;
  nkk=1;
  n1=n;
while err>1e-11&&nkk<10
nkk=nkk+1;
%%



   
    [L,R]=LRH(dtt,n,n0,G_e,G_k,G_et,G_kt,nu_e,nu_k,r_k,dr_k,A,MM);
    dh=-L\R;

  
    h=dh+h;
    n(3)=n(3)+dh(1).*M(3,3); %
    n(5)=n(5)+dh(2).*M(5,5);
    n=kinetic_reaction_sub_sub(Pl,Pg,h,n,G_e,G_et,nu_e,Tk,MM);
   
err=max(abs(dh./h));
    [A,a,loga]=aloga(n,v,Pg,Tk);
   [r_k,dr_k]= rkdrk(nu_k,a,A,Tk,Pg,n,Sr);
end
 n2=n;
nk=nk+1;
   t=t+dtt;
    [A,a,loga]=aloga(n,v,Pg,Tk);
   [r_k,dr_k]= rkdrk(nu_k,a,A,Tk,Pg,n,Sr);
       %dtt=min(1e11.*min(abs(n(3:6)./(2.*r_k+dr_k(3:6).*dt))),dt-t);
 dtt=min(1e11.*min(abs(n(3:6)./(2.*r_k))),dt-t);
 dtt=max(dtt,dt/5);
 dtt=min(dtt,dt-t);
      % dtt(dtt<1e-2)=1e-2;
    
       n0=n;
 %test_co2(nk-1)=n(7);
  %test_calcite(nk-1)=n(3);
  %plot(1:length( test_calcite), test_calcite,'-')
  %plot(1:length(  test_co2),  test_co2)
  %plot(1:length(  DTT),  DTT)
%h=G_t*MM*n;
%% REACTION RATE 

 end
 
else
   % ntt=n;
   n=reaction_kinetic_equilibrium(Pl,Pg,nu,G_t,n,H,Tk); 
    
end
end
%(sum(MM*ntt)-sum(MM*n))/sum(MM*n)
%n=nt;%%%%%%%%%%%%%%%%%%%%%%%%%
 n=MM*n;
h=G_t*n;
n=n';

end
function [A,a,loga]=aloga(n,v,Pg,Tk)
%% born coefficient and effective electrostatic radius
m=n./n(1).*55.51;
r_na=1.91;r_cl=1.81;r_h=3.08;r_ca=2.87; r_hco3=2.10;r_co3=2.81;r_oh=1.4;
w_na=0.3305*10^5;w_cl=1.456*10^5;w_h=0.5387*10^5;w_ca=2.314*10^5;w_hco3=0.7906*10^5;w_co3=2.3634*10^5;w_oh=1.1859*10^5;
ns=9;
z5=2;z6=-1;z8=-2;z9=-1;
%%

P0=1; %reference pressure 1 [bar]

Tc=Tk-273.15;%temperature in [kelvin]
R=83.144598; %universal gas constant []
bh2o=18.18; % parameter to calculate the fugacity [cm^3/mol]
bco2=27.8; % parameter to calculate the fugacity [cm^3/mol]
ach=7.89e7; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
aco2=7.54e7-4.13e4*Tk; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
am=aco2; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
bm=bco2; % parameter to calculate the fugacity [cm^3/mol]
%% calculate a1 and dlnaa;H2O(l)
nsm=0.5; %molality of NaCl.BE CAREFUL nsm=0, but we should not use 0. and instead we should use 0.00001>0
nsn=[n(1)*18.02/1000*nsm];%
m_na=nsm;
m_cl=nsm;
n_na=nsn;
n_cl=nsn;
I_m=0.5*(m_na*(1)^2+m_cl*(-1)^2+m(4)*(1)^2+m(5)*(2)^2+m(6)*(-1)^2+m(8)*(-2)^2+m(9)*(-1)^2);
dI_m=zeros(9,1);
dI_m(4)=0.5*1*55.51/n(1);dI_m(5)=0.5*4*55.51/n(1);dI_m(6)=0.5*1*55.51/n(1);dI_m(8)=0.5*4*55.51/n(1);dI_m(9)=0.5*1*55.51/n(1);
dI_m(1)=-(dI_m(4)*n(4)+dI_m(5)*n(5)+dI_m(6)*n(6)+dI_m(8)*n(8)+dI_m(9)*n(9))/n(1);
x_w=n(1)/(n(1)+sum(n(4:end))+n_na+n_cl);% should we consider the mole of ion. it does not matter because the number of ion is small
dx_w=zeros(ns,1);
dx_w(1)=(sum(n(4:9))+n_na+n_cl)./((sum(n(4:9))+n_na+n_cl+n(1))).^2;
dx_w(4:9)=-n(1)./((sum(n(4:9))+n_na+n_cl+n(1))).^2;
%a_gamma=input('from table B1 of allan');
a_gamma=0.5281;
%b_nacl=input('from table B.3 of allan');
b_nacl=4.783e-7;
%b_na_cl=input('from table B.4 of allan');
b_na_cl=-5.603e-2;
%b_gamma=input('from table B2 of allan');
b_gamma=0.5281e8;



a_o=(m_na*r_na+m_cl*r_cl)/(m_na+m_cl);
AA=1+a_o*b_gamma*(I_m)^0.5;
sigma=3/(a_o*b_gamma*(I_m)^0.5)^3*(AA-1/AA-2*log(AA));
psi_na=a_gamma*(1)^2*(I_m)^0.5/3*sigma+x_w/(1-x_w)*log10(x_w)-0.5*(w_na*b_nacl+b_na_cl-0.19*(1-1))*I_m;
psi_cl=a_gamma*(-1)^2*(I_m)^0.5/3*sigma+x_w/(1-x_w)*log10(x_w)-0.5*(w_cl*b_nacl+b_na_cl-0.19*(1-1))*I_m;
loga1=2.303/55.508*(m_na*psi_na+m_cl*psi_cl);
a1=exp(loga1);
dx_w_2=(n(4)+n_na+n_cl)/(n(2)+n(4)+n_na+n_cl)^2;
dx_w_4=-n(2)/(n(2)+n(4)+n_na+n_cl)^2;

da1_2=2.303/55.508*((dx_w_2/(1-x_w)^2*log10(x_w)+x_w/(1-x_w)*log10(2.7183)/x_w*dx_w_2)*m_na+...
    (dx_w_2/(1-x_w)^2*log10(x_w)+x_w/(1-x_w)*log10(2.7183)/x_w*dx_w_2)*m_cl);

da1_4=2.303/55.508*((dx_w_4/(1-x_w)^2*log10(x_w)+x_w/(1-x_w)*log10(2.7183)/x_w*dx_w_4)*m_na+...
    (dx_w_4/(1-x_w)^2*log10(x_w)+x_w/(1-x_w)*log10(2.7183)/x_w*dx_w_4)*m_cl);

da1=zeros(ns,1);
da1(1)=da1_2; da1(4:end)=da1_4;
 %% calculate a2 and dlna2;CO2(g)
 logphi2=log(v/(v-bm))+bco2/(v-bm)+am*bco2/R/Tk^1.5/bm^2*(log((v+bco2)/v)...
-bm/(v+bm))-2*aco2/R/Tk^1.5/bm*log((v+bm)/v)-log(Pg*v/R/Tk);

phi2=exp(logphi2);

a2=phi2*Pg/P0;
da2=zeros(9,1);

%% calculate a3 and dloga3;CaCO3(s)
a3=1;
da3=zeros(9,1);
%% calculate a4 and dlna4;H+

a_o=(m_na*r_na+m_cl*r_cl+m(4)*r_h+m(5)*r_ca+m(6)*r_hco3+m(8)*r_co3+m(9)*r_oh)/(m_na+m_cl+m(4)+m(5)+m(6)+m(8)+m(9));
AA=1+a_o*b_gamma*(I_m)^0.5;

log10a4=-a_gamma*(1)^2*sqrt(I_m)/AA+log10(x_w)+(w_h*b_nacl+b_na_cl-0.19*(abs(1)-1))*I_m;
a4=m(4)*10^(log10a4);
d_AA=a_o.*b_gamma.*0.5./sqrt(I_m).*dI_m;
da4=log(10)*((-a_gamma*(1)^2).*(0.5./sqrt(I_m).*dI_m.*AA-d_AA.*sqrt(I_m))./AA.^2+1./log(10)./x_w.*dx_w+(w_h*b_nacl+b_na_cl-0.19*(abs(1)-1))*dI_m);
da4(4)=da4(4)+1/n(4);
da4(1)=da4(1)-1/n(1);
%% calculate a5 and dlna5;ca2+
log10a5=-a_gamma*(z5)^2*sqrt(I_m)/AA+log10(x_w)+(w_ca*b_nacl+b_na_cl-0.19*(abs(z5)-1))*I_m;
a5=m(5)*10^(log10a5);
da5=log(10)*((-a_gamma*(z5)^2).*(0.5./sqrt(I_m).*dI_m.*AA-d_AA.*sqrt(I_m))./AA.^2+1./log(10)./x_w.*dx_w+(w_ca*b_nacl+b_na_cl-0.19*(abs(z5)-1))*dI_m);
da5(5)=da5(5)+1/n(5);
da5(1)=da5(1)-1/n(1);
%% calculate a6 and dlna6;hcao3-
log10a6=-a_gamma*(z6)^2*sqrt(I_m)/AA+log10(x_w)+(w_hco3*b_nacl+b_na_cl-0.19*(abs(z6)-1))*I_m;
a6=m(6)*10^(log10a6);
da6=log(10)*((-a_gamma*(z6)^2).*(0.5./sqrt(I_m).*dI_m.*AA-d_AA.*sqrt(I_m))./AA.^2+1./log(10)./x_w.*dx_w+(w_hco3*b_nacl+b_na_cl-0.19*(abs(z6)-1))*dI_m);
da6(6)=da6(6)+1/n(6);
da6(1)=da6(1)-1/n(1);
%% calculate a7 and dlna7;co2(aq)
ap1=-1.0312;
ap2=1.2806e-3;
ap3=255.9;
ap4=0.4445;
ap5=-1.606e-3;
I=0.5*(m_na*(1)^2+m_cl*(-1)^2);
loggamma7=(ap1+ap2*Tk+ap3/Tk)*I-(ap4+ap5*Tk)*I/(I+1);
gamma7=exp(loggamma7);
a7=55.508*n(7)/n(1)*gamma7;
da7=zeros(ns,1);
da7(1)=-1./n(1); da7(7)=1./n(7);

%% calculate a8 and dlna8;co3--
log10a8=-a_gamma*(z8)^2*sqrt(I_m)/AA+log10(x_w)+(w_co3*b_nacl+b_na_cl-0.19*(abs(z8)-1))*I_m;
a8=m(8)*10^(log10a8);
da8=log(10)*((-a_gamma*(z8)^2).*(0.5./sqrt(I_m).*dI_m.*AA-d_AA.*sqrt(I_m))./AA.^2+1./log(10)./x_w.*dx_w+(w_co3*b_nacl+b_na_cl-0.19*(abs(z8)-1))*dI_m);
da8(8)=da8(8)+1/n(8);
da8(1)=da8(1)-1/n(1);
%% calculate a9 and dlna9;oh-
log10a9=-a_gamma*(z9)^2*sqrt(I_m)/AA+log10(x_w)+(w_oh*b_nacl+b_na_cl-0.19*(abs(z9)-1))*I_m;
a9=m(9)*10^(log10a9);
da9=log(10).*((-a_gamma*(z9)^2).*(0.5./sqrt(I_m).*dI_m.*AA-d_AA.*sqrt(I_m))./AA.^2+1./log(10)./x_w.*dx_w+(w_oh*b_nacl+b_na_cl-0.19*(abs(z9)-1)).*dI_m);
da9(9)=da9(9)+1/n(9);
da9(1)=da9(1)-1/n(1);

A=[da1';da2';da3';da4';da5';da6';da7';da8';da9'];
a=[a1;a2;a3;a4;a5;a6;a7;a8;a9];
loga=log([a1;a2;a3;a4;a5;a6;a7;a8;a9]);
% id=(n<10^(-16))&(n>0);
% A(id,id)=1./n(id);


end
function [L,R]= LRH(dt,n,n0,G_e,G_k,G_et,G_kt,nu_e,nu_k,r_k,dr_k,A,MM)
inde=true(9,1);
inde([3 5])=false;
%MM=[0.0180150000000000,0,0,0,0,0,0,0,0;0,0.0440100000000000,0,0,0,0,0,0,0;0,0,0.100086900000000,0,0,0,0,0,0;0,0,0,0.00100780000000000,0,0,0,0,0;0,0,0,0,0.0400780000000000,0,0,0,0;0,0,0,0,0,0.0610170000000000,0,0,0;0,0,0,0,0,0,0.0440100000000000,0,0;0,0,0,0,0,0,0,0.0600090000000000,0;0,0,0,0,0,0,0,0,0.0170080000000000];
h_e=G_et*MM*n;
h_e(3)=[];
h_k=[n(3)*MM(3,3);n(5)*MM(5,5)];
h_e0=G_et*MM*n0;
h_e0(3)=[];
h_k0=[n0(3)*MM(3,3);n0(5)*MM(5,5)];
R=[h_k-h_k0;h_e-h_e0]./dt;
g_et=G_et;
g_et(3,:)=[];
ge=G_e;
ge(3,:)=[];
R=R-[nu_k(3)'*MM(3,3);nu_k(5)'*MM(5,5);g_et*MM*nu_k']*r_k;
%% jacobian
L=eye(5,5)./dt;
drK=zeros(1,5);
nu_e(:,[3 5])=[];
A([3 5],:)=[];
A(:,[3 5])=[];
ge1=ge;
ge1(:,[3,5])=[];
g_et1=g_et;
g_et1(:,[3,5])=[];
mm=MM;
mm(:,[3 5])=[];
mm([3 5],:)=[];
j1=[nu_e*A;g_et1*mm];
f1=zeros(7,3);
f1(5:7,:)=eye(3,3);
dne_h=j1\f1;
%drK(1)=-dr_k(5);
drK(2)=dr_k(5)./MM(5,5);

drK(3:5)=dr_k(inde)'*dne_h;
L=L-[MM(3,3)*nu_k(:,3)';MM(5,5)*nu_k(:,5)';g_et*MM*nu_k']*drK;%
end


function [r_k,dr_k]= rkdrk(nu_k,a,A,Tk,Pg,n,AS)

k1=10^(10.329-8.48);
I_s=exp(nu_k*log(a)-log(k1));
[rk10,rk20,rk30]=deal(10^(-0.3),10^(-5.81),10^(-3.48));
[Ea1,Ea2,Ea3]=deal(14.4e3,23.5e3,35.4e3);
R=8.314; %in j/K/mol
rk1=rk10*exp(-Ea1/R*(1./Tk-1/298.15))*a(4); %rk10 in mol/m^2/s 
rk2=rk20*exp(-Ea2/R*(1./Tk-1/298.15));%Ea in j
rk3=rk30*exp(-Ea3/R*(1./Tk-1/298.15))*Pg;%pressure in bar.
r_k=AS*(1-I_s)*(rk1+rk2+rk3);%AS is the surface area. in m^2
dr_k=zeros(9,1);
dr_k(4)=AS*(1-I_s)*rk10*exp(-Ea1/R*(1./Tk-1/298.15))./n(1).*55.51;

% nu_kk=nu_k;
% nu_kk([1 2 3 4 6 7 8 9])=0;
% dr_k(5)=AS*(rk1+rk2+rk3).*(-exp(nu_k*log(a)-log(k1))).*nu_k*(A*(nu_kk'));
% nu_kk=nu_k;
% nu_kk([1 2 3 5 6 7 8 9])=0;
% dr_k(4)=dr_k(4)+AS*(rk1+rk2+rk3).*(-exp(nu_k*log(a)-log(k1))).*nu_k*(A*(nu_kk'));
% nu_kk=nu_k;
% nu_kk([1 2 3 4 5 7 8 9])=0;
% dr_k(6)=AS*(rk1+rk2+rk3).*(-exp(nu_k*log(a)-log(k1))).*nu_k*(A*(nu_kk'));
% 
% dr_k=zeros(9,1);
nu_kk=ones(size(nu_k));
nu_kk([1 2 3 4 6 7 8 9])=0;
dr_k(5)=AS*(rk1+rk2+rk3).*(-exp(nu_k*log(a)-log(k1))).*nu_k*(A*(nu_kk'));
nu_kk=ones(size(nu_k));
nu_kk([1 2 3 5 6 7 8 9])=0;
dr_k(4)=dr_k(4)+AS*(rk1+rk2+rk3).*(-exp(nu_k*log(a)-log(k1))).*nu_k*(A*(nu_kk'));
nu_kk=ones(size(nu_k));
nu_kk([1 2 3 4 5 7 8 9])=0;
dr_k(6)=AS*(rk1+rk2+rk3).*(-exp(nu_k*log(a)-log(k1))).*nu_k*(A*(nu_kk'));

end

function v=VCO2S(Pg,Tk)
%% CO2 volume
P0=1; %reference pressure 1 [bar]
R=83.144598; %universal gas constant []
bh2o=18.18; % parameter to calculate the fugacity [cm^3/mol]
bco2=27.8; % parameter to calculate the fugacity [cm^3/mol]
ach=7.89e7; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
aco2=7.54e7-4.13e4*Tk; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
am=aco2; % parameter to calculate the fugacity [bar*cm^6*K^0.5/mol]
bm=bco2; % parameter to calculate the fugacity [cm^3/mol]
    pv=[1, -R*Tk/Pg, -(R*Tk*bm/Pg-am/Pg/Tk^0.5+bm^2), -am*bm/Pg/Tk^0.5];
v=roots(pv);
v=v(imag(v)==0);
v=sort(v);
vg=v(end);
vl=v(1);
w1=Pg*(vg-vl);
w2=R*Tk*log((vg-bm)/(vl-bm))+am/Tk^0.5/bm*log((vg+bm)*vl/(vl+bm)/vg);
if w2<w1
    v=vl;
else
    v=vg;
end
end
