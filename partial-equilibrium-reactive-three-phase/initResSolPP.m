function resSol = initResSolPP(G, p0, s0,frac0, rhoref,muref,fluid,m_NaCl,Tk,rock,property,varargin)
% initialize the system.
%
% SYNOPSIS:
% resSol = initResSolPP(G, p0, s0,frac0, rhoref,muref,fluid,m_NaCl,Tk,rock,property,varargin)
%
% DESCRIPTION:
%   This function serves to initialize the system.
%
% REQUIRED PARAMETERS:
%   G  - grid system
%  p0  - initial pressure at the top of the aquifer
%   s0  - initial saturation at the top of the aquifer
% frac0  - intial mass fraction at the top of the aquifer
%rhoref, muref  - reference density and viscosity; not used if the density
%and viscosity change with pressure.
%fluid  - relative permeability and retention curves;
% m_NaCl,Tk  -  (constant) NaCl concentration and (constant) temperature
%rock  - an object including the porosity and intrinsic permeability
%property  - indicates what phases and reactions do we have.


% RETURNS:
%   reSol   - an object containing the initialized variables.
%
% SEE ALSO:
%   initStatePP

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

%the following line is commented by Y.
%    if nargin <= 2, s0 = 0; end
[nc, nf] = deal(G.cells.num, G.faces.num);
if size(frac0,1) == 1
    frac0 = repmat(frac0, [nc, 1]);
end
if size(rhoref,1) == 1
    rhoref0 = repmat(rhoref, [nc, 1]);
end
if size(muref,1) == 1
    muref0 = repmat(muref, [nc, 1]);
end

if numel(p0) == 2
    p00 = repmat(p0, [nc, 1]);
end
if size(s0, 1) == 1 
    s0 = repmat(s0(1,:), [nc, 1]);  
elseif size(s0, 1) ~= nc    
    error(msgid('InitSat:Inconsistent'), ...
        ['Initial saturation must either be 1-by-np ', ...
        'or (G.clls.num (=%d))-by-np.'], nc);   
end
state.Tk=Tk;
if (p0(:,2)>p0(:,1)) && property==4 % two phases; compressible soluble& reaction with rock
    state.pressure=p00;
    s0=fluid.S(state);
    state.m_NaCl=m_NaCl;
    state.frac=frac0;
    state.poro=rock.poro;
    state.poro0=rock.poro;
    for i=1:(G.cartDims(3)-1)
        state=rock_frac(state);
        frac0=state.frac;
        p00=state.pressure;
        [rho0,mu0]=rhomu_p_frac_rock(p00,frac0,m_NaCl);
        rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),1)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),1)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),2)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),2)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state=rock_frac(state);
        frac0=state.frac;
        p00=state.pressure;
        [rho0,mu0]=rhomu_p_frac_rock(p00,frac0,m_NaCl);
        s0=fluid.S(state);
        rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),1)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),1)...
            +0.5.*(rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))+rhoeff0(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2)))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),2)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),2)...
            +0.5.*(rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))+rhoeff0(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2)))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
    end
    state=rock_frac(state);
    frac0=state.frac;
    p00=state.pressure;
    s0=fluid.S(state);
    [rho0,mu0]=rhomu_p_frac_rock(p00,frac0,m_NaCl);
end

if (p0(:,2)>p0(:,1)) && property==3 % two phases; compressible soluble
    state.pressure=p00;
    s0=fluid.S(state);
    state.m_NaCl=m_NaCl;
    state.frac=frac0;
    for i=1:(G.cartDims(3)-1)
        state=EOSPP(state);
        frac0=state.frac;
        p00=state.pressure;
        [rho0,mu0]=rhomu_p_frac(p00,frac0,m_NaCl,Tk);
        rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),1)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),1)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),2)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),2)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state=EOSPP(state);
        frac0=state.frac;
        p00=state.pressure;
        [rho0,mu0]=rhomu_p_frac(p00,frac0,m_NaCl,Tk);
        s0=fluid.S(state);
        rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),1)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),1)...
            +0.5.*(rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))+rhoeff0(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2)))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),2)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),2)...
            +0.5.*(rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))+rhoeff0(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2)))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);        
    end
    state=EOSPP(state);
    frac0=state.frac;
    p00=state.pressure;
    s0=fluid.S(state);
    [rho0,mu0]=rhomu_p_frac(p00,frac0,m_NaCl,Tk);
end

if (p0(:,2)>p0(:,1)) && property==2 % two phases; compressible, insoluble
    state.pressure=p00;
    s0=fluid.S(state);
    state.m_NaCl=m_NaCl;
    state.frac=frac0;
    for i=1:(G.cartDims(3)-1)   
        p00=state.pressure;
        [rho0,mu0]=rhomu_p_frac(p00,frac0,m_NaCl,Tk);
        rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),1)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),1)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),2)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),2)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        
        p00=state.pressure;
        [rho0,mu0]=rhomu_p_frac(p00,frac0,m_NaCl,Tk);
        s0=fluid.S(state);
        rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),1)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),1)...
            +0.5.*(rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))+rhoeff0(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2)))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),2)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),2)...
            +0.5.*(rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))+rhoeff0(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2)))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
    end
    p00=state.pressure;
    s0=fluid.S(state);
    [rho0,mu0]=rhomu_p_frac(p00,frac0,m_NaCl,Tk);
end


if (p0(:,2)>p0(:,1)) && property==1 % two phases; incompressible, insoluble
    state.pressure=p00;
    s0=fluid.S(state);
    state.m_NaCl=m_NaCl;
    state.frac=frac0;
    rho0=rhoref0;
    mu0=muref0;
    %   [rho0,mu0]=rhomu_p_frac(p00,frac0,m_NaCl);
    for i=1:(G.cartDims(3)-1) 
        p00=state.pressure;
        rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),1)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),1)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),2)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),2)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);  
    end
    p00=state.pressure;
    s0=fluid.S(state);
end

if (p0(:,2)<=p0(:,1))  && property==4% one phase
    state.pressure=p00;
    s0=fluid.S(state);
    state.m_NaCl=m_NaCl;
    state=rock_frac(state);
    state.poro=rock.poro;
    state.poro0=rock.poro;
    for i=1:(G.cartDims(3)-1)
          state=rock_frac(state);
        frac0=state.frac;
        p00=state.pressure;
        [rho0,mu0]=rhomu_p_frac_rock(p00,frac0,m_NaCl);
        rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),1)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),1)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),2)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),2)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state=rock_frac(state);
        frac0=state.frac;
        p00=state.pressure;
        [rho0,mu0]=rhomu_p_frac_rock(p00,frac0,m_NaCl);
        rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),1)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),1)...
            +0.5.*(rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))+rhoeff0(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2)))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),2)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),2)...
            +0.5.*(rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))+rhoeff0(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2)))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
    end
    state=rock_frac(state);
    frac0=state.frac;
    p00=state.pressure;
    [rho0,mu0]=rhomu_p_frac_rock(p00,frac0,m_NaCl);
end


if (p0(:,2)<=p0(:,1))  && property==3% one phase
    state.pressure=p00;
    s0=fluid.S(state);
    state.m_NaCl=m_NaCl;
    state=EOSPP(state);
    for i=1:(G.cartDims(3)-1)
        state=EOSPP(state);
        frac0=state.frac;
        p00=state.pressure;
        [rho0,mu0]=rhomu_p_frac(p00,frac0,m_NaCl,Tk);
        rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),1)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),1)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),2)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),2)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state=EOSPP(state);
        frac0=state.frac;
        p00=state.pressure;
        [rho0,mu0]=rhomu_p_frac(p00,frac0,m_NaCl,Tk);
        rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),1)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),1)...
            +0.5.*(rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))+rhoeff0(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2)))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),2)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),2)...
            +0.5.*(rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))+rhoeff0(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2)))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
    end
    state=EOSPP(state);
    frac0=state.frac;
    p00=state.pressure;
    [rho0,mu0]=rhomu_p_frac(p00,frac0,m_NaCl,Tk);
end

if (p0(:,2)<=p0(:,1))  && property==2% one phase compressible, insoluble
    state.pressure=p00;
    s0=fluid.S(state);
    state.m_NaCl=m_NaCl;
    state.frac=frac0;
    for i=1:(G.cartDims(3)-1) 
        p00=state.pressure;
       
        [rho0,mu0]=rhomu_p_frac(p00,frac0,m_NaCl,Tk);
        rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),1)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),1)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),2)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),2)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        
        frac0=state.frac;
        p00=state.pressure;
        [rho0,mu0]=rhomu_p_frac(p00,frac0,m_NaCl,Tk);
        rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),1)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),1)...
            +0.5.*(rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))+rhoeff0(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2)))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),2)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),2)...
            +0.5.*(rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))+rhoeff0(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2)))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
    end
    frac0=state.frac;
    p00=state.pressure;
    [rho0,mu0]=rhomu_p_frac(p00,frac0,m_NaCl,Tk);   
end
if (p0(:,2)<=p0(:,1))  && property==1% one phase incompressible, insoluble
    state.pressure=p00;
    s0=fluid.S(state);
    state.m_NaCl=m_NaCl;
    state.frac=frac0;
    rho0=rhoref0;
    mu0=muref0;
    for i=1:(G.cartDims(3)-1) 
        p00=state.pressure;
        rhoeff0=effective_rho(rho0,mu0,s0,fluid);
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),1)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),1)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);
        
        state.pressure(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),2)=...
            state.pressure((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),2)...
            +rhoeff0((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2))...
            .* sum((G.cells.centroids(i*G.cartDims(1)*G.cartDims(2)+1:(i+1)*G.cartDims(1)*G.cartDims(2),:)-G.cells.centroids((i-1)*G.cartDims(1)*G.cartDims(2)+1:(i)*G.cartDims(1)*G.cartDims(2),:)).*gravity(),2);  
    end
    frac0=state.frac;
    p00=state.pressure;
 end

resSol = struct('pressure', state.pressure,             ...
    'flux',     zeros([nf, 1]), ...
    's',        s0,...
    'rho',rho0,...
    'mu', mu0,...
    'frac0',frac0,...
    'rhoref',rhoref,...
    'dp',[],...
    'ds',[],...
    'pressure0',state.pressure,...
    'mu0',mu0,...
    'rho0',rho0,...
    'drho0',[],...
    'dmu0', [],...
    's0',s0,...
    'frac',frac0,...
    'dfrac',zeros(size(s0,1),4),...
    'Tk',Tk,...
    'm_NaCl',m_NaCl);
end



function rhoeff = effective_rho(rho0,mu0,s0, fluid)
kr        = fluid.relperm(s0);
mob    = bsxfun(@rdivide, kr, mu0);
totmob = sum(mob, 2);
rhoeff = sum(bsxfun(@times, mob, rho0(:,1:2)), 2) ./ totmob;
end
