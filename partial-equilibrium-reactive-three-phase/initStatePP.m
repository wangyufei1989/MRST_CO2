function state = initStatePP(G, W, p0, frac,rhoref,muref,fluid,varargin)
% initialize the system.
%
% SYNOPSIS:
% state = initStatePP(G, W, p0, frac,rhoref,muref,fluid,varargin)
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
%   state   - an object containing the initialized variables. 
%
% SEE ALSO:
%   initReSolPP

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
y.pressure=p0;
s0=fluid.S(y);
   state = initResSolPP(G, p0,s0, frac,rhoref,muref,fluid,varargin{:});

   if ~isempty(W)
      state.wellSol = initWellSol(W, p0);
   end
end
