function [k,dk]=kdk(rock,state)
 %  Közeny-Carman model for iintrinsic permeability
%
% SYNOPSIS:
%  [k,dk]=kdk(rock,state)
%
% DESCRIPTION:
%   This function serves to calculate the intrinsic permeaility
% based on the initial porosity, initial permeability and new porosity. 
%Description is already included in the user's guide. 

%
% REQUIRED PARAMETERS:
% rock       -an object containing initial permeability and porosity
%state      -  state variables containing new porosity etc.



% RETURNS:
%   k   -  new intrinsic permeability
% dk   -dk/dphi

% SEE ALSO:
% er3p

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


poro=state.poro;
poro0=rock.poro;
k0=rock.perm;
gamma=10;
k=k0.*poro.^gamma./(1-poro).^2.*(1-poro0).^2./poro0.^gamma;
dk=k0.*(gamma.*poro.^(gamma-1).*(1-poro).^2-2.*(poro-1).*poro.^gamma)./(1-poro).^(4).*(1-poro0).^2./poro0.^gamma;
end