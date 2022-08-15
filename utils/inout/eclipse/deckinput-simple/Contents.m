% DECKINPUT-SIMPLE
%   Simplified input of grid-related data keywords in ECLIPSE format
%
% Files
%   applyOperatorSimple - Apply ECLIPSE/FrontSim operator to input array.
%   cutGrdecl           - Extract logically Cartesian subset of a corner-point description
%   readDefaultedKW     - Read data, possibly containing default designators, for a single keyword.
%   readDefaultedRecord - Read data, possibly containing default designators, for a single record.
%   readGRDECL          - Read subset of ECLIPSE GRID file

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

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
