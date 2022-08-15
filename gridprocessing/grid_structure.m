%Grid structure used in the MATLAB Reservoir Simulation Toolbox.
%
% SYNOPSIS:
%   1) Construct Cartesian grid.
%        G = cartGrid(...);
%        G = computeGeometry(G);
%
%   2) Read corner point grid specification from file.
%        grdecl = readGRDECL(...);
%        G      = processGRDECL(grdecl);
%        G      = computeGeometry(G);
%
% RETURNS:
%   G - representation of grids on an unstructured format.  A master
%   structure having the following fields:
%    - cells --
%        A structure specifying properties for each individual cell in the
%        grid.  See CELLS below for details.
%
%    - faces --
%        A structure specifying properties for each individual face in the
%        grid.  See FACES below for details.
%
%    - nodes --
%        A structure specifying properties for each individual node
%        (vertex) in the grid.  See NODES below for details.
%
%    - type --
%        A cell array of strings describing the history of grid constructor
%        and modifier functions through which a particular grid structure
%        has been defined.
%
%    - griddim --
%        The dimension of the grid which in most cases will equal
%        size(G.nodes.coords,2).
%
%   CELLS - Cell structure G.cells:
%    - num --
%        Number of cells in the global grid.
%
%    - facePos --
%        Indirection map of size [num+1,1] into the 'cells.faces' array.
%        Specifically, the face information of cell 'i' is found in the
%        submatrix
%
%            G.cells.faces(facePos(i) : facePos(i+1)-1, :)
%
%        The number of faces of each cell may be computed using the
%        statement DIFF(facePos).
%
%    - faces --
%        A (G.cells.facePos(end)-1)-by-2 array of global faces connected to
%        a given cell.  Specifically, if 'cells.faces(i,1)==j', then global
%        face 'cells.faces(i,2)' is connected to global cell `j'.
%
%        To conserve memory, only the second column is actually stored in
%        the grid structure.  The first column may be re-constructed using
%        the statement
%
%            rldecode(1 : G.cells.num, ...
%                   diff(G.cells.facePos), 2) .'
%
%        A grid constructor may, optionally, append a third column to this
%        array.  In this case 'cells.faces(i,3)' often contains a tag by
%        which the cardinal direction of face 'cells.faces(i,2)' within
%        cell 'cells.faces(i,1)' may be distinguished.  Some ancillary
%        utilities within the toolbox depend on this specific semantics of
%        'cells.faces(i,3)', e.g., to easily specify boundary conditions
%        (functions 'pside' and 'fluxside').
%
%    - indexMap --
%        Maps internal to external grid cells (i.e., active cell numbers to
%        global cell numbers).  In the case of Cartesian grids, indexMap ==
%        (1 : G.cells.num)' .
%
%        For grids with a logically Cartesian topology of dimension 'dims'
%        (a curvilinear grid, a corner-point grid, etc), a map of cell
%        numbers to logical indices may be constructed using the following
%        statement in 2D
%
%             [ij{1:2}] = ind2sub(dims, G.cells.indexMap(:));
%             ij        = [ij{:}];
%
%        and likewise in 3D
%
%             [ijk{1:3}] = ind2sub(dims, G.cells.indexMap(:));
%             ijk        = [ijk{:}];
%
%        In the latter case, ijk(i,:) is global (I,J,K) index of cell 'i'.
%
%    - volumes --
%        A G.cells.num-by-1 array of cell volumes.
%
%    - centroids --
%        A G.cells.num-by-d array of cell centroids in R^d.
%
%   FACES - Face structure G.faces:
%    - num --
%        Number of global faces in grid.
%
%    - nodePos --
%        Indirection map of size [num+1,1] into the 'faces.nodes' array.
%        Specifically, the node information of face 'i' is found in the
%        submatrix
%
%              G.faces.nodes(nodePos(i) : nodePos(i+1)-1, :)
%
%        The number of nodes of each face may be computed using the
%        statement DIFF(nodePos).
%
%    - nodes --
%        A (G.faces.nodePos(end)-1)-by-2 array of vertices in the grid.
%        Specifically, if 'faces.nodes(i,1)==j', then local vertex `i' is
%        part of global face number `j' and corresponds to global vertex
%        'faces.nodes(i,2)'.  To conserve memory, only the last column is
%        stored.  The first column can be constructed using the statement
%
%           rldecode(1:G.faces.num, diff(G.faces.nodePos), 2) .'
%
%    - neighbors --
%        A G.faces.num-by-2 array of neighboring information. Global face
%        number `i' is shared by global cells 'neighbors(i,1)' and
%        'neighbors(i,2)'. One of 'neighbors(i,1)' or 'neighbors(i,2)', but
%        not both, may be zero, meaning that face `i' is an external face
%        shared only by a single cell (the nonzero one).
%
%    - tag --
%        A G.faces.num-by-1 array of face tags.  A tag is a scalar.  The
%        exact semantics of this field is currently undecided and subject
%        to change in future releases of MRST.
%
%    - areas --
%        A G.faces.num-by-1 array of face areas.
%
%    - normals --
%        A G.faces.num-by-d array of *AREA WEIGHTED*, directed face normals
%        in R^d.  The normal on face `i' points from cell
%        'G.faces.neighbors(i,1)' to cell 'G.faces.neighbors(i,2)'.
%
%    - centroids --
%        A G.faces.num-by-d array of face centroids in R^d.
%
%   NODES - Node (vertex) structure G.nodes:
%    - num --
%        Number of global nodes (vertices) in grid.
%
%    - coords --
%        A G.nodes.num-by-d array of physical nodal coordinates in R^d.
%        Global node `i' is at physical coordinate [coords(i,:)]
%
% REMARKS:
%   The grid is constructed according to a right-handed coordinate system
%   where the Z coordinate is interpreted as depth.  Consequently, plotting
%   routines such as plotGrid display the grid with a reversed Z axis.
%
% SEE ALSO:
%   In case the description scrolled off screen too quickly, you may access
%   the information at your own pace using the command
%
%             more on, help grid_structure, more off

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

