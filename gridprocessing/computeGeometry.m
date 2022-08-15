function G = computeGeometry(G, varargin)
%Compute geometry of grid.
%
% SYNOPSIS:
%   G = computeGeometry(G)
%   G = computeGeometry(G, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure.
%
%   'pn'/pv - List of 'key'/value pairs for supplying optional parameters.
%             The supported options are
%               - verbose --
%                   Whether or not to display informational messages during
%                   the computational process.
%                   Logical.  Default value: Verbose = false (don't display
%                   any informational messages).
%
%               - hingenodes --
%                   Structure with fields 'faces' and 'nodes'.  A hinge
%                   node is an extra center node for a face, that is used
%                   to triangulate the face geometry.  For each face number
%                   F in 'faces' there is a row in 'nodes' which holds the
%                   node coordinate for the hinge node belonging to face F.
%
%                   Default value: hingenodes = [] (no additional center
%                   nodes).
%
%               - MaxBlockSize --
%                   Maximum number of grid cells to process in a single
%                   pass.  Increasing this number may reduce overall
%                   computational time, but will increase total memory use.
%                   If empty (i.e., if ISEMPTY(MaxBlockSize) is TRUE) or
%                   negative, process all grid cells in a single pass.
%
%                   Numeric scalar.  Default value: MaxBlockSize = 20e3.
%
% RETURNS:
%   G - Grid structure with added fields:
%         - cells
%             - volumes   -- A G.cells.num-by-1 array of cell volumes.
%
%             - centroids -- A G.cells.num-by-SIZE(G.nodes.coords, 2) array
%                            of (approximate) cell centroids.
%
%         - faces
%             - areas     -- A G.faces.num-by-1 array of face areas.
%
%             - normals   -- A G.faces.num-by-G.griddim array
%                            of face normals.
%
%             - centroids -- A G.faces.num-by-SIZE(G.nodes.coords, 2) array
%                            of (approximate) face centroids.
%
% COMMENTS:
%   Individual face normals have length (i.e., Euclidian norm) equal to
%   the corresponding face areas.  In other words, subject to numerical
%   round-off, the identity
%
%         NORM(G.faces.normals(i,:), 2) == G.faces.areas(i)
%
%   holds for all faces i=1:G.faces.num .
%
%   In three space dimensions, i.e., when G.griddim == 3,
%   function 'computeGeometry' assumes that the nodes on a given face, f,
%   are ordered such that the face normal on f is directed from cell
%   G.faces.neighbors(f,1) to cell G.faces.neighbors(f,2).
%
% SEE ALSO:
%   grid_structure.

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

% Setup
assert(size(G.faces.nodes, 2)==1);
opt     = struct('verbose',              mrstVerbose, ...
                 'findNeighbors',        false,       ...
                 'hingenodes',           [],          ...
                 'MaxBlockSize',         20e3);
opt     = merge_options(opt, varargin{:});

assert(isempty(opt.hingenodes) || G.griddim == 3, ...
   'Hinge nodes are only supported for 3D grids.');

% Possibly find neighbors
if opt.findNeighbors,
   G.faces.neighbors = findNeighbors(G);
   G = findNormalDirections(G);
else
   if ~isfield(G.faces, 'neighbors'),
      warning(msgid('GridType:Incomplete'), ...
         ['No field ''faces.neighbors'' found. ',...
         'Adding plausible values... proceed with caution!']);
      G.faces.neighbors = findNeighbors(G);
      G = findNormalDirections(G);
   end
end

% Main part
if G.griddim == 3,

   [faceAreas, faceNormals, faceCentroids, ...
      cellVolumes, cellCentroids] = geom_3d(G, opt);

elseif (G.griddim ==2) && (size(G.nodes.coords,2)==2)

   [faceAreas, faceNormals, faceCentroids, ...
      cellVolumes, cellCentroids] = geom_2d2(G, opt);

elseif (G.griddim ==2) && (size(G.nodes.coords,2)==3)

   [faceAreas, faceNormals, faceCentroids, ...
      cellVolumes, cellCentroids] = geom_2d3(G, opt);

else
  assert(false);
end

% Update grid
G.faces.areas     = faceAreas;
G.faces.normals   = faceNormals;
G.faces.centroids = faceCentroids;

G.cells.volumes   = cellVolumes;
G.cells.centroids = cellCentroids;

if ~isfield(G, 'type'),
   warning(msgid('GridType:Unknown'),                            ...
          ['Input grid has no known type. ',                     ...
           'I''ll assume it arose from the primordial soup...']);
   G.type = { 'Primordial Soup' };
end

G.type = [G.type, { mfilename }];
end

%--------------------------------------------------------------------------

function [faceAreas, faceNormals, faceCentroids, ...
      cellVolumes, cellCentroids] = geom_3d(G, opt)

   assert (size(G.nodes.coords, 2) == 3, ...
           'Internal error: 3D geometry on non-3D coordinates');

   geom = allocate_geometry_3d(G);

   if opt.verbose,
      nDigits = 1 + floor(log10(G.cells.num));

      t0 = tic;
   end

   for block = partition_cells(G, opt),
      if opt.verbose,
         fprintf('Processing Cells %*d : %*d of %*d ... ', ...
                 nDigits, block(1), ...
                 nDigits, block(2), ...
                 nDigits, G.cells.num);

         t1 = tic;
      end

      geom = geom_3d_cell_block(geom, G, block, opt);

      if opt.verbose,
         t1   = toc(t1);
         rate = (1 + diff(block)) / t1;

         fprintf('done (%.2f [s], %.2e cells/second)\n', t1, rate);
      end
   end

   if opt.verbose,
      fprintf('Total 3D Geometry Processing Time = %.3f [s]\n', toc(t0));
   end

   faceAreas     = geom.face.Areas;
   faceCentroids = geom.face.Centroids;
   faceNormals   = geom.face.Normals;

   cellVolumes   = geom.cell.Volumes;
   cellCentroids = geom.cell.Centroids;
end

%--------------------------------------------------------------------------

function [faceAreas, faceNormals, faceCentroids, ...
      cellVolumes, cellCentroids] = geom_2d2(G, opt)

   quadArea = @(a, b) abs(a(:,1).*b(:,2) - a(:,2).*b(:,1));

   [faceAreas, faceNormals, faceCentroids, ...
      cellVolumes, cellCentroids] = geom_2d_impl(G, opt, quadArea);
end

%--------------------------------------------------------------------------

function [faceAreas, faceNormals, faceCentroids, ...
      cellVolumes, cellCentroids] = geom_2d3(G, opt)

   dispif(opt.verbose, ...
      'Experimental implementation only available for surface grids\n');

   quadArea = @(a, b) sqrt(sum(cross(a, b) .^ 2, 2));

   [faceAreas, faceNormals, faceCentroids, ...
      cellVolumes, cellCentroids] = geom_2d_impl(G, opt, quadArea);
end

%--------------------------------------------------------------------------

function geom = geom_3d_cell_block(geom, G, block, opt)
   cells = block(1) : block(2);

   [cf, cft] = cell_face_triangle_connectivity(G, cells);

   % Geometric primitives for those faces that aren't already valid.
   geom = ensureValidFaceGeometry(geom, G, cf(:,2), opt);

   % Divide each cell into sub-tetrahedra according to sub-triangles
   % geom.sub, all having one node as
   %
   %     cCenter = sum(faceCentroids) / #faceCentroids

   cCenter = sparse(double(cf(:,1)), 1 : size(cf, 1), 1) ...
      * [ geom.face.Centroids(cf(:,2), :), ones([size(cf, 1), 1]) ];

   cCenter = bsxfun(@rdivide, cCenter(:, 1 : end - 1), cCenter(:, end));

   relSubC = geom.sub.Centroids(cft(:,3), :) - cCenter(cft(:,1), :);

   % Outward normal on cell-face triangles with respect to cell.
   cfsign = 2*double(G.faces.neighbors(cft(:,2), 1) == ...
                     reshape(cells(cft(:,1)), [], 1)) - 1;

   outNormals = bsxfun(@times, geom.sub.Normals(cft(:,3), :), ...
                       geom.sub.NormalSigns(cft(:,3)) .* cfsign);

   tVolumes   = (1/3) * sum(relSubC .* outNormals, 2);
   tCentroids = (3/4) * relSubC;

   relCentroid = sparse(double(cft(:,1)), 1 : size(cft, 1), tVolumes) ...
      * [ tCentroids, ones([size(cft, 1), 1]) ];

   volume      = relCentroid(:, end);
   relCentroid = bsxfun(@rdivide, relCentroid(:, 1 : end - 1), volume);
   centroid    = relCentroid + cCenter;

   geom.cell.Volumes(cells)     = volume;
   geom.cell.Centroids(cells,:) = centroid;
end

%--------------------------------------------------------------------------

function geom = allocate_geometry_3d(G)
   dim = 3;

   [cell.Volumes, cell.Centroids] = ...
      preallocate(G.cells.num, [1, dim]);

   [face.Areas, face.Normals, face.Centroids] = ...
      preallocate(G.faces.num, [1, dim, dim]);

   face.valid = false(size(face.Areas));

   [sub.Centroids, sub.Normals, sub.NormalSigns] = ...
      preallocate(size(G.faces.nodes, 1), [dim, dim, 1]);

   geom = struct('cell', cell, 'face', face, 'sub', sub);
end

%--------------------------------------------------------------------------

function id_blocks = partition_cells(G, opt)
   B = opt.MaxBlockSize;
   M = G.cells.num;

   if isempty(B) || (B < 0),
      B = M;
   end

   B = min(B, floor(M / ceil(M / B)));
   L = floor(M / B);
   R = mod(M, B);

   p = cumsum([ 1, repmat(B, [1, L]) + ...
                   [ones([1, R]), zeros([1, L - R])] ]);

   id_blocks = [ p(1 : (end - 1)) ; p(2 : end) - 1 ];
end

%--------------------------------------------------------------------------

function [cf, cft] = cell_face_triangle_connectivity(G, cells)
   [fp1, fp2] = deal(G.cells.facePos(cells + 0), ...
                     G.cells.facePos(cells + 1));

   cf = [rldecode(1 : numel(fp1), fp2 - fp1, 2) .', ...
         G.cells.faces(mcolon(fp1, fp2 - 1), 1)];

   [np1, np2] = deal(G.faces.nodePos(cf(:,2) + 0), ...
                     G.faces.nodePos(cf(:,2) + 1));

   % One triangle associated with each node of each cell-face.
   cft = [ rldecode(cf, np2 - np1), mcolon(np1, np2 - 1).' ];
end

%--------------------------------------------------------------------------

function [faceAreas, faceNormals, faceCentroids, ...
      cellVolumes, cellCentroids] = geom_2d_impl(G, opt, quadArea)

   [edges, faceAreas, faceNormals, faceCentroids] = face_geom2d(G, opt);

   dispif(opt.verbose, 'Computing cell volumes and centroids...\t\t');
   t0 = ticif (opt.verbose);

   numfaces = diff(G.cells.facePos);

   [cCenter, cellno] = ...
      averageCoordinates(numfaces, faceCentroids(G.cells.faces(:,1), :));

   subArea     = face_geom2d_subarea(G, edges, cCenter, cellno, quadArea);

   subCentroid = (cCenter(cellno, :) + ...
                  2 * faceCentroids(G.cells.faces(:,1), :)) / 3;

   [cellCentroids, cellVolumes, cellVolumes] = ...
      averageCoordinates(numfaces, subCentroid, subArea);       %#ok<ASGLU>

   tocif(opt.verbose, t0)
end

%--------------------------------------------------------------------------

function geom = ensureValidFaceGeometry(geom, G, cfaces, opt)

   ufaces = unique_faces(G.faces.num, cfaces);
   inval  = ~ geom.face.valid(ufaces);

   if any(inval),
      ifaces = ufaces(inval);

      [fA, fN, fC, sC, sN, sNSigns] = face_geom3d(G, ifaces, opt);

      geom.face.Areas    (ifaces)    = fA;
      geom.face.Centroids(ifaces, :) = fC;
      geom.face.Normals  (ifaces, :) = fN;

      six = mcolon(G.faces.nodePos(ifaces + 0), ...
                   G.faces.nodePos(ifaces + 1) - 1);

      geom.sub.Centroids  (six, :) = sC;
      geom.sub.Normals    (six, :) = sN;
      geom.sub.NormalSigns(six)    = sNSigns;

      geom.face.valid(ifaces) = true;
   end
end

%--------------------------------------------------------------------------

function ufaces = unique_faces(numfaces, cfaces)
   present         = false([numfaces, 1]);
   present(cfaces) = true;

   ufaces = find(present);
end

%--------------------------------------------------------------------------

function [fA, fN, fC, sC, sN, sNSigns] = face_geom3d(G, faces, opt)

   % Divide each face into sub-triangles all having one node as
   %
   %   pCenter = sum(node coordinates, 1) / #nodes
   %
   % Compute area-weighted normals, and add to obtain approximate
   % face-normals.  Compute resulting areas and centroids.

   [nodePos, faceNodes] = ...
      indexSubSet(G.faces.nodePos, G.faces.nodes, faces);

   numNodes = diff(nodePos);

   [pCenters, faceNo, Accum, Accum] = ...
      averageCoordinates(numNodes, G.nodes.coords(faceNodes, :));      %#ok

   % Use hinge nodes for selected faces if present.
   if ~isempty(opt.hingenodes) && isstruct(opt.hingenodes) && ...
         all(isfield(opt.hingenodes, { 'faces', 'nodes' })),

      pCenters = insertHingeNodes(pCenters, G, faces, opt.hingenodes);
   end

   [sN, sA, sC] = ...
      face_geom3d_subface(G, nodePos, faceNodes, pCenters(faceNo, :));

   fN      = Accum * sN;
   sNSigns = sign(sum(sN .* fN(faceNo, :), 2));

   [fC, fA, fA] = averageCoordinates(numNodes, sC, sA);         %#ok<ASGLU>

   % Computation above does not make sense for faces with zero area
   i = find(~ (fA > 0));
   if ~ isempty(i),
      warning(msgid('computeGeometry:faceAreas'), ...
             ['%d faces with non-positive area detected.\n', ...
              'Such faces should be removed before calling %s'], ...
              numel(i), mfilename);

      fC(i,:) = pCenters(i,:);
   end
end

%--------------------------------------------------------------------------

function [edges, faceAreas, faceNormals, faceCentroids] = ...
      face_geom2d(G, opt)

   dispif(opt.verbose, 'Computing normals, areas, and centroids...\t');
   t0 = ticif(opt.verbose);

   edges = reshape(G.faces.nodes, 2, []) .';

   [n1, n2] = deal(G.nodes.coords(edges(:,1), :), ...
                   G.nodes.coords(edges(:,2), :));

   edgeLength    = n2 - n1;
   faceAreas     = sqrt(sum(edgeLength .^ 2, 2));

   faceCentroids = (n1 + n2) ./ 2;

   faceNormals   = [edgeLength(:,2), -edgeLength(:,1)];

   tocif(opt.verbose, t0)
end

%--------------------------------------------------------------------------

function [subNormals, subAreas, subCentroids] = ...
      face_geom3d_subface(G, nodePos, faceNodes, pCenters)

   p = nodePos;

   next                 = (2 : size(faceNodes, 1) + 1) .';
   next(p(2 : end) - 1) = p(1 : end-1);

   a = G.nodes.coords(faceNodes      , :);
   b = G.nodes.coords(faceNodes(next), :);

   subNormals   = cross(b - a, pCenters - a) ./ 2;
   subAreas     = sqrt(sum(subNormals .^ 2, 2));
   subCentroids = (a + b + pCenters) ./ 3;
end

%--------------------------------------------------------------------------

function subArea = ...
      face_geom2d_subarea(G, edges, cCenter, cellno, quadArea)
   cellEdges      = edges(G.cells.faces(:,1), :);
   r              = G.faces.neighbors(G.cells.faces(:,1), 2) == cellno;
   cellEdges(r,:) = cellEdges(r, [2, 1]);

   cc = cCenter(cellno, :);
   a  = G.nodes.coords(cellEdges(:,1), :) - cc;
   b  = G.nodes.coords(cellEdges(:,2), :) - cc;

   subArea = quadArea(a, b) ./ 2;
end

%--------------------------------------------------------------------------'

function varargout = preallocate(m, n)
   nout = nargout;

   assert (any(numel(m) == [1, nout]), ...
           'Row size array must be scalar or match number of items');

   assert (any(numel(n) == [1, nout]), ...
           'Column size array must be scalar or match number of items');

   if numel(m) ~= nout, m = repmat(m, [1, nout]); end
   if numel(n) ~= nout, n = repmat(n, [1, nout]); end

   varargout = arrayfun(@NaN, m, n, 'UniformOutput', false);
end

%--------------------------------------------------------------------------

function N = findNeighbors(G)

   % Internal faces
   cellNo         = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
   [cellfaces, j] = sort(G.cells.faces(:,1));
   cellNo         = cellNo(j);
   hf             = find(cellfaces(1:end-1) == cellfaces(2:end));

   N                       = zeros(G.faces.num, 1);
   N(cellfaces(hf,1),   1) = cellNo(hf);
   N(cellfaces(hf + 1), 2) = cellNo(hf + 1);

   % Boundary faces
   isboundary         = true(numel(cellNo), 1);
   isboundary(hf)     = false;
   isboundary(hf + 1) = false;
   %hf                 = find(isboundary);
   N(cellfaces(isboundary), 1) = cellNo(isboundary);
end

%--------------------------------------------------------------------------

function G = findNormalDirections(G)
   assert (all([G.griddim, size(G.nodes.coords, 2)] == 3), ...
          ['Detecting neighbourship based on normal directions ', ...
           'is only supported in 3D grids.']);

   % Assume convex faces.   Compute average of node coordinates.
   fcenters = ...
      averageCoordinates(diff(G.faces.nodePos), ...
                         G.nodes.coords(G.faces.nodes, :));

   % Assume convex cells.   Compute average of face centre coordinates.
   [ccenters, cellno] = ...
      averageCoordinates(diff(G.cells.facePos), ...
                         fcenters(G.cells.faces(:,1), :));

   % Compute triple product v1 x v2 · v3 of vectors v1 = fc-cc, v2 = n1-fc,
   % and v3 = n2-n1 --- cc and fc being cell centres and face centres, n1
   % and n2 being the first and second node of the face.  Triple product
   % should be positive for half-faces with positive sign.

   n1 = G.nodes.coords(G.faces.nodes(G.faces.nodePos(1:end-1)    ), :);
   n2 = G.nodes.coords(G.faces.nodes(G.faces.nodePos(1:end-1) + 1), :);

   v1 = fcenters(G.cells.faces(:,1), :) - ccenters(cellno, :);
   v2 = n1(G.cells.faces(:,1), :) - fcenters(G.cells.faces(:,1), :);
   v3 = n2(G.cells.faces(:,1), :) - n1(G.cells.faces(:,1), :);

   a   = sum(cross(v1, v2) .* v3, 2);
   sgn = 2*(G.faces.neighbors(G.cells.faces(:,1), 1) == cellno) - 1;

   i   = accumarray(G.cells.faces(:,1), a .* sgn) < 0;
   G.faces.neighbors(i, :) = G.faces.neighbors(i, [2, 1]);
end

%--------------------------------------------------------------------------

function [p, i] = indexSubSet(p, i, s)
   [p1, p2] = deal(p(s), p(s + 1));

   p = cumsum([1 ; p2 - p1]);
   i = i(mcolon(p1, p2 - 1));
end

%--------------------------------------------------------------------------

function pCenters = insertHingeNodes(pCenters, G, faces, hingenodes)
   i  = zeros([G.faces.num, 1]);  i(faces) = 1 : numel(faces);
   ix = i(hingenodes.faces);
   p  = ix > 0;

   pCenters(ix(p), :) = hingenodes.nodes(p, :);
end

%--------------------------------------------------------------------------

function [c, no, w, accum] = averageCoordinates(n, c, w)
   if nargin < 3,
      w = 1;
   end

   no    = rldecode(1 : numel(n), n, 2) .';
   accum = sparse(no, 1 : numel(no), w);

   c = accum * [ c, ones([size(c, 1), 1]) ];
   w = c(:, end);
   c = bsxfun(@rdivide, c(:, 1 : end - 1), w);
end
