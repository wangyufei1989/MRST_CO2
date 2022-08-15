function grdecl = permeabilityConverter(grdecl,varargin)
%
%  add tensor permeability to the grdecl.permtensor
%
% SYNOPSIS:
%  grdecl = permeabilityConverter(grdecl,varargin)
%
% PARAMETERS:
%   grdecl  - struct representing the grdecl file with
%             permeability modified or possibly new extra
%             fields PERMXY PERMZX PERMYZ
%
%  OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%  'K_system' - which define the system the PERM tesor is defiend valid values are
%                   - 'xyz' the coordinate system DEFAULT
%                   - 'cell' the cell system defined by direction between centroids of the logical cell
%                   - 'bedding_plane' the same above for xy direction but
%                   use the cross product of these direction for z
%                   direction
%                   - 'bedding_plane_normal' the same as above but now the
%                   xy system is made orthogonal
%
% RETURNS:
%   grdecl - retrun grdecl with modified  PERM*
%
% COMMENTS:
%
% SEE ALSO:
%   grdecl2rock
%
opt = struct('K_system','xyz');
opt = merge_options(opt, varargin{:});
if(strcmp(opt.K_system,'xyz'))
    disp('Using permeability system xyz nothing to do');
    return;
end


[X,Y,Z]  = buildCornerPtNodes(grdecl);
CCOORD = {X,Y,Z};
PERM =[grdecl.PERMX,grdecl.PERMY,grdecl.PERMZ];
numcells = prod(grdecl.cartDims);
for i = 1:3
    COORD{i}{1} = CCOORD{i}(1:2:end,1:2:end,1:2:end);
    COORD{i}{2} = CCOORD{i}(2:2:end,1:2:end,1:2:end);
    COORD{i}{3} = CCOORD{i}(1:2:end,2:2:end,1:2:end);
    COORD{i}{4} = CCOORD{i}(2:2:end,2:2:end,1:2:end);
    COORD{i}{5} = CCOORD{i}(1:2:end,1:2:end,2:2:end);
    COORD{i}{6} = CCOORD{i}(2:2:end,1:2:end,2:2:end);
    COORD{i}{7} = CCOORD{i}(1:2:end,2:2:end,2:2:end);
    COORD{i}{8} = CCOORD{i}(2:2:end,2:2:end,2:2:end);
end
for i=1:3
    for j=1:8
        COORD{i}{j} = reshape(COORD{i}{j},prod(grdecl.cartDims),1);
    end
end
for i=1:3
    dir{i} = zeros(prod(grdecl.cartDims),3);
end
for i=1:3
    dir{1}(:,i) = sum( [COORD{i}{[2,4,6,8]}]-[COORD{i}{[1,3,5,7]}], 2);
    dir{2}(:,i) = sum( [COORD{i}{[3,4,7,8]}] -[COORD{i}{[1,2,5,6]}], 2);
    dir{3}(:,i) = sum( [COORD{i}{[5:8]}]-[COORD{i}{[1:4]}], 2);
end

if(strcmp(opt.K_system,'cell'))
    disp(['Using cell_xyz system for permeabiltiy']);
elseif(strcmp(opt.K_system,'bedding_plane'))
    dir{3} = cross(dir{1},dir{2});
elseif(strcmp(opt.K_system,'bedding_plane_normal'))
    dir{3} = cross(dir{1},dir{2},2);
    dir{2} = dir{2} - repmat(dot(dir{1},dir{2},2)./sum(dir{1}.^2,2),1,3).*dir{1};
else
    error('Unknown coordinate system for permeability');
end


for i = 1:3;
    dir{i} = dir{i}./repmat(sqrt(sum(dir{i}.*dir{i},2)),1,3);
end



permtensor = zeros(numcells,9);
for i=1:3
    dirtmp= dir{i}';
    dir1 = rldecode(dirtmp(:),repmat(3,3*numcells,1));
    dir2 = rldecode(dir{i},repmat(3,numcells,1));
    dir2 = dir2';
    A=reshape(dir2(:).*dir1(:),9,numcells)';
    %indi=rldecode([1:numcells]',repmat(3,numcells,1));
    %indj=[1:3*numcells];
    %A=sparse(indi,indj,dir{i}');
    %A=A'*A;
    permtensor = permtensor+repmat(PERM(:,i),1,9).*A;
end
ofdiagonal = {'PERMYZ','PERMXY','PERMZX'};
%if(isfield(grdecl,'PERMYZ') || isfield(grdecl,'PERMXY') ||
%isfield(grdecl,'PERMZX'))
if(any(isfield(grdecl,ofdiagonal)))
    disp('Warning: You are using non-diagonal tensor for non orthogonal coordinate system');
    for field = ofdiagonal
        switch field
            case 'PERMXY',
                i=1;
                j=2;
            case 'PERMYZ',
                i=2;
                j=3
            case 'PERMZX',
                i=1;
                j=3;
        end
    end
    dirtmp= dir{i}';
    dir1 = rldecode(dirtmp(:),repmat(3,3*numcells,1));
    dir2 = rldecode(dir{j},repmat(3,numcells,1));
    dir2 = dir2';
    A=reshape(dir2(:).*dir1(:),9,numcells)';
    A=(A(:,1:end)+A(:,end:-1:1))/2;
    permtensor = permtensor+repmat(PERM(:,i),1,9).*A;
end


%grdecl.permtensor = [permtensor(:,1:3),permtensor(:,4:6),pemtensor(:,9)];
grdecl.PERMX = permtensor(:,1); grdecl.PERMXY = permtensor(:,2);grdecl.PERMZX = permtensor(:,3);
grdecl.PERMY = permtensor(:,5); grdecl.PERMYZ = permtensor(:,6);grdecl.PERMZ = permtensor(:,9);