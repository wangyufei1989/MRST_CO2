function fluid = initFluidPPVG(varargin)

% define the relative permeability and retention curves based on Van Genuchten model.
%
% SYNOPSIS:
% fluid = initFluidPPVG(varargin)
%
% DESCRIPTION:
%   This function serves to define the .relative permeability and retention curves
%   based on Van Genuchten model.
%Description is already included in the user's guide.

%
% REQUIRED PARAMETERS:
%    'sr', residual saturation, 'kwm', maximum relative permeability,
%'pc_scale', scaling parameter, 'alpham', shape factor ;


% RETURNS:
%   fluid   - an object containing relative permeability and retention curves based on Van Genuchten model.
%
% SEE ALSO:
%   er2p, er3p,initFluidPPVGM

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


opt = struct( 'sr', [], 'kwm',[], 'pc_scale',[],'alpham',[]);
opt = merge_options(opt, varargin{:});
kr   = @(s,varargin) relperm(s, opt, varargin{:});
S   = @(state) S_funct(state,opt);
fluid = struct(   'relperm'   , kr,...
    'S'        , S);
end
%-------------------------------------------------------------------------
function varargout = S_funct(state,opt)

pc=(state.pressure(:,2)-state.pressure(:,1));
pc(pc<0)=0;
p_d=0;
ps=opt.pc_scale;
ps=ps./10^5*opt.alpham(1);
m=opt.alpham(2); n=1/(1-m);
varargout{1}   =[(1+(ps.*pc).^n).^(-m).*(1-opt.sr(1)-opt.sr(2))+opt.sr(1)-ps.*p_d.*pc,1-((1+(ps.*pc).^n).^(-m).*(1-opt.sr(1)-opt.sr(2))+opt.sr(1)-ps.*p_d.*pc)];

if nargout > 1, varargout{2} = -m.*(1+(ps.*pc).^n).^(-m-1).*ps.^n.*n.*pc.^(n-1).*(1-opt.sr(1)-opt.sr(2))-ps.*p_d;
end
end

%--------------------------------------------------------------------------

function varargout = relperm(s, opt, varargin)
[s1, s2, den] = modified_saturations(s, opt);
s_d=2e-2;
s_d=0;
kwm = opt.kwm;
m=0.8;
varargout{1}    = [ kwm(1) * s1 .^ 0.5.*(1-(1-s1.^(1/m)).^m).^2+s_d.*s(:,1), kwm(2) * s2 .^ (1/2).*(1-s1.^(1/m)).^(2*m)-s_d.*s(:,1)+s_d];
if nargout > 1
    dkr11=zeros(size(s1,1),1);
    dkr22=zeros(size(s1,1),1);
    ins1=(s1>0&s1<1);
    dkr11(ins1)=kwm(1)./den.* 0.5.*s1(ins1).^(-0.5).*(1-(1-s1(ins1).^(1/m)).^m).^2+...
        kwm(1)./den.*s1(ins1).^(0.5).*2.*(1-(1-s1(ins1).^(1/m)).^m).*(-m.*(1-s1(ins1).^(1./m)).^(m-1)).*(-1./m.*s1(ins1).^(1./m-1))+s_d;
    dkr22(ins1)=-kwm(2)./den./2.*(1-s1(ins1)).^(-1/2).*(1-s1(ins1).^(1/m)).^(2*m)...
        +kwm(2)./den.*(1-s1(ins1)).^(1/2).*2.*m.*(1-s1(ins1).^(1/m)).^(2*m-1).*(-1/m).*s1(ins1).^(1/m-1)-s_d;
    dkr11(~ins1)=0;
    dkr22(~ins1)=0;
    ins2=(s1==1);
    dkr11(ins2)=0;
    varargout{2} = [ dkr11,    dkr22];
end

if nargout > 2
    a = n .* (n - 1);
    varargout{3} = [ kwm(1) * a(1) .* s1 .^ (n(1) - 2), ...
        kwm(2) * a(2) .* s2 .^ (n(2) - 2)] ./ den;
end
end


%--------------------------------------------------------------------------

function [s1, s2, den] = modified_saturations(s, opt)
den = 1 - opt.sr(3)-opt.sr(4);
s1  = (    s(:,1) - opt.sr(3)) ./ den;  s1(s1 < 0) = 0;  s1(s1 > 1) = 1;
s2  = (1 - s(:,1) -opt.sr(4)) ./ den;  s2(s2 < 0) = 0;  s2(s2 > 1) = 1;
end
