function fluid = initCoreyFluidPP_3p_benchmark(varargin)
%Initialize incompressible two-phase fluid model with capillary forces
%
% SYNOPSIS:
%   fluid = initSimpleFluidPc('pn1', pv1, ...)
%
% PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining specific fluid
%             characteristics.  The following parameters must be defined
%             with one value for each of the two fluid phases:
%               - mu  -- Phase viscosities in units of Pa*s.
%               - rho -- Phase densities in units of kilogram/meter^3.
%               - n   -- Phase relative permeability exponents.
%               - pc_scale -- Constant multiplying the linear capillary
%                        term, i.e., pc(S) = pc_scale*(1-S)
%
% RETURNS:
%   fluid - Fluid data structure as described in 'fluid_structure'
%           representing the current state of the fluids within the
%           reservoir model.
%
% EXAMPLE:
%   fluid = initSimpleFluidPc('mu' , [   1,  10]*centi*poise     , ...
%                             'rho', [1014, 859]*kilogram/meter^3, ...
%                             'n'  , [   2,   2], ...
%                             'pc_scale', 2*barsa);
%
%   s = linspace(0, 1, 101).'; kr = fluid.relperm(s);
%   subplot(1,2,1), plot(s, kr), legend('kr_1(S)', 'kr_2(S)')
%   x.s = [s 1-s]; pc = fluid.pc(x);
%   subplot(1,2,2), plot(s, pc); legend('P_c(S)');
%
% SEE ALSO:
%   fluid_structure, initSimpleFluid, solveIncompFlow.

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
ps=opt.pc_scale;
ps=ps./10^5.*opt.alpham(1);
m=opt.alpham(2); n=1/(1-m);
  varargout{1}   =[(1+(ps.*pc).^n).^(-m).*(1-opt.sr(1))+opt.sr(1),1-((1+(ps.*pc).^n).^(-m).*(1-opt.sr(1))+opt.sr(1))];
 %varargout{1}   ( pc>5e5,1) =(1+(ps(pc>5e5).*5e5).^n).^(-m).*(1-sum(opt.sr))+opt.sr(1); 
 %varargout{1}   ( pc>5e5,2) =1-(1+(ps(pc>5e5).*5e5).^n).^(-m).*(1-sum(opt.sr))-opt.sr(1);
%varargout{1}(varargout{1}   ( :,1)<2.08e-4,1) =0; 
%varargout{1}(varargout{1}   ( :,1)<2.08e-4,2) =1; 
  
    if nargout > 1, varargout{2} = -m.*(1+(ps.*pc).^n).^(-m-1).*ps.^n.*n.*pc.^(n-1).*(1-opt.sr(1)); 
    varargout{2}   ( pc==0) =-2e-10;
    %varargout{2}(varargout{1}   ( :,1)<2.08e-4) =0; 

    end
    
    
    
   
end



% function varargout = S_funct(state,opt)
% 
% pc=(state.pressure(:,2)-state.pressure(:,1));
% pc(pc<0)=0; pc(pc>1e8)=1e8;
% % ps=opt.pc_scale;
% % ps=ps./10^5.*0.37*10;
% % n=4.37;m=0.77;
% %   varargout{1}   =[(1+(ps.*pc).^n).^(-m).*(1-opt.sr(1))+opt.sr(1),1-((1+(ps.*pc).^n).^(-m).*(1-opt.sr(1))+opt.sr(1))];
%  %varargout{1}   ( pc>5e5,1) =(1+(ps(pc>5e5).*5e5).^n).^(-m).*(1-sum(opt.sr))+opt.sr(1); 
%  %varargout{1}   ( pc>5e5,2) =1-(1+(ps(pc>5e5).*5e5).^n).^(-m).*(1-sum(opt.sr))-opt.sr(1);
% %varargout{1}(varargout{1}   ( :,1)<2.08e-4,1) =0; 
% %varargout{1}(varargout{1}   ( :,1)<2.08e-4,2) =1; 
%   varargout{1}   =[1-pc.*1e-8,pc.*1e-8];
%     if nargout > 1, varargout{2} = -1e-8.*ones(size(pc)); 
%     %varargout{2}   ( pc>5e5) =0;
%     %varargout{2}(varargout{1}   ( :,1)<2.08e-4) =0; 
% 
% varargout{2}(pc<0)=0;varargout{2}(pc>1e8)=0;
%     end
%     
%     
%     
%    
% end

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------

function varargout = relperm(s, opt, varargin)
   [s1, s2, den] = modified_saturations(s, opt);
   kwm = opt.kwm;
   m=opt.alpham(2);
   varargout{1}    = [ kwm(1) * s1 .^ 0.5.*(1-(1-s1.^(1/m)).^m).^2, kwm(2) * s2 .^ (1/2).*(1-s1.^(1/m)).^(2*m)];
   if nargout > 1
       dkr11=zeros(size(s1,1),1);
        dkr22=zeros(size(s1,1),1);
     ins1=(s1>0&s1<1);
      %ins2=(s2>0&s2<1);
      dkr11(ins1)=kwm(1)./den.* 0.5.*s1(ins1).^(-0.5).*(1-(1-s1(ins1).^(1/m)).^m).^2+...
          kwm(1)./den.*s1(ins1).^(0.5).*2.*(1-(1-s1(ins1).^(1/m)).^m).*(-m.*(1-s1(ins1).^(1./m)).^(m-1)).*(-1./m.*s1(ins1).^(1./m-1));
      dkr22(ins1)=-kwm(2)./den./2.*(1-s1(ins1)).^(-1/2).*(1-s1(ins1).^(1/m)).^(2*m)...
         +kwm(2)./den.*(1-s1(ins1)).^(1/2).*2.*m.*(1-s1(ins1).^(1/m)).^(2*m-1).*(-1/m).*s1(ins1).^(1/m-1);
       dkr11(~ins1)=0;
       %dkr11(s1==1)=1e-15;
      dkr22(~ins1)=0;
      % dkr22(s1==1)=-1e-15;
      varargout{2} = [ dkr11,    dkr22];
   end

   if nargout > 2,
      a = n .* (n - 1);
      varargout{3} = [ kwm(1) * a(1) .* s1 .^ (n(1) - 2), ...
                       kwm(2) * a(2) .* s2 .^ (n(2) - 2)] ./ den;
   end
end





%--------------------------------------------------------------------------

function [s1, s2, den] = modified_saturations(s, opt)
   den = 1 - opt.sr(1);
   s1  = (    s(:,1) - opt.sr(1)) ./ den;  s1(s1 < 0) = 0;  s1(s1 > 1) = 1;
   s2  = (1 - s(:,1) ) ./ den;  s2(s2 < 0) = 0;  s2(s2 > 1) = 1;
end
