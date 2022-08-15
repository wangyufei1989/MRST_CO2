function varargout = ticif(bool)
%Evaluate function TIC if input is true.
%
% SYNOPSIS:
%            ticif(bool)
%   tstart = ticif(bool)
%
% PARAMETERS:
%   bool - Boolean variable (true/false).
%
% RETURNS:
%   tstart - Save start time to output argument.  The numeric value of
%            'tstart' is only useful as an input argument for a subsequent
%            call to functions 'toc' or 'tocif'.
%
%            If 'bool' is FALSE, then 'tstart' is the emtpy array ([]).
%
% COMMENTS:
%   Function used for making code cleaner where verbose option is used.
%
% SEE ALSO:
%   tic, tocif, dispif.

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


if bool,
   if nargout > 0,
      varargout{1} = tic;
   else
      tic;
   end
elseif nargout > 0,
   varargout{1} = [];
end
end
