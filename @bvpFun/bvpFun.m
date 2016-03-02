classdef bvpFun < skpObject
%bvpFun is the base class for the BVP functions.
%
% obj = bvpFun(D, N, phi)
%   Base class for SKPrime BVP based functions. Not a truly functional
%   class.
%
% Provides properties:
%   phiFun
%   trunctation
%
% See also skpObject.

% Everett Kropf, 2015
% 
% This file is part of SKPrime.
% 
% SKPrime is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SKPrime is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SKPrime.  If not, see <http://www.gnu.org/licenses/>.

properties(SetAccess=protected)
    phiFun                      % Unknown function object
    truncation = 64             % Series truncation setting
end

methods
    function bvp = bvpFun(D, N, phi)
        if ~nargin
            sarg = {};
        elseif isa(D, 'bvpFun')
            phi = D.phiFun;
            if nargin < 2
                N = D.truncation;
            end
            sarg = {D.domain};
        else
            sarg = {D};
            if nargin < 2
                N = [];
            end
            if nargin < 3
                phi = [];
            end
        end
        bvp = bvp@skpObject(sarg{:});
        if ~nargin
            return
        end
        
        if ~isempty(N)
            if ~isnumeric(N) || N <= 0 || N ~= floor(N)
                error('SKPrime:invalidArgument', ...
                    'N must be a non-zero, positive integer.')
            end
            bvp.truncation = N;
        end
        
        if isempty(phi)
            bvp.phiFun = schwarz(bvp.domain, bvp.truncation);
        else
            if ~isa(phi, 'schwarz')
                error('SKPrime:invalidArgument', ...
                    '"phi" must be a "schwarz" object.')
            end
            bvp.phiFun = phi;
        end
    end
    
    function dw = diff(bvp)
        %gives derivative of BVP function via DFT and Cauchy continuation.
        %
        %   dw = dftDeriv(bvp)
        %      Returns a function handle dw to the derivative of the bvpFun
        %      object with repsect to the complex variable via
        %      bvpFun.feval(). The derivative is restricted to the unit disk.
        
        dw = dftDerivative(bvp, @bvp.feval);
    end
end

methods(Access=protected)
    function dw = dftDerivative(bvp, F)
        %gives derivative via DFT and continuation.
        %
        %   dw = dftDeriv(bvp, F)
        %      Returns the derivative of function handle F using the DFT
        %      and Cauchy continuation. The derivative is restricted to the
        %      unit disk.
        
        nf = 256;
        [d, q, m] = domainDataB(bvp.domain);
        zf = bsxfun(@plus, d.', ...
            bsxfun(@times, q', exp(2i*pi/nf*(0:nf-1)')));
        dmult = [0:nf/2-1, 0, -nf/2+1:-1]';
        ipos = nf/2:-1:1;
        ineg = nf/2+1:nf;

        dk = complex(nan(nf, m+1));
        p = cell(1, m+1);
        for j = 1:m+1
            dk(:,j) = 1i*dmult.*fft(F(zf(:,j)))/nf;
            p{j} = @(z) polyval(dk(ipos,j), z) ...
                + polyval([dk(ineg,j); 0], 1./z);
        end
        
        function [val, onBdry] = dwBdry(z)
            val = complex(nan(size(z)));
            if nargout > 1
                onBdry = false(size(z));
            end
            for i = 1:m+1
                onCj = abs(q(i) - abs(z - d(i))) < eps(2);
                if any(onCj(:))
                    eij = (z(onCj) - d(i))/q(i);
                    val(onCj) = p{i}(eij)./(1i*q(i)*eij);
                    if nargout > 1
                        onBdry = onBdry | onCj;
                    end
                end
            end
        end
        
        dwCont = SKP.bmcCauchy(@dwBdry, bvp.domain, bvp.truncation);
        
        function val = dwEval(z)
            [val, onBdry] = dwBdry(z);
            val(~onBdry) = dwCont(z(~onBdry));
        end
        
        dw = @dwEval;
    end
end

end
