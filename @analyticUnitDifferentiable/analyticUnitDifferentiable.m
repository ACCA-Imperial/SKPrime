classdef(Abstract) analyticUnitDifferentiable
%analyticUnitDifferentiable is the protocol to provide differentiability
%for functions analytic in the unit domain.

% Everett Kropf, 2016
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

properties(Access=protected)
    derivativeCollocationPoints = 128
    derivativeFourierPoints = 256
end

methods(Access=protected)
    function dfun = nthOrderDftDerivative(obj, afun, n)
        %provides the nth order derivative of a given analytic function.
        
        % Catch case that subclass hasn't defined domain property from
        % somewhere.
        try
            if ~isa(obj.domain, 'skpDomain')
                error('SKPrime:runtimeError', ...
                    ['Expected "skpDomain" object, but "%s.domain" is ' ...
                    'a "%s" instead.'], class(obj), class(obj.domain))
            end
        catch err
            switch err.identifier
                case {'MATLAB:noSuchMethodOrField', ...
                        'SKPrime:invalidReference'}
                    error('SKPrime:logicError', ...
                        ['Classes implementing protocol %s must have ' ...
                        'defined the "domain" property.'], ...
                        mfilename('class'))
                    
                otherwise
                    rethrow(err)
            end
        end
        
        dfun = recursiveDftDerivative(obj, afun, n);
    end
end
   
methods(Access=private)
    function dfun = recursiveDftDerivative(obj, afun, n)
        %nth order derivative via recursion.
        
        if n > 1
            afun = recursiveDftDerivative(obj, afun, n - 1);
        end
        dfun = dftDerivative(obj, afun);
    end
    
    function dF = dftDerivative(obj, F)
        %gives derivative via DFT and continuation.
        %
        %  dF = dftDeriv(bvp, F)
        %    Returns the derivative of function handle F using the DFT
        %    and Cauchy continuation for values of F on the boundary of
        %    the domain. The derivative is restricted to the bounded unit
        %    domain. It is assumed that F represents the boundary values of
        %    a function analytic at all points in the bounded unit domain.
        
        % FIXME: Check that obj is of class 'skpObject' so we know it has
        % the domain property.
        
        nf = obj.derivativeFourierPoints;
        D = obj.domain;
        
        [d, q, m] = domainDataB(D);
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
        
        function [val, onBdry] = dfBdry(z)
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
        
        dfCont = SKP.bmcCauchy(@dfBdry, D, obj.derivativeCollocationPoints);
        
        function val = dfEval(z)
            [val, onBdry] = dfBdry(z);
            val(~onBdry) = dfCont(z(~onBdry));
        end
        
        dF = @dfEval;
    end
end

end
