classdef skprime < bvpFun
%SKPRIME is the Shottky-Klein prime function.

% E. Kropf, 2015
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
    parameter    
    vjFuns = {}
end

properties(Access=protected)
    refdata
    logXhatBdry
    logXhatCont
end

properties(Dependent, Access=protected)
    primeCorrect
    normFactor
    logXhatOutBdry
    logXhatOutCont
end

methods
    function skp = skprime(alpha, D, N)
        vjfuns = {};
        phi = [];
        if ~nargin
            sargs = {};
        else
            if iscell(D) && all(cellfun(@(c) isa(c, 'vjCauchy'), D(:)))
                vjfuns = D;
                D = D{1}.domain;
            elseif isa(D, 'skprime')
                vjfuns = D.vjFuns;
                phi = D.phiFun;
                D = D.domain;
            end
            if nargin > 2
                sargs = {D, N, phi};
            else
                sargs = {D, [], phi};
            end
        end
        skp = skp@bvpFun(sargs{:});
        if ~nargin
            return
        end
        
        % Initialize handle data.
        skp.refdata = SKP.primedata;
                
        param = skpParameter(alpha, skp.domain);
        skp.parameter = param;
        
        m = skp.domain.m;
        if m == 0
            return
        end
        
        if ~isempty(param.ison) && param.ison > 0 && param == 0
            dj = skp.domain.centers(param.ison);
            qj = skp.domain.radii(param.ison);
            if abs(abs(dj) - qj) < eps(2)
                warning('SKPrime:unstable', ...
                    ['The case when the parameter is at the origin or\n' ...
                    'at infinity and on a boundary is unstable. Do not\n'...
                    'trust the result of the computation.'])
            end
        end
        
        if isempty(vjfuns)
            skp.vjFuns = cell(m, 1);
            for j = 1:m
                skp.vjFuns{j} = vjCauchy(j, skp);
            end
        elseif numel(vjfuns) ~= m
                error('SKPrime:invalidArgument', ...
                    'Incorrect number of v_j functions given.')
        else
            skp.vjFuns = vjfuns;
        end
        
        skp = solveBVPs(skp);
        
        % Does the prime function need a correction to stay on the proper branch?
        % FIXME: This should check the test point is not in a circle! This should
        % have a maximum retry count on test point selection!
        testpt = 1 - 1e-6;
        while abs(testpt - param) < 1e-2
            testpt = (1 - 1e-6)*exp(2i*pi*rand(1));
        end
        skp.primeCorrect = abs(imag(skp.logXhatCont(testpt) ...
            - conj(skp.logXhatOutCont(1/conj(testpt))))) > pi/4;
    end % ctor
    
    function skp = invParam(skp)
        skp = skpinvparam(skp);
    end
    
    function dw = diff(skp)
        %DIFF Derivative of skprime via DFT and Cauchy continuation.
        %   dw = DIFF(skp), returns a function handle dw to the derivative
        %      of skp with repsect to the complex variable, restricted to
        %      the unit disk.
        
        nf = 256;
        [d, q, m] = domainDataB(skp.domain);
        zf = bsxfun(@plus, d.', ...
            bsxfun(@times, q', exp(2i*pi/nf*(0:nf-1)')));
        dmult = [0:nf/2-1, 0, -nf/2+1:-1]';
        ipos = nf/2:-1:1;
        ineg = nf/2+1:nf;

        dk = complex(nan(nf, m+1));
        p = cell(1, m+1);
        for j = 1:m+1
            dk(:,j) = 1i*dmult.*fft(feval(skp, zf(:,j)))/nf;
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
        
        dwCont = bmcCauchy(@dwBdry, skp.domain, skp.truncation);
        
        function val = dwEval(z)
            [val, onBdry] = dwBdry(z);
            val(~onBdry) = dwCont(z(~onBdry));
        end
        
        dw = @dwEval;
    end
    
    function w = feval(skp, z)
        if skp.domain.m == 0
            w = z - skp.parameter;
            return
        end
        
        [w, L] = evalLogX(skp, z);
        w = exp(w/2);
        
        if skp.primeCorrect
            if abs(skp.parameter) < 1 + eps(2)
                L = ~L;
            end
            w(L) = -w(L);
        end
    end
    
    function X = Xeval(skp, z)
        if skp.domain.m == 0
            X = (z - skp.parameter).^2;
            return
        end
        X = exp(evalLogX(skp, z));
    end
end

methods(Access=protected)
    function skp = bvpInDomain(skp)
        alpha = skp.parameter;
        D = skp.domain;
        N = skp.truncation;

        imLogXhat = primeRHS(alpha, skp.vjFuns);
        phi = solve(skp.phiFun, imLogXhat);
        skp.logXhatBdry = @(z) phi(z) + 1i*imLogXhat(z);
        skp.logXhatCont = bmcCauchy(skp.logXhatBdry, D, N);
        if abs(alpha) < 1 + eps(2)
            skp.normFactor = skp.logXhatCont(alpha);
        end
        
        imLogXhatConj = primeRHS(inv(alpha), skp.vjFuns);
        phiConj = solve(skp.phiFun, imLogXhatConj);
        skp.logXhatOutBdry = @(z) phiConj(z) + 1i*imLogXhatConj(z);
        skp.logXhatOutCont = bmcCauchy(skp.logXhatOutBdry, D, N);
        if abs(alpha) >= 1 + eps(2)
            skp.normFactor = conj(skp.logXhatOutCont(1/conj(alpha)));
        end
    end
    
    function skp = bvpIsUnit(skp)
        alpha = skp.parameter;
        D = skp.domain;
        N = skp.truncation;

        imLogXhat = primeRHS(alpha, skp.vjFuns);
        phi = solve(skp.phiFun, imLogXhat);
        skp.normFactor = phi(alpha) + 1i*imLogXhat(alpha);
        skp.logXhatBdry = @(z) phi(z) + 1i*imLogXhat(z) - skp.normFactor;
        skp.logXhatCont = bmcCauchy(skp.logXhatBdry, D, N);
        
        skp.logXhatOutBdry = skp.logXhatBdry;
        skp.logXhatOutCont = bmcCauchy(skp.logXhatOutBdry, D, N);
    end
    
    function skp = bvpOnInner(skp)
        alpha = skp.parameter;
        D = skp.domain;
        N = skp.truncation;
        [d, q] = domainDataB(D);
        
        imLogXhat = primeRHS(alpha, skp.vjFuns);
        j = alpha.ison;
        vj = skp.vjFuns{j};
        
        phi = solve(skp.phiFun, imLogXhat);
        skp.normFactor = phi(alpha) + 1i*imLogXhat(alpha);
        skp.logXhatBdry = @(z) phi(z) + 1i*imLogXhat(z) - skp.normFactor;
        skp.logXhatCont = bmcCauchy(skp.logXhatBdry, D, N);
        
        if alpha ~= 0
            logConjPart = @(z) 2*log((z - alpha)./(z - 1/conj(alpha))) ...
                + 4i*pi*(real(vj(alpha)) - vj(z)) ...
                - 2*log(q(j+1)/(1 - conj(d(j+1)/alpha)));
            skp.logXhatOutBdry = @(z) logConjPart(z) + skp.logXhatBdry(z);
            skp.logXhatOutCont = @(z) logConjPart(z) + skp.logXhatCont(z);
        else
            imLogXhatConj = primeRHS(inv(alpha), skp.vjFuns);
            phiConj = solve(skp.phiFun, imLogXhatConj);
            skp.logXhatOutBdry = @(z) phiConj(z) + 1i*imLogXhatConj(z);
            skp.logXhatOutCont = bmcCauchy(skp.logXhatOutBdry, D, N);
        end
    end
    
    function skp = bvpOnOuter(skp)
        alpha = skp.parameter;
        D = skp.domain;
        N = skp.truncation;
        [d, q] = domainDataB(D);
        
        imLogXhatConj = primeRHS(inv(alpha), skp.vjFuns);
        j = alpha.ison;
        vj = skp.vjFuns{j};
        
        phi = solve(skp.phiFun, imLogXhatConj);
        skp.normFactor = phi(1/conj(alpha)) + 1i*imLogXhatConj(1/conj(alpha));
        skp.logXhatOutBdry = @(z) phi(z) + 1i*imLogXhatConj(z) ...
            - skp.normFactor;
        skp.logXhatOutCont = bmcCauchy(skp.logXhatOutBdry, D, N);
        
        logConjPart = @(z) 2*log((z - 1/conj(alpha))./(z - alpha)) ...
            + 4i*pi*(real(vj(alpha)) - vj(z)) ...
            - 2*log(q(j+1)/(1 - conj(d(j+1))*alpha));
        skp.logXhatBdry = @(z) logConjPart(z) + skp.logXhatOutBdry(z);
        skp.logXhatCont = @(z) logConjPart(z) + skp.logXhatOutCont(z);
    end
    
    function skps = copyProperties(skp)
        mco = ?skprime;
        for k = 1:numel(mco.PropertyList)
            pname = mco.PropertyList(k).Name;
            skps.(pname) = skp.(pname);
        end
    end
    
    function [logX, inUnit] = evalLogX(skp, z)
        logXhat = complex(nan(size(z)));
        
        inUnit = abs(z) < 1 + eps(2);
        if any(inUnit(:))
            onB = onBoundary(skp, z);
            L = inUnit & onB;
            if any(L(:))
                logXhat(L) = skp.logXhatBdry(z(L));
            end
            L = inUnit & ~onB;
            if any(L(:))
                logXhat(L) = skp.logXhatCont(z(L));
            end
        end
        if any(~inUnit(:))
            onB = onBoundary(skp, 1./conj(z));
            L = ~inUnit & onB;
            if any(L(:))
                logXhat(L) = conj(skp.logXhatOutBdry(1./conj(z(L))));
            end
            L = ~inUnit & ~onB;
            if any(L(:))
                logXhat(L) = conj(skp.logXhatOutCont(1./conj(z(L))));
            end
        end
        
        logX = logXhat;
        if skp.parameter.state ~= paramState.atInf
            logX = logX + 2*log(z - skp.parameter);
        end
        if isempty(skp.parameter.ison)
            logX = logX - skp.normFactor;
        end
    end
    
    function skp = solveBVPs(skp)
        switch skp.parameter.state
            case {paramState.isZero, paramState.atInf, ...
                    paramState.innerFD, paramState.outerFD}
                skp = bvpInDomain(skp);
                
            case paramState.isUnit
                skp = bvpIsUnit(skp);
                
            case paramState.onInnerBdry
                skp = bvpOnInner(skp);
                
            case paramState.onOuterBdry
                skp = bvpOnOuter(skp);
        end
    end
end

methods % Setting and getting
    function nf = get.normFactor(skp)
        nf = skp.refdata.normFactor;
    end
    
    function skp = set.normFactor(skp, nf)
        skp.refdata.normFactor = nf;
    end
    
    function pc = get.primeCorrect(skp)
        pc = skp.refdata.primeCorrect;
    end
    
    function skp = set.primeCorrect(skp, pc)
        skp.refdata.primeCorrect = pc;
    end
    
    function lxh = get.logXhatOutBdry(skp)
        lxh = skp.refdata.logXhatOutBdry;
    end
    
    function skp = set.logXhatOutBdry(skp, lxh)
        skp.refdata.logXhatOutBdry = lxh;
    end
    
    function lxh = get.logXhatOutCont(skp)
        lxh = skp.refdata.logXhatOutCont;
    end
    
    function skp = set.logXhatOutCont(skp, lxh)
        skp.refdata.logXhatOutCont = lxh;
    end
end

end % skprime
