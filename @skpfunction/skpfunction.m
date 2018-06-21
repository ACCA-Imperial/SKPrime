classdef skpfunction < bvpFun
%SKPFUNCTION is the Shottky-Klein prime function.
%
% skpfunction(alpha, dv, qv)
% skpfunction(alpha, D)
%   Solve the boundary value problem for the SK-prime function given the
%   parameter alpha and the vectors dv and qv, respectively the centers and
%   radii of the circles bounded by the unit circle. The domain may also be
%   given as a skpDomain object D.
%
% w = skpfunction(...)
% w2 = skpfunction(param, w)
%   Use a previously computed prime function to accelerate the computation
%   of a new prime function on the same domain with a different parameter.
%
% skpfunction(..., N)
%   Specify the truncation level of the Fourier series on each boundary
%   circle (same for all circles). See bvpFun.truncation for default.
%
% Garbage parameters are ignored, i.e., if the constructor has enough
% information to build a valid object, there will be no warning given about
% unused extra parameters.
%
% See also: skpDomain, bvpFun

% Everett Kropf, 2015, 2016
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
    parameter               % Prime function "parameter"
    vjFuns = {}             % Cell array of first kind integral objects
end

properties(Access=protected)
    refdata                 % Shared memory state storage
    logXhatBdry             % Boundary function handle
    logXhatCont             % Continuation function handle
end

properties(Dependent, Access=protected)
    primeCorrect            % Root branch correction flag
    normFactor              % Normalisation factor
    logXhatOutBdry          % Outer boundary function handle
    logXhatOutCont          % Outer domain function handle
end

methods
    function skp = skpfunction(alpha, varargin) %(alpha, dv, qv, N)
        if nargin
            [D, N, phi, vjfuns] = skpfunction.parseArguments(varargin{:});
            sargs = {D, N, phi};
        else
            sargs = {};
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
        % We're done if the domain is simply connected.
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
                skp.vjFuns{j} = vjFirstKind(j, skp);
            end
        elseif numel(vjfuns) ~= m
                error('SKPrime:invalidArgument', ...
                    'Incorrect number of v_j functions given.')
        else
            skp.vjFuns = vjfuns;
        end
        
        skp = solveBVPs(skp);
    end % ctor
    
    function dw = diff(skp, n)
        %derivative of prime function with respect to variable z.
        %
        % dw = diff(skp)
        % dw = diff(skp, n)
        
        if skp.domain.m == 0
            dw = @(z) complex(ones(size(z)));
            return
        end
        
        % FIXME: validate 1 <= n <= 3.
        if nargin < 2
            n = 1;
        end
        if n > 3
            error('SKPrime:invalidArgument', ...
                'Derivatives only supported to 3rd order.')
        end
        
        din = dftDerivative(skp, @skp.feval, n);
        
        wi = invParam(skp);        
        alpha = skp.parameter;
        switch n
            case 1
                dwi = dftDerivative(skp, @wi.feval);
                if alpha.state == paramState.isZero ...
                        || alpha.state == paramState.atInf
                else
                    dout = @(z) -alpha*conj(wi.feval(1./conj(z))) ...
                        + alpha*conj(dwi(1./conj(z)))./z;
                end
                
            case 2
                d2wi = dftDerivative(skp, @wi.feval, 2);
                if alpha.state == paramState.isZero ...
                        || alpha.state == paramState.atInf
                else
                    dout = @(z) -alpha*conj(d2wi(1./conj(z)))./z.^3;
                end
                
            case 3
                d2wi = dftDerivative(skp, @wi.feval, 2);
                d3wi = dftDerivative(skp, d2wi);
                if alpha.state == paramState.isZero ...
                        || alpha.state == paramState.atInf
                else
                    dout = @(z) 3*alpha*conj(d2wi(1./conj(z)))./z.^4 ...
                        + alpha*conj(d3wi(1./conj(z)))./z.^5;
                end
        end
        
        function v = dwEval(z)
            v = complex(nan(size(z)));
            
            mask = abs(z) <= 1 + eps;
            if any(mask(:))
                v(mask) = din(z(mask));
            end
            if any(~mask(:))
                v(~mask) = dout(z(~mask));
            end
        end
        
        dw = @dwEval;
    end
    
    function dwh = diffh(skp, n)
        %derivative of prime "hat" function wrt z.
        %
        % dwh = diffh(skp)
        % dwh = diffh(skp, n)
        
        % FIXME: Validate 1 <= n <= 3.
        if nargin < 2
            n = 1;
        end
        if n > 3
            error('SKPrime:invalidArgument', ...
                'Derivatives only supported to 3rd order.')
        end
        
        if skp.domain.m == 0
            dwh = @(z) complex(zeros(size(z)));
            return
        end
        
        din = dftDerivative(skp, @skp.hat, n);
        
        wi = invParam(skp);
        dwi = dftDerivative(skp, @wi.hat);
        switch n
            case 1
                dout = @(z) -conj(dwi(1./conj(z)))./z.^2;
                
            case 2
                d2wi = dftDerivative(skp, dwi);
                dout = @(z) 2*conj(dwi(1./conj(z)))./z.^3 ...
                            + conj(d2wi(1./conj(z)))./z.^4;
                
            case 3
                d2wi = dftDerivative(skp, dwi);
                d3wi = dftDerivative(skp, d2wi);
                dout = @(z) -6*conj(dwi(1./conj(z)))./z.^4 ...
                            - 6*conj(d2wi(1./conj(z)))./z.^5 ...
                            - conj(d3wi(1./conj(z)))./z.^6;
        end
        
        function v = dwhEval(z)
            v = complex(nan(size(z)));
            
            inUnit = abs(z) <= 1 + eps(2);
            if any(inUnit(:))
                v(inUnit) = din(z(inUnit));
            end
            if any(~inUnit(:))
                v(~inUnit) = dout(z(~inUnit));
            end
        end
        
        dwh = @dwhEval;
    end
    
    function dwp = diffp(skp)
        %derivative of the prime function wrt the parameter.
        %Note the returned function is a function of the variable
        %with the parameter *fixed*.
        %
        %See also skpfunction/diffXp.
        
%         warning('SKPrime:underDevelopment', ...
%             'This method is under development. Do not trust its output.')
        
        dXp = diffXp(skp);
        dwp = @(z) dXp(z)./feval(skp,z)/2;
    end
    
    function dX = diffX(skp, n)
        %Variable derivative of X to nth order (up to 3).
        %
        % dX = diffX(w)
        % dX = diffX(w, n)
        
        if nargin < 2
            n = 1;
        end
        
        dw = diff(skp);
        if n == 1
            dX = @(z) 2*skp.feval(z).*dw(z);
        else
            d2w = diff(skp, 2);
            if n < 3
                dX = @(z) 2*(dw(z).^2 + skp.feval(z).*d2w(z));
            else
                d3w = diff(skp, 3);
                dX = @(z) 2*(3*dw(z).*d2w(z) + skp.feval(z).*d3w(z));
            end
        end
    end
    
    function dXh = diffXh(skp, n)
        %Variable derivative of X.hat to nth order (up to 3).
        %
        % dXh = diffXh(w)
        % dXh = diffXh(w, n)
        
        if nargin < 2
            n = 1;
        end
        
        dwh = diffh(skp);
        if n == 1
            dXh = @(z) 2*skp.hat(z).*dwh(z);
        else
            d2wh = diffh(skp, 2);
            if n < 3
                dXh = @(z) 2*(dwh(z).^2 + skp.hat(z).*d2wh(z));
            else
                d3wh = diffh(skp, 3);
                dXh = @(z) 2*(3*dwh(z).*d2wh(z) + skp.hat(z).*d3wh(z));
            end
        end
    end
    
    function dXp = diffXp(skp)
        %derivative of the square of the prime function wrt the parameter.
        %Note the returned function is a function of the variable
        %with the parameter *fixed*.
        %
        %See also skpfunction/diffp.
        
        alpha = skp.parameter;
        if alpha == 0
            error('SKPrime:undefinedState', ...
                ['The case when alpha is at the origin is undefined ' ...
                'as currently implemented.'])
        end
        if ~isempty(alpha.ison)
            error('SKPrime:undefinedState', ...
                ['The case when alpha is on a boundary is undefined ' ...
                'as currently implemented.'])
        end
        
        dag0 = greensC0Dp(alpha, skp);
        din = @(z) (4i*pi*dag0(z) + 1/alpha).*skp.X(z);
        
        dig0 = greensC0Dp(1/conj(alpha), skp);
        wi = invParam(skp);
        dXip = @(z) (4i*pi*dig0(z) + conj(alpha)).*wi.X(z);
        dout = @(z) z.^2.*(2*alpha*conj(wi.X(1./conj(z))) ...
            - conj(dXip(1./conj(z))));
        
        function v = dXpEval(z)
            v = complex(nan(size(z)));
            
            mask = abs(z) <= 1;
            if any(mask(:))
                v(mask) = din(z(mask));
            end
            if any(~mask(:))
                v(~mask) = dout(z(~mask));
                
                warning('SKPrime:undefinedState', ...
                    ['Evaluation of points outside the unit disk ' ...
                    'is currently under development. Results of such ' ...
                    'evaluation is not to be trusted.'])
            end
        end
        
        dXp = @dXpEval;
    end
    
    function skp = invParam(skp)
        %gives the prime function with inverted parameter
        %
        %  w = skpfunction(...);
        %  wi = invParam(w);
        %  Then isa(wi, 'skpfunction') is a true statememt.
        
        skp = skpinvparam(skp);
    end
    
    function w = feval(skp, z)
        %provides prime function evaluation
        %
        % w = feval(skp, z)

        if skp.parameter.state ~= paramState.atInf
            w = z - skp.parameter;
        else
            w = 1;
        end
        if skp.domain.m == 0
            return
        end
        
        w = w.*primeHat(skp, z);
    end
    
    function wh = hat(skp, z)
        %evaluates the prime function with the zero factored out.
        
        wh = primeHat(skp, z);
    end
    
    function X = X(skp, z)
        %returns the square of the prime function
        %
        % w = skpfunction(...);
        % Xval = X(w, z);
        
        if skp.parameter.state ~= paramState.atInf
            X = (z - skp.parameter).^2;
        else
            X = 1;
        end
        if skp.domain.m == 0
            return
        end
        
        X = X.*primeXhat(skp, z);
    end
    
    function Xh = Xhat(skp, z)
        %evaluate Xhat at z.
        
        Xh = primeXhat(skp, z);
    end
end

methods(Access=protected, Sealed)
    function wh = primeHat(skp, z)
        %prime function with zero factored out.
        %
        %Provided here so this functionality can't be overriden.
        %
        %See skpinvparam.
        
        if skp.domain.m == 0
            wh = complex(zeros(size(z)));
            return
        end
        
        [wh, L] = evalLogXhat(skp, z);
        wh = exp(wh/2);
        
        if ~isempty(skp.primeCorrect) && skp.primeCorrect
            if abs(skp.parameter) < 1 + eps(2)
                L = ~L;
            end
            wh(L) = -wh(L);
        end
    end
    
    function Xh = primeXhat(skp, z)
        %prime function square with zeros factored out.
        %
        %Provided here so this functionality can't be overriden.
        %
        %See skpinvparam.
        
        if skp.domain.m == 0
            Xh = complex(zeros(size(z)));
            return
        end
        
        Xh = exp(evalLogXhat(skp, z));
    end
end

methods(Access=protected)
    function skps = copyProperties(skp)
        %copies all object properties to structure
        %
        % Gives external access to a copy of object internals.
        
        mco = ?skpfunction;
        for k = 1:numel(mco.PropertyList)
            pname = mco.PropertyList(k).Name;
            skps.(pname) = skp.(pname);
        end
    end
    
    function [logXhat, inUnit] = evalLogXhat(skp, z)
        %evaluate the log of Xhat.
        
        logXhat = complex(nan(size(z)));
        
        inUnit = abs(z) < 1 + eps(2);
        if any(inUnit(:))
            if ~isempty(skp.logXhatOutCont) && isempty(skp.primeCorrect)
                setCorrection(skp);
            end
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
            if isempty(skp.logXhatOutCont)
                solveBVPouter(skp);
            end
            if isempty(skp.primeCorrect)
                setCorrection(skp);
            end
            
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
            
        if isempty(skp.parameter.ison)
            logXhat = logXhat - skp.normFactor;
        end
    end
    
    function skp = setCorrection(skp)
        %does square root branch correction
        %
        % Does the prime function need a correction to stay on the proper branch?
        % FIXME: This should check the test point is not in a circle! This should
        % have a maximum retry count on test point selection!
        
        testpt = 1 - 1e-6;
        while abs(testpt - skp.parameter) < 1e-2
            testpt = (1 - 1e-6)*exp(2i*pi*rand(1));
        end
        skp.primeCorrect = abs(imag(skp.logXhatCont(testpt) ...
            - conj(skp.logXhatOutCont(1/conj(testpt))))) > pi/4;
    end
end

methods(Access=protected) % BVP stuff
    function skp = bvpInDomain(skp)
        %solves BVP for parameter in unit disk
        
        alpha = skp.parameter;
        psi = calcPsi(skp, alpha);

        imLogXhat = primeRHS(alpha, skp.vjFuns);
        phi = solve(skp.phiFun, imLogXhat);
        skp.logXhatBdry = @(z) phi(z) + 1i*imLogXhat(z) + 2*log(psi(z));
        skp.logXhatCont = ...
            SKP.bmcCauchy(skp.logXhatBdry, skp.domain, skp.truncation);
        if abs(alpha) < 1 + eps(2)
            skp.normFactor = skp.logXhatCont(alpha);
        else
            skp = bvpInDomainOuter(skp);
        end
    end
    
    function skp = bvpInDomainOuter(skp)
        %solves BVP for parameter in unit disk for outer evaluation
        
        alpha = skp.parameter;
        psi = calcPsi(skp, inv(alpha));

        imLogXhatInvp = primeRHS(inv(alpha), skp.vjFuns);
        phiInvp = solve(skp.phiFun, imLogXhatInvp);
        skp.logXhatOutBdry = ...
            @(z) phiInvp(z) + 1i*imLogXhatInvp(z) + 2*log(psi(z));
        skp.logXhatOutCont = ...
            SKP.bmcCauchy(skp.logXhatOutBdry, skp.domain, skp.truncation);
        if abs(alpha) >= 1 + eps(2)
            skp.normFactor = conj(skp.logXhatOutCont(1/conj(alpha)));
        end
    end
    
    function skp = bvpIsUnit(skp)
        %solves BVP for parameter on unit circle
        
        alpha = skp.parameter;

        imLogXhat = primeRHS(alpha, skp.vjFuns);
        phi = solve(skp.phiFun, imLogXhat);
        skp.normFactor = phi(alpha) + 1i*imLogXhat(alpha);
        skp.logXhatBdry = @(z) phi(z) + 1i*imLogXhat(z) - skp.normFactor;
        skp.logXhatCont = ...
            SKP.bmcCauchy(skp.logXhatBdry, skp.domain, skp.truncation);
    end
    
    function skp = bvpIsUnitOuter(skp)
        %sovles BVP for unit parameter and outer evaluation
        
        skp.logXhatOutBdry = skp.logXhatBdry;
        skp.logXhatOutCont = SKP.bmcCauchy(skp.logXhatOutBdry, ...
            skp.domain, skp.truncation);
    end
    
    function skp = bvpOnInner(skp)
        %solves BVP for parameter on inner boundary
        
        alpha = skp.parameter;
        D = skp.domain;
        N = skp.truncation;
        [d, q] = domainDataB(D);
        j = alpha.ison;
        vj = skp.vjFuns{j};
        
        imLogXhat = primeRHS(alpha, skp.vjFuns);        
        phi = solve(skp.phiFun, imLogXhat);
        skp.normFactor = phi(alpha) + 1i*imLogXhat(alpha);
        skp.logXhatBdry = @(z) phi(z) + 1i*imLogXhat(z) - skp.normFactor;
        skp.logXhatCont = SKP.bmcCauchy(skp.logXhatBdry, D, N);
                
        if alpha ~= 0
            logInvpPart = @(z) 2*log((z - alpha)./(z - 1/conj(alpha))) ...
                + 4i*pi*(real(vj(alpha)) - vj(z)) ...
                - 2*log(q(j+1)/(1 - conj(d(j+1)/alpha)));
            skp.logXhatOutBdry = @(z) logInvpPart(z) + skp.logXhatBdry(z);
            skp.logXhatOutCont = @(z) logInvpPart(z) + skp.logXhatCont(z);
        else
            imLogXhatInvp = primeRHS(inv(alpha), skp.vjFuns);
            phiInvp = solve(skp.phiFun, imLogXhatInvp);
            skp.logXhatOutBdry = @(z) phiInvp(z) + 1i*imLogXhatInvp(z);
            skp.logXhatOutCont = SKP.bmcCauchy(skp.logXhatOutBdry, D, N);
        end        
    end
    
    function skp = bvpOnOuter(skp)
        %solves BVP for parameter on outer boundary
        
        alpha = skp.parameter;
        [d, q] = domainDataB(skp.domain);
        j = alpha.ison;
        vj = skp.vjFuns{j};
        
        imLogXhatInvp = primeRHS(inv(alpha), skp.vjFuns);        
        phi = solve(skp.phiFun, imLogXhatInvp);
        skp.normFactor = phi(1/conj(alpha)) + 1i*imLogXhatInvp(1/conj(alpha));
        skp.logXhatOutBdry = @(z) phi(z) + 1i*imLogXhatInvp(z) ...
            - skp.normFactor;
        skp.logXhatOutCont = ...
            SKP.bmcCauchy(skp.logXhatOutBdry, skp.domain, skp.truncation);
        
        logInvpPart = @(z) 2*log((z - 1/conj(alpha))./(z - alpha)) ...
            + 4i*pi*(real(vj(alpha)) - vj(z)) ...
            - 2*log(q(j+1)/(1 - conj(d(j+1))*alpha));
        skp.logXhatBdry = @(z) logInvpPart(z) + skp.logXhatOutBdry(z);
        skp.logXhatCont = @(z) logInvpPart(z) + skp.logXhatOutCont(z);
        
        skp = setCorrection(skp);
    end
    
    function psi = calcPsi(skp, param)
        %computes parameter near boundary correction
        
        [d, q] = domainData(skp.domain);
        
        psi = @(z) 1;
        for j = isclose(skp.domain, param);
            thj = @(z) d(j) + q(j)^2*z./(1 - conj(d(j))*z);
            if abs(param) == 0
                thjp = d(j);
            elseif isinf(param)
                thjp = d(j) - q(j)^2/conj(d(j));
            else
                thjp = thj(param);
            end
            psi = @(z) psi(z).*(z - thjp)./(z - thj(z));
        end
    end
    
    function skp = solveBVPs(skp)
        %calls proper BVP solver method
        
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
    
    function skp = solveBVPouter(skp)
        %calls proper BVP solver method for outer evaluation
        
        switch skp.parameter.state
            case {paramState.isZero, paramState.atInf, ...
                    paramState.innerFD, paramState.outerFD}
                bvpInDomainOuter(skp);
                
            case paramState.isUnit
                skp = bvpIsUnitOuter(skp);
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

methods(Static,Hidden)
    function [D, N, phi, vjfuns] = parseArguments(varargin)
        %Assume a minimum of one argument given.
        %
        %Cases:
        %  1) dv, qv, [N]
        %  2) skpfunction, [N]
        %  3) vjFunctions, [N]
        %  4) skpDomain, [N]
        %  5) bvp object, [N]
        
        phi = [];
        vjfuns = {};
        Ndx = 2;
        
        D = varargin{1};
        if isa(D, 'skpfunction')
            vjfuns = D.vjFuns;
            phi = D.phiFun;
            D = D.domain;
        elseif iscell(D) ...
                && all(cellfun(@(c) isa(c, 'vjFirstKind'), D(:)))
            vjfuns = D;
            phi = D.phiFun;
            D = D.domain;
        elseif isa(D, 'bvpFun')
            phi = D.phiFun;
            D = D.domain;
        elseif isa(D, 'skpDomain') || isa(D, 'circleRegion')
            D = skpDomain(D);
        else
            qv = varargin{2};
            D = skpDomain(D, qv);
            Ndx = 3;
        end
        
        if nargin == Ndx
            N = varargin{Ndx};
        else
            N = [];
        end
    end
end

end % skpfunction
