classdef skpParameter < double
%skpParameter is the S-K prime function parameter.
%
% aobj = skpParameter(alpha, D),
%    alpha is a scalar parameter in the fundamental domain,
%    D is a skpDomain object.
%
% Properties:
%    state is the parameter state enumerated in paramState.
%    ison is empty unless the parameter is on a boundary, in which case
%        this is the integer number of the boundary, 0:m.

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

% FIXME: Until subsref automatically reads properties from metadata, it
% must be modified to match property definitions.
properties(SetAccess=protected)
    state               % paramState object
    ison                % boundary number where parameter is located
    indisk              % disk containing parameter
end

properties(Dependent)
    inUnitDomain            % boolean value for abs(alpha) <= 1
end

methods
    function aobj = skpParameter(alpha, D)
        if ~nargin
            alpha = [];
        end
        alpha = skpParameter.normalizeParameter(alpha);
        aobj = aobj@double(alpha);
        if ~nargin
            return
        end
        
        if isa(D, 'skpParameter')
            aobj.state = D.state;
            aobj.ison = D.ison;
            return
        end
        
        if ~isa(D, 'skpDomain')
            error('SKPrime:invalidArgument', 'Expected a "skpDomain" object.')
        end
        
        [d, q] = domainDataB(D);
        
        function ja = whichHoleNoFail(a)
            ja = find(abs(a - d(2:end)) < q(2:end) + eps(2));
            assert(~isempty(ja) && numel(ja) == 1, ...
                'SKPrime:logicError', ...
                'Code logic problem, submit bug report.')
        end
        function bool = isinFD(a)
            bool = isin(D, a) | isin(D, 1/conj(a));
        end
        function errorOutsideDomain()
            error('SKPrime:invalidArgument', ...
                ['The given parameter is outside the ' ...
                'acceptable region. See help.'])
        end
        
        ja = [];
        if abs(1 - abs(alpha)) < eps(2)
            aloc = paramState.isUnit;
            aon = 0;
        elseif abs(alpha) >= 1 + eps(2)
            % Outside unit disk.
            aon = abs(q - abs(1/conj(alpha) - d)) < eps(2);
            if any(aon)
                aon = find(aon) - 1;
                aloc = paramState.onOuterBdry;
            else
                aon = [];

                % Is it in a hole or the domain?
                if ~isin(D, 1/conj(alpha))
                    % In a 1st level reflection?
                    ja = whichHoleNoFail(1/conj(alpha));
                    if isinFD(D.theta(ja, alpha))
                        aloc = paramState.outerDisk;
                    else
                        errorOutsideDomain()
                    end
                else
                    if isinf(alpha)
                        aloc = paramState.atInf;
                    else
                        aloc = paramState.outerFD;
                    end
                end
            end
        else
            % Inside unit disk.
            aon = abs(q - abs(alpha - d)) < eps(2);
            if any(aon)
                aon = find(aon) - 1;
                aloc = paramState.onInnerBdry;
            else
                aon = [];
                if ~isin(D, alpha)
                    % Check if it's in a 1st level reflection.
                    ja = whichHoleNoFail(alpha);
                    if isinFD(D.theta(ja, 1/conj(alpha)))
                        aloc = paramState.innerDisk;
                    else
                        errorOutsideDomain()
                    end
                else
                    if abs(alpha) == 0
                        aloc = paramState.isZero;
                    else
                        aloc = paramState.innerFD;
                    end
                end
            end
        end
        
        aobj.state = aloc;
        aobj.ison = aon;
        aobj.indisk = ja;
    end
    
    function disp(aobj)
        disp(double(aobj))
    end
    
    function aobj = inv(aobj)
        aobj.state = inv(aobj.state);
        aobj = skpParameter(1/conj(aobj), aobj);
    end

    function out = subsref(aobj, S)
        if S(1).type == '.'
            switch S(1).subs
                case 'state'
                    out = aobj.state;
                    return
                    
                case 'ison'
                    out = aobj.ison;
                    return
                    
                case 'indisk'
                    out = aobj.indisk;
                    return
                    
                case 'inUnitDomain'
                    out = aobj.inUnitDomain;
                    return
                    
                otherwise
                    error('SKPrime:invalidReference', ...
                        'Parameter ''%s'' not a member of skpParameter.', ...
                        S(1).subs)
            end
        elseif strcmp(S(1).type, '()') && ~isempty(S(1).subs)
            out = subsref(double(aobj), S);
            return
        end
        error('SKPrime:invalidReference', ...
            'Not a supported subscripted reference.')
    end
    
    function out = horzcat(varargin)
        varargin = cellfun(@double, varargin, 'UniformOutput', false);
        out = horzcat(varargin{:});
    end
    
    function out = vertcat(varargin)
        varargin = cellfun(@double, varargin, 'UniformOutput', false);
        out = vertcat(varargin{:});
    end
    
    function bool = get.inUnitDomain(alpha)
        %true if abs(parameter) <= 1.
        
        bool = alpha.state <= 0;
    end
end

methods(Static)
    function param = normalizeParameter(param)
        if isempty(param)
            return
        end
        
        if ~(isnumeric(param) && numel(param) == 1)
            error('SKPrime:invalidArgument', ...
                'The parameter must be a single, numeric scalar.')
        end

        if 0 < abs(param) && abs(param) < eps
            param = 0;
            warning('SKPrime:normalising', ...
                'Treating "small" alpha as zero.')
        elseif ~isinf(param) && 1/eps < abs(param)
            param = inf;
            warning('SKPrime:normalising', ...
                'Treating "large" alpha as infinity.')
        end
    end
end

end
