function result = skptests(subset, selector)
%SKPTESTS runs unit tests for SKPrime.
%
% Unit tests (used as verification tests) are defined as subclasses of
% matlab.unittest.TestCase in the skpUnitTest package.
%
% skptests prime
%   Runs all the available tests for the prime function in the test suite.
%
% skptests greens
%   Runs all the available tests for the Green's functions in the test
%   suite.
%
% skptests
% skptests all
%   Runs all available tests in the SKPrime test suite. This is the
%   default.
%
% skptests <string>
%   Runs a test in the test suite specified by <string>. For example,
%   using string 'skpUnitTest.prime*' is the same as calling
%   `skptests prime`. See the TestSuite.fromPackage documentation
%   for more information on string format.
%
% skptests <string> <string>
% The first string may be appended by another which indicates the set of
% tests to run with `parameterAt` set explicitly by the second string. For
% example
%
%    skptests prime nearCirc1
%
% runs all of the prime tests with `parameterAt=nearCirc1`.
%
% See also: matlab.unittest, matlab.unittest.TestSuite.fromPackage

% Copyright Everett Kropf, 2015
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

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.selectors.HasParameter

runner = TestRunner.withTextOutput('Verbosity', 2);

suiteArgs = {};
if nargin > 0
    switch subset
        case 'prime'
            suiteArgs = {'Name', 'skpUnitTest.prime*'};
        case {'greens'}
            suiteArgs = {'Name', 'skpUnitTest.greens*'};
        case 'all'
            % This is the default above: suiteArgs = {};
        otherwise
            suiteArgs = {'Name', ['skpUnitTest.' subset]};
    end
end

if nargin > 1
    s = HasParameter('Property', 'parameterAt', 'Name', selector);
    suiteArgs = [{s}, suiteArgs];
end

tests = TestSuite.fromPackage('skpUnitTest', suiteArgs{:});

rng('shuffle')
result = run(runner, tests);

fprintf('\nTest run summary:\n\n')
disp(table(result))

if ~nargout
    clear result
end
