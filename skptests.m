function result = skptests(select)
%SKPTESTS runs unit tests for SKPrime.
%
% skptests
%   Runs all tests in the SKPrime test suite.
%
% skptests prime
%   Runs all the tests for the prime function in the test suite.
%
% skptests g0
%   Runs all the tests for the G0 function in the test suite.
%
% skptests <string>
%   Runs a test in the test suite specified by <string>.

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

runner = TestRunner.withTextOutput('Verbosity', 2);

suiteArgs = {};
if nargin
    switch select
        case 'prime'
            suiteArgs = {'Name', 'skpUnitTest.prime*'};
        case {'g0', 'G0'}
            suiteArgs = {'Name', 'skpUnitTest.G0*'};
        otherwise
            suiteArgs = {'Name', ['skpUnitTest.' select]};
    end
end
tests = TestSuite.fromPackage('skpUnitTest', suiteArgs{:});

rng('shuffle')
result = run(runner, tests);

fprintf('\nTest run summary:\n\n')
disp(table(result))

if ~nargout
    clear result
end
