function testSuite = createTestSuiteFromData(sys, out_ddra)
% createTestSuiteFromData - converts DDRA trajectory data into CORA testSuite

dt = sys.dt;
x_all = out_ddra.X_model{1}.center; % X0 center
n = sys.dims.n;
m = sys.dims.m;

% Generate dummy trajectory data
steps = size(out_ddra.X_model, 1) - 1;
u_nom = repmat(10 * ones(1,m), steps, 1); % matches DDRA input center
y_nom = zeros(steps+1, n); % dummy, not used in conformance reach
y_nom(1,:) = x_all'; % initial state

for i = 2:steps+1
    % Just put placeholder output if not needed
    y_nom(i,:) = y_nom(i-1,:) + 0.1 * randn(1,n);
end

% Create test case
testSuite{1} = testCase(y_nom, u_nom, x_all, dt, 'linearSysDT');
end
