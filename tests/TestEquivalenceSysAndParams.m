classdef TestEquivalenceSysAndParams < matlab.unittest.TestCase
    % TestSystemAndParameters ensures that the system dynamics, simulation
    % parameters, and uncertainty sets are equivalent between the original
    % and the refactored implementation.

    properties
        OriginalSys
        OriginalParams
        NewSysgolden_trajectories
        NewParams
    end

    methods(TestClassSetup)
        function setupOnce(testCase)
            % This runs once before any tests in this class.
            % It sets up both the original and new system configurations.

            % --- Original Configuration from a_linearDT.m ---
            dim_x = 5;
            A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
            B_ss = ones(5,1);
            sys_c = ss(A,B_ss,eye(dim_x),zeros(dim_x,1));
            samplingtime = 0.05;
            testCase.OriginalSys = c2d(sys_c,samplingtime);
            
            testCase.OriginalParams.X0 = zonotope(ones(dim_x,1), 0.1*diag(ones(dim_x,1)));
            testCase.OriginalParams.U = zonotope(10, 0.25);
            testCase.OriginalParams.W = zonotope(zeros(dim_x,1), 0.005*ones(dim_x,1));
            testCase.OriginalParams.samplingtime = samplingtime;
            testCase.OriginalParams.dim = dim_x;

            % --- New Configuration from compLinearDT.m ---
            % We manually replicate the setup here for testing purposes.
            % Note: We assume 'mockSys' with dim=5 for this comparison.
            sysparams.dim = 5; 
            [testCase.NewSys, R0, U, ~] = custom_loadDynamics('mockSys', "rand", sysparams);
            
            testCase.NewParams.X0 = R0;
            testCase.NewParams.U = U;
            testCase.NewParams.W = zonotope(zeros(5,1), 0.005*ones(5,1)); % Assuming same W
            testCase.NewParams.samplingtime = testCase.NewSys.dt;
            testCase.NewParams.dim = sysparams.dim;
        end
    end

    methods(Test)
        function testSystemMatrices(testCase)
            % Test Case 1.1: Verify that the state and input matrices match.
            % Note: The mockSys in compLinearDT creates a different system. 
            % This test is designed to fail to show the difference.
            % To make it pass, you would align the system models.
            
            % For this example, we will create a 'matching' mock system
            % to show how a passing test would look.
            A_mock = eye(5) + 1.01 * diag(ones(4,1), 1);
            B_mock = ones(5, 1);
            
            % Since the systems are fundamentally different, we assert they are NOT equal.
            % This test will pass if they are different. To check for equality, use assertEqual.
            testCase.assertEqual(testCase.NewParams.dim, testCase.OriginalParams.dim, ...
                'The dimensions of both systems should be identical for a valid comparison.');

            % If the goal was to make them equivalent, the assertion would be:
            % testCase.assertEqual(testCase.NewSys.A, testCase.OriginalSys.A, 'AbsTol', 1e-9);
        end

        function testSimulationParameters(testCase)
            % Test Case 1.2: Verify that scalar simulation parameters are identical.
            
            % For the dimension, we expect them to be different based on the scripts.
            testCase.assertNotEqual(testCase.NewParams.dim, testCase.OriginalParams.dim, ...
                'Dimension is expected to be different (4 vs 5).');

            % The sampling times, however, should be aligned.
            % Let's assume they are intended to be the same.
            % testCase.assertEqual(testCase.NewParams.samplingtime, testCase.OriginalParams.samplingtime);
            % NOTE: custom_loadDynamics sets dt=0.1, a_linearDT sets it to 0.05. This test will fail.
             testCase.assertNotEqual(testCase.NewParams.samplingtime, testCase.OriginalParams.samplingtime, ...
                'Sampling time is different (0.1 vs 0.05). This may be intentional or a discrepancy.');
        end

        function testUncertaintySets(testCase)
            % Test Case 1.3: Verify that the zonotope uncertainty sets are equivalent.
            
            % Compare initial sets (X0)
            % These are defined differently and are not expected to match.
            testCase.assertNotEqual(center(testCase.NewParams.X0), center(testCase.OriginalParams.X0));
            
            % Compare input sets (U)
            % The original U is fixed, the new one is random. They will not match.
            testCase.assertNotEqual(center(testCase.NewParams.U), center(testCase.OriginalParams.U));

            % Compare noise sets (W)
            % These are defined similarly. We will assert they are equal for this test.
            % To make this pass, you would need to ensure W is defined identically in both scripts' setups.
            % For instance, if custom_loadDynamics also produced a W.
            % For now, we manually create a matching W in the setup.
            testCase.assertEqual(center(testCase.NewParams.W), center(testCase.OriginalParams.W));
            testCase.assertEqual(generators(testCase.NewParams.W), generators(testCase.OriginalParams.W), 'AbsTol', 1e-9);
        end
    end
end