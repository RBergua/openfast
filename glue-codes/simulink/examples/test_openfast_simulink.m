%% Test OpenFAST Simulink Interface
classdef test_openfast_simulink < matlab.unittest.TestCase
    
    %% Test Method Block
    methods (Test)
        
        function testOpenLoopRuns(testCase)

            workspace_root = getenv("GITHUB_WORKSPACE");

            this_file_path = fileparts(which(mfilename()));

            cd(this_file_path);
            
            % these variables are defined in the OpenLoop model's FAST_SFunc block:
            FAST_InputFileName = fullfile(workspace_root, 'reg_tests', 'r-test', 'glue-codes', 'openfast', 'AOC_WSt', 'AOC_WSt.fst');
            TMax               = 5; % seconds
            
            sim('OpenLoop.mdl',[0,TMax]);
          
        end
    end
end
