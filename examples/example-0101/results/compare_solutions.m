clc
clear all
close all

%% FILE DEFINITIONS
% 
CASES_LOADING = {'UNIAX'}; % SHEAR; UNIAX
CASES_REFINEMENT = {'COARSE', 'MEDIUM', 'FINE'};
CASES_DIMENSION = {'2D', '3D'};
CASES_INTERPOLATION = [1, 2];
CASES_CONTROL = {'DISPLACEMENT'}; % DISPLACEMENT, FORCE, GRAVITY
% should also test solver types: 0, 1

tight_tol = 1.0e-10;
num_tests = 0;
num_failed_tests = 0;

addpath('../../scripts/compare_abaqus_opencmiss/')

fid_failed  = fopen('failed.tests', 'w');
fid_summary = fopen('results.summary', 'w');


for loadingIdx = 1:size(CASES_LOADING, 2)
    for refinementIdx = 1:size(CASES_REFINEMENT, 2)
        for dimensionIdx = 1:size(CASES_DIMENSION, 2)
            for interpolationIdx = 1:size(CASES_INTERPOLATION, 2)
                for controlIdx = 1:size(CASES_CONTROL, 2)
                    % set current test parameters
                    LOADING = CASES_LOADING{loadingIdx};
                    REFINEMENT = CASES_REFINEMENT{refinementIdx};
                    DIMENSION = CASES_DIMENSION{dimensionIdx};
                    INTERPOLATION = CASES_INTERPOLATION(interpolationIdx);
                    CONTROL = CASES_CONTROL{controlIdx};
                    % run comparison
                    compare_solutions_main
                    % evaluate results
                    num_tests = num_tests + 1;
                    if(tight_tol<l2_norm)
                        num_failed_tests = num_failed_tests + 1;
                        if(num_failed_tests==1)
                            fprintf(fid_failed, ['Failed tests:']);
                        end
                        fprintf(fid_failed, [' ' int2str(num_failed_tests) ' =====']);
                        fprintf(fid_failed, ['    LOADING       = ' LOADING '\n']);
                        fprintf(fid_failed, ['    REFINEMENT    = ' REFINEMENT '\n']);
                        fprintf(fid_failed, ['    DIMENSION     = ' DIMENSION '\n']);
                        fprintf(fid_failed, ['    INTERPOLATION = ' int2str(INTERPOLATION) '\n']);
                        fprintf(fid_failed, ['    CONTROL       = ' CONTROL '\n']);
                    end
                end
            end
        end
    end
end


if(num_failed_tests==0)
    fprintf(fid_failed, 'No failed tests.');
end
fprintf(fid_summary, sprintf('Passed tests: %i / %i', num_tests-num_failed_tests, num_tests));

fclose(fid_failed);
fclose(fid_summary);

%quit
