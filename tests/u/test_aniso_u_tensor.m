% test_aniso_u_tensor.m
clc;
clear;
close all;
%% main:
% Example usage
U = [0.01, 0.002, 0.001;
     0.002, 0.015, 0.003;
     0.001, 0.003, 0.02]; % Example ADP tensor (symmetric, positive-definite)

random_displacement = generate_random_displacement(U);
disp('Random displacement:');
disp(random_displacement);

%% function:
function displacement = generate_random_displacement(U)
    % Ensure the ADP tensor is symmetric
    if ~isequal(U, U')
        error('The ADP tensor must be symmetric.');
    end
    
    % Diagonalize the ADP tensor
    [V, D] = eig(U);
    
    % Ensure eigenvalues are positive (numerical stability)
    eigenvalues = diag(D);
    if any(eigenvalues <= 0)
        error('All eigenvalues should be positive for a valid ADP tensor.');
    end
    
    % Generate random displacements in the eigenvector basis
    random_displacements = sqrt(eigenvalues) .* randn(3, 1);
    
    % Transform back to the original basis
    displacement = V * random_displacements;
end
