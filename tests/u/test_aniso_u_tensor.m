% test_aniso_u_tensor.m
clc;
clear;
close all;
%% main:
% Example usage
U = [0.050, 0.025, 0.000;
     0.025, 0.050, 0.000;
     0.000, 0.000, 0.030]; % Example ADP tensor (symmetric, positive-definite)

random_displacement = generate_random_displacement(U, 100000);
figure;
scatter(random_displacement(1, :), random_displacement(2, :), 'filled');
axis equal;
axis tight;
% disp('Random displacement:');
% disp(random_displacement);

%% function:
function displacement = generate_random_displacement(U, n)
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
    random_displacements = sqrt(eigenvalues) .* randn(3, n);
    
    % Transform back to the original basis
    displacement = V * random_displacements;
end
