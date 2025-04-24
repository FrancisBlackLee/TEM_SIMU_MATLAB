function displacement = GenerateRandomDisplacement(U, n, tol)
    if nargin == 2
        tol = 1.0e-10;
    end
    % Ensure the ADP tensor is symmetric
%     if ~isequal(U, U')
%         error('The ADP tensor must be symmetric.');
%     end
    if(any(abs(U - U') > tol))
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
