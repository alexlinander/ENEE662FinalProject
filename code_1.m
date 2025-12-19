% Define the parameters
s = 0.0000001;
p=0.5;



for n = [5 7]
    fprintf('For a %d-cycle:\n', n);
    WY_X = zeros(n, n);
    for x = 1:n
        WY_X(x, x) = 1-p;
        y_next = mod(x, n) + 1; 
        WY_X(x, y_next) = p;
    end
    for m = 1:5
        lambda_result = calculation(WY_X, m, s);
        fprintf('  Letter size m = %d: largest lambda ≈ %.7f\n', m, lambda_result);
    end
end



function lambda = calculation(WY_X_base, m, s)
    n = size(WY_X_base, 1);
    N = n^m;
    WY_Xm = WY_X_base;
    for i = 2:m
        WY_Xm = kron(WY_Xm, WY_X_base);
    end
    pXm = (1/N) * ones(N, 1);
    pXYm = diag(pXm) * WY_Xm;
    pYm = sum(pXYm, 1)'; 
    pX_Ym = zeros(N, N);
    non_zero_py_indices = pYm > eps; 
    pX_Ym(:, non_zero_py_indices) = pXYm(:, non_zero_py_indices) ./ pYm(non_zero_py_indices)';
    KX = diag(pXm);
    KY = pXYm * pX_Ym';
    
    lambda = 1.0; 
    while true
        H = (lambda - 1) * KX - lambda * KY;
        max_eigenvalue = max(eig(H));
        if max_eigenvalue > 1e-16
            break;
        else 
            lambda = lambda + s;
        end
    end
    lambda = lambda - s;
end
