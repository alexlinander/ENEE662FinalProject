% Define the parameters
p=0.5;



for n = [5 7]
    fprintf('For a %d-cycle:\n', n);
    WY_X = zeros(n, n);
    for x = 1:n
        WY_X(x, x) = 1-p;
        y_next = mod(x, n) + 1; 
        WY_X(x, y_next) = p;
    end
    singular_value = calculation(WY_X);
    for m = 1:5
        min_singular_value = singular_value^(m);
        lambda_result = 1/(1-min_singular_value);
        fprintf('  Letter size m = %d: smallest singular value of Q_n ≈ %.7f, corresponding lambda ≈ %.7f\n', m, min_singular_value, lambda_result);
    end
end



function min_singular_value = calculation(WY_X_base)
    n = size(WY_X_base, 1);
    N = n^2;
    WY_Xm = WY_X_base;
    for i = 2:2
        WY_Xm = kron(WY_Xm, WY_X_base);
    end
    pXm = (1/N) * ones(N, 1);
    pXYm = diag(pXm) * WY_Xm;
    pYm = sum(pXYm, 1)';
    min_singular_value=calculate_singular_value_of_Q(pXYm, pXm, pYm);
end

function min_sigma_Qm = calculate_singular_value_of_Q(pXYm, pXm, pYm)
    N = size(pXYm, 1);
    Denom = sqrt(pXm * pYm');
    Qm = zeros(N, N);
    non_zero_denom_indices = abs(Denom) > eps; 
    Qm(non_zero_denom_indices) = pXYm(non_zero_denom_indices) ./ Denom(non_zero_denom_indices);
    singular_values = svd(Qm);
    min_sigma_Qm = singular_values(end);
    
end
