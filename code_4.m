% Define the parameters
s = 0.0000001;
p_set = [
    1/2; 
    1/3; 
    1/4; 
    1/5; 
    1/10; 
    1/100; 
    1/1000;
    1/1000000;
];
n=7;




for m=3:4
    fprintf('For a 7-cycle with letter size m = %d:\n', m);
    for p = p_set'
        WY_X = zeros(n, n);
        for x = 1:n
            WY_X(x, x) = 1-p;
            y_next = mod(x, n) + 1; 
            WY_X(x, y_next) = p;
        end
        [lambda_result, upper_bound] = calculation(WY_X, m, s);
        normalized_upper_bound = upper_bound^(1/m);
        fprintf('  Probability q = %.7f: Largest lambda ≈ %.7f, Upper bound ≈ %.7f, Normalized upper bound ≈ %.7f\n', p, lambda_result, upper_bound, normalized_upper_bound);
    end
end



function [lambda, upper_bound] = calculation(WY_X_base, m, s)
    min_singular_value=calculate_singular_value_of_Q(WY_X_base);
    min_singular_value_power = min_singular_value^(m);
    lambda = 1/(1-min_singular_value_power);
    
    lambda = lambda - s;
        n = size(WY_X_base, 1);
    N = n^(m);
    WY_Xm = WY_X_base;
    for i = 2:m
        WY_Xm = kron(WY_Xm, WY_X_base);
    end
    pXm = (1/N) * ones(N, 1);
    pXYm = diag(pXm) * WY_Xm;
    pYm = sum(pXYm, 1)'; 
    pX_Ym = pXYm ./ pYm;
    H_Xm = -sum(pXm .* log2(pXm), "all", "omitnan");
    H_Xm_Ym = -sum(pXYm .* log2(pX_Ym), "all", "omitnan");
    upper_bound = 2^((H_Xm-lambda*H_Xm_Ym));
end

function min_sigma_Qm = calculate_singular_value_of_Q(WY_X_base)
    n = size(WY_X_base, 1);
    m=2;
    N = n^(m);
    WY_Xm = WY_X_base;
    for i = 2:m
        WY_Xm = kron(WY_Xm, WY_X_base);
    end
    pXm = (1/N) * ones(N, 1);
    pXYm = diag(pXm) * WY_Xm;
    pYm = sum(pXYm, 1)'; 
    N = size(pXYm, 1);
    Denom = sqrt(pXm * pYm');
    Qm = zeros(N, N);
    non_zero_denom_indices = abs(Denom) > eps; 
    Qm(non_zero_denom_indices) = pXYm(non_zero_denom_indices) ./ Denom(non_zero_denom_indices);
    singular_values = svd(Qm);
    min_sigma_Qm = singular_values(end);
    
end
