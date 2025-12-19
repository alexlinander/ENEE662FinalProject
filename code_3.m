p=0.5;
lambda_matrix = [
    1.1055728, 1.0092025, 1.0008715, 1.0000831, 1.0000079;
    1.0520950, 1.0024578, 1.0001214, 1.0000060, 1.0000002
];

for n = [5 7]
    fprintf('For a %d-cycle:\n', n);
    WY_X = zeros(n, n);
    if n==5
        i=1;
    else
        i=2;
    end
    for x = 1:n
        WY_X(x, x) = 1-p;
        y_next = mod(x, n) + 1; 
        WY_X(x, y_next) = p;
    end
    for m = 1:5
        upper_bound = calculation(WY_X, m, lambda_matrix(i, m));
        normalized_upper_bound = upper_bound^(1/m);
        fprintf('  Letter size m = %d: upper bound ≈ %.7f, normalized upper bound ≈ %.7f\n', m, upper_bound, normalized_upper_bound);
    end
end

function upper_bound= calculation(WY_X_base, m, lambda_star)
    n = size(WY_X_base, 1);
    N = n^m;
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
    upper_bound = 2^((H_Xm-lambda_star*H_Xm_Ym));
end
