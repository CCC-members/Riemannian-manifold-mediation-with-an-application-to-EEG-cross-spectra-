function [regularizedB, DOF] = RegMatrixCVestimationDOFcovOpt(totalX, totalY, p, q, s, n, lambda)

    % Preallocate variables
    B = cell(1, 1);
    B{1} = zeros(p, s * q);
    B0 = zeros(p, s * q);
    nrm = (norm(totalX, 2))^2;
    delta = n / nrm;
    alpha = cell(1, 1);
    alpha{1} = 1;
    solution = cell(1, 1);
    temphvalue = cell(1, 1);

    t = 1;

    while t == 1 || t == 2 || abs(temphvalue{t - 1} - temphvalue{t - 2}) / abs(temphvalue{t - 2}) > 1e-4
        if t == 1
            solution{t} = B{t} + (0 - 1) / alpha{t} * (B{t} - B0);
        else
            solution{t} = B{t} + (alpha{t - 1} - 1) / alpha{t} * (B{t} - B{t - 1});
        end

        % Btemp = zeros(p, s * q);
        Atemp = solution{t} - delta * dlossmatrixdividenOpt(totalY, totalX, solution{t}, s, n);
        subAtemp = reshape(Atemp, p, q, s);
        subBtemp = zeros(p, q, s);

        for j = 1:s
            [U, S, V] = svd(subAtemp(:, :, j), 'econ');
            avector = diag(S);
            bvector = max(avector - lambda * delta, 0);
            subBtemp(:, :, j) = U * diag(bvector) * V';
        end

        Btemp = reshape(subBtemp, p, s * q);
        temphvalue{t} = hmatrixfunctiondividen(totalY, totalX, B{t}, lambda, s, n);
        B{t + 1} = Btemp;
        alpha{t + 1} = (1 + sqrt(1 + 4 * alpha{t}^2)) / 2;
        t = t + 1;
    end

    % DOF calculation
    Aspectrum = svd(subAtemp(:, :, 1), 'econ');
    if p ~= q
        Aspectrum(max(p, q)) = 0;
    end
    DOF = 0;
    for i = 1:nnz(bvector)
        DOF = DOF + 1 ...
            + sum(Aspectrum(i) * (Aspectrum(i) - delta * lambda) ...
            ./ (Aspectrum(i)^2 - [Aspectrum(1:i-1); Aspectrum(i+1:p)].^2)) ...
            + sum(Aspectrum(i) * (Aspectrum(i) - delta * lambda) ...
            ./ (Aspectrum(i)^2 - [Aspectrum(1:i-1); Aspectrum(i+1:q)].^2));
    end

    regularizedB = B{t};
end
