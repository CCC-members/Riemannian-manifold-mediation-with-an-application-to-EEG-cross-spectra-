function [dlossvalue] = dlossmatrixdividenOpt(Y, X, B, s, n)
    [row, col] = size(B);
    subcol = col / s;

    % Preallocate arrays for efficiency
    subB = cell(1, s);
    subX = cell(1, s);
    tempsum = zeros(row * size(X, 1), subcol, s);

    % Use a single loop and preallocate memory
    for j = 1:s
        subB{j} = B(:, ((j - 1) * subcol + 1):(j * subcol));
        subX{j} = X(:, j);
        tempsum(:, :, j) = kron(subX{j}, subB{j});
    end

    % Sum along the third dimension
    fit = sum(tempsum, 3);
    residual = Y - fit;

    % Preallocate tempdloss1 and use vectorized operations
    tempdloss1 = zeros(row, subcol, s, n);

    for j = 1:s
        for i = 1:n
            tempdloss1(:, :, j, i) = -2 * X(i, j) * residual((i - 1) * row + 1:i * row, :);
        end
    end

    % Sum along the fourth dimension
    tempdloss = sum(tempdloss1, 4);
    dlossvalue = reshape(tempdloss, row, col) / n / 2;
end
