function [X, M1, y, p, q] = DataGenerator_complexCovariatesOptimized(n, w)
tic
    % error correlation between voxels
    arnumber = 0.5;

    % error variance
    errorvariance = 1;

    load('matrixB3.mat', 'matrixB3');
    load('matrixB4.mat', 'matrixB4');

    [p, q] = size(matrixB3);
    
    % true coefficient matrices
    matrixA1 = zeros(p, q);
    matrixA1(31:34, 17:30) = 1;
    matrixA1(17:48, 31:34) = 1;
    matrixA1(31:34, 35:48) = 1;

    matrixA2 = zeros(p, q);
    matrixA2(17:48, 17:48) = 1;

    shape1 = imread('cross.gif');
    shape1 = imresize(shape1, [32, 32]); % Use imresize instead of array_resize
    matrixB5 = zeros(2 * size(shape1));
    matrixB5((size(matrixB5, 1) / 4):(size(matrixB5, 1) / 4) + size(shape1, 1) - 1, ...
        (size(matrixB5, 2) / 4):(size(matrixB5, 2) / 4) + size(shape1, 2) - 1) = shape1;
    matrixB5(matrixB5 ~= 0) = 1;

    shape1 = imread('Asterisk.gif');
    shape1 = array_resize(shape1,[32, 32]); % 32-by-32
    matrixB6 = zeros(2*size(shape1));
    matrixB6((size(matrixB6,1)/4):(size(matrixB6,1)/4)+size(shape1,1)-1, ...
        (size(matrixB6,2)/4):(size(matrixB6,2)/4)+size(shape1,2)-1) = shape1;
    matrixB6(matrixB6~=0) = 1.4142+1.4142i;

    % Complex Matrix Figures for simulations
    trianguleFigureCmplx    = matrixB3;
    butterflyFigureCmplx    = matrixB4;
    myCrossFigureCmplx      = matrixB5;
    myAsteriskCmplx         = matrixB6;

    trianguleFigureCmplx(trianguleFigureCmplx ~= 0) = 1.4142 + 1.4142i;
    butterflyFigureCmplx(butterflyFigureCmplx ~= 0) = 1.4142 + 1.4142i;
    myCrossFigureCmplx(myCrossFigureCmplx ~= 0)     = 1.4142 + 1.4142i;
    myAsteriskCmplx(myAsteriskCmplx ~= 0)           = 1.4142 + 1.4142i; 

    truecoefficientA = zeros(p, q, w);
    truecoefficientA(:, :, 1) = myAsteriskCmplx;
    truecoefficientA(:, :, 2) = trianguleFigureCmplx;
    truecoefficientA(:, :, 3) = butterflyFigureCmplx;
    truecoefficientA(:, :, 4) = myCrossFigureCmplx;

    truecoefficientB = myAsteriskCmplx;

    % Generate correlation matrix among X
    sigma = 0.5 .^ abs((1:w)' - (1:w));
    RR = chol(sigma);

    % Generate correlation matrix between pixels of the random error
    pq = p * q;
    [I, J] = meshgrid(1:pq, 1:pq);
    errorsigma = arnumber .^ (abs(floor((I - 1) / p) - floor((J - 1) / p)) + abs(mod(I - 1, p) - mod(J - 1, p)));
    sigmaRR = chol(errorsigma);

    % M = Ax + D1*w
    RandStream.setGlobalStream(RandStream('mt19937ar', 'seed', 2018));

    M = zeros(n * p, q);
    x = zeros(n, 1);
    x(randperm(n, ceil(n / 2))) = 1;

    X = zeros(n, w);
    tempcomponent = zeros(p, q, w);
    % arrayM = zeros(n, p, q);
    % sumtempM = zeros(p, q, n);

    for i = 1:n
        XX = randn(1, w);
        Xvector = XX * RR;
        Xvector(1) = x(i); % Give values of zeros and ones to x
        X(i, :) = Xvector;

        for j = 1:w
            tempcomponent(:, :, j) = Xvector(j) * truecoefficientA(:, :, j);
        end
        fitsum = sum(tempcomponent, 3);
        error = reshape(reshape(errorvariance * randn(p, q), 1, pq) * sigmaRR, p, q);
        M((i - 1) * p + 1:i * p, :) = fitsum + error;
        % arrayM(i, :, :) = M((i - 1) * p + 1:i * p, :);
        % sumtempM(:, :, i) = M((i - 1) * p + 1:i * p, :);
    end

    % y = BM + d1*w
    M1 = permute(reshape(M.', q, p, []), [2, 1, 3]); % Convert from 2D to 3D

    % True coefficients for regular (non-array) covariates
    b0 = ones(w, 1);
    M2 = tensor(randn(p, q, n)); % p1-by-p2-by-n matrix variates

    if isreal(M1)
        M_real_im = M1;
    else
        M_real_im = [real(M1); imag(M1)];
        truecoefficientB = [real(truecoefficientB); imag(truecoefficientB)];
    end

    Mt = tensor(M_real_im);
    mu = X * b0 + double(ttt(tensor(truecoefficientB), Mt, 1:2));

    sigma = 1; % noise level
    error = sigma * randn(n, 1);
    y = mu + error;
toc
end
