function [lambdaMax, nullMSE] = computeLambdaMax(X, Y, weights, alpha, standardize)
%
% lambdaMax is the penalty term (lambda) beyond which coefficients
% are guaranteed to be all zero.
%
% nullMse is the mse of the fit using just a constant term.
% It is provided in this function as a convenience, because it needs 
% to be calculated in the same context as lambdaMax whenever
% lambdaMax is calculated.

if ~isempty(weights)
    observationWeights = true;
    weights = weights(:)';        
    % Normalized weights are used for standardization and calculating lambdaMax.
    normalizedweights = weights / sum(weights);
else
    observationWeights = false;
end

[N,~] = size(X);

% If we were asked to standardize the predictors, do so here because
% the calculation of lambdaMax needs the predictors as we will use
% them to perform fits.

if standardize
    % If X has any constant columns, we want to protect against
    % divide-by-zero in normalizing variances.
    constantPredictors = (range(X)==0);

    if ~observationWeights
        % Center and scale
        [X0,~,~] = zscore(X,1);
    else
        % Weighted center and scale
        muX = normalizedweights * X;
        X0 = bsxfun(@minus,X,muX);
        sigmaX = sqrt( normalizedweights * (X0.^2) );
        % Avoid divide by zero with constant predictors
        sigmaX(constantPredictors) = 1;
        X0 = bsxfun(@rdivide, X0, sigmaX);
    end
else
    if ~observationWeights
        % Center
        muX = mean(X,1);
        X0 = bsxfun(@minus,X,muX);
    else
        % Weighted center
        muX = normalizedweights(:)' * X;
        X0 = bsxfun(@minus,X,muX);
    end
end

% If using observation weights, make a weighted copy of the 
% predictor matrix, for use in weighted dot products.

if observationWeights
    wX0 = bsxfun(@times, X0, weights');
end

if ~observationWeights
    muY = mean(Y);
else
    muY = weights*Y;
end
% Y0 = bsxfun(@minus,Y,muY);
Y0 = Y - muY;

% Calculate max lambda that permits non-zero coefficients
%
if ~observationWeights
    dotp = abs(X0' * Y0);
    lambdaMax = max(dotp) / (N*alpha);
else
    dotp = abs(sum(bsxfun(@times, wX0, Y0)));
    lambdaMax = max(dotp) / alpha;
end

if ~observationWeights
    nullMSE = mean(Y0.^2);
else
    % This works because weights are normalized and Y0 is already
    % weight-centered.
    nullMSE = weights * (Y0.^2);
end
end