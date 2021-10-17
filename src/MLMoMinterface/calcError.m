function [meanSquaredError, relNormPercentError] = calcError(ref, pred)
    refSquaredSum = sum(abs(ref(:)).^2);
    squaredErrorSum = sum(abs(ref(:) - pred(:)).^2);
    relNormPercentError = 100* sqrt(squaredErrorSum/refSquaredSum);
    meanSquaredError = squaredErrorSum / numel(ref);
end