function frobNorm = calcFrobNorm(zMatrix)
    zMatrix = zMatrix(:);
    sum = 0;
    for k = 1:numel(zMatrix)
        sum = sum + abs(zMatrix(k))^2;
    end
    frobNorm = sqrt(sum);
end