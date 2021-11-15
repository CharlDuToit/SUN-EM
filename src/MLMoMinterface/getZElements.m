function z = getZElements(indices, zMat, clusterInd)
    % ind: num x 2, each row = m,n
    % z = getZElements(mlmom.indicesStruct.fourUniqueIndices,zMatrices.values,13);
    % avg = sum(z)/numel(z);
    %ind = find(imag(z) > 0);
    %z(ind)
    [ind, ~] = find(indices(:,3) == clusterInd);
    indices = indices(ind, 1:2);
    
    [~,~, numCol] = size(zMat); 
    num = numel(indices(:,1));
    z = zeros(num, numCol);
    for k = 1:num
        z(k, :) = zMat(indices(k,1), indices(k,2), :);
    end
end