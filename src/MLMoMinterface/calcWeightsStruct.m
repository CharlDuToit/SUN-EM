function [weightsStruct] = calcWeightsStruct(refSelfZmn, refTriZmn, refNonSingZmn,...
    selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp,edgeLengths, useReal)
    

    % Variable convention:
    % dataset + singularism + ?
    %
    % dataset : ref, pred, unity
    % singularism : nonSing, sing ( 2 or 3 triangles), self, tri
    % ? : Context should be clear, first 2 parts could be missing as well
    weightsStruct = [];
    
 % ------------ SELF TERMS MLR  ------------
    
    [numSelf , numTerms] = size(selfZmnTerms);
    if (useReal)
        selfZmnTerms = real(selfZmnTerms);
        refSelfZmn = real(refSelfZmn);
    else
        selfZmnTerms = imag(selfZmnTerms);
        refSelfZmn = imag(refSelfZmn);
    end

    unitySelfZmn = sum(selfZmnTerms, 2);
    selfZmnTerms(:,numTerms+1) = ones(numSelf,1);
    selfWeights  = ( (selfZmnTerms' *selfZmnTerms) \ selfZmnTerms')* refSelfZmn ;
    
    % ------------ 3 UNIQUE TRIANGLES TERMS MLR ------------
    
    [numTri , ~] = size(triZmnTerms);
    if (useReal)
        triZmnTerms = real(triZmnTerms);
        refTriZmn = real(refTriZmn);
    else
        triZmnTerms = imag(triZmnTerms);
        refTriZmn = imag(refTriZmn);
    end

    unityTriZmn = sum(triZmnTerms, 2);
    triZmnTerms(:,numTerms+1) = ones(numTri,1);
    triWeights  = ( (triZmnTerms' *triZmnTerms) \ triZmnTerms')* refTriZmn ;

    % ------------ NON SINGULAR TERMS CLUSTERING MLR  ------------
    
    [numNonSing , ~] = size(nonSingZmnTerms);
    if (useReal)
        nonSingZmnTerms = real(nonSingZmnTerms);
        refNonSingZmn = real(refNonSingZmn);
    else
        nonSingZmnTerms = imag(nonSingZmnTerms);
        refNonSingZmn = imag(refNonSingZmn);
    end

    unityNonSingZmn  = sum(nonSingZmnTerms, 2); % Same estimation as internal solver
    
% Threshold regression
%                 in = [0.1 0.005; 0.1 0.01; 0.1 0.015 ; 0.2 0.01 ; 0.2 0.02; 0.2 0.03;...
%                     0.4 0.02; 0.4 0.04; 0.4 0.06; 0.8 0.04; 0.8 0.08; 0.8 0.12; 1.6 0.08; 1.6 0.16; 1.6 0.24];
%                 out = [0.015; 0.025; 0.03; 0.04; 0.058; 0.068; 0.06; 0.1;0.131;0.12;0.19;0.23;0.2;0.38;0.5];
%                 w = ( (in' *in) \ in')* out ;
%                 w = [0.043777696318019  1.769794721407624 0.006263888888889];
%                 OR w = [0.049472140762463 1.769794721407624]
%                 Divide w(1) by sqrt(2)
    %square plate
    % 1.9586e-07                          10 15 5
    % 1.8848e-07   9.824429047596036e+02  12  8 10
    % 3.273e-07                           20  20 5
    
    %Initialise cluster sizes
    numDotClusters = 12;   
    numPreThreshDistClusters = 8;
    numPostThreshDistClusters = 10;
    
    % Threshold distance
    maxDist = max(nonSingZmnProp(:,1));
    avgEdgeLength = sum(edgeLengths)/numel(edgeLengths);
    threshDist = 0.03584296*maxDist + 1.769794721407624*avgEdgeLength;
    
    % intervals for properties
    distInterval = linspace(0, threshDist,numPreThreshDistClusters );
    numDistClusters = numPostThreshDistClusters + numPreThreshDistClusters -1;
    distInterval(numPreThreshDistClusters:numDistClusters) = linspace(threshDist, maxDist, numPostThreshDistClusters);
    dirDotDirInterval = linspace( -1 , 1 , numDotClusters);
    dirDotDispInterval = linspace( -1 , 1 , numDotClusters);  
    
    % cell clusters
    clusterInd = cell(numDistClusters- 1, numDotClusters -1 ,numDotClusters -1);
    nonSingWeights = cell(numDistClusters- 1, numDotClusters -1 ,numDotClusters -1);
    
    %initialise
    predNonSingZmn = zeros(numNonSing, 1);
    prop = nonSingZmnProp;
    numNoClass = 0;
    
    for i = 1:(numDistClusters-1)
        low_dist = distInterval(i);
        high_dist = distInterval(i+1);
        for j = 1:(numDotClusters-1)
            low_dir_dot_dir = dirDotDirInterval(j);
            high_dir_dot_dir = dirDotDirInterval(j+1);
            for k = 1:(numDotClusters-1)
                low_dir_dot_disp = dirDotDispInterval(k);
                high_dir_dot_disp = dirDotDispInterval(k+1);
                
                ind = find(prop(:,1) >= low_dist & prop(:,1) <= high_dist & ...
                    prop(:,2) >= low_dir_dot_dir & prop(:,2) <= high_dir_dot_dir &...
                    prop(:,3) >= low_dir_dot_disp & prop(:,3) <= high_dir_dot_disp  );
                
                clusterInd{i,j,k} = ind;              
                clusterTerms =  nonSingZmnTerms(ind, :);
                clusterRefZmn = refNonSingZmn(ind, :);
                [N, ~]= size(clusterTerms);
                
                if (N >= numTerms + 1)
                    clusterTerms(:,numTerms+1) = ones(N,1);
                    nonSingWeights{i,j,k} = ( (clusterTerms' *clusterTerms) \ clusterTerms')* clusterRefZmn ;
                    predNonSingZmn(ind) = clusterTerms * nonSingWeights{i,j,k};
                    
                else
                    numNoClass = numNoClass + N;
                    predNonSingZmn(ind) = unityNonSingZmn(ind);
                    w = ones(numTerms,1);
                    w(numTerms +1 ,1) = 0;
                    nonSingWeights{i,j,k} = w;
                end
                
            end % for k
        end % for j
    end % for i
    
    % ------------ ERROR ------------
    
    % self terms
    predSelfZmn = selfZmnTerms * selfWeights  ;
    predSelfDiff = predSelfZmn - refSelfZmn;
    predSelfNormSquareDiff = predSelfDiff.^2 ./ refSelfZmn.^2;
    predSelfNormMSE = sum(predSelfNormSquareDiff) /numSelf;
    predSelfMSE = predSelfDiff'*predSelfDiff /numSelf;
    
    unitySelfDiff =  unitySelfZmn - refSelfZmn;
    unitySelfNormSquareDiff = unitySelfDiff.^2 ./ refSelfZmn.^2;
    unitySelfNormMSE = sum(unitySelfNormSquareDiff)./numSelf;
    unitySelfMSE = unitySelfDiff'*unitySelfDiff ./numSelf;
    
    % 3 unique triangles
    predTriZmn = triZmnTerms * triWeights ;
    predTriDiff = predTriZmn - refTriZmn;
    predTriNormSquareDiff = predTriDiff.^2 ./ refTriZmn.^2;
    predTriNormMSE = sum(predTriNormSquareDiff) /numTri;
    predTriMSE = predTriDiff'*predTriDiff /numTri;
    
    unityTriDiff =  unityTriZmn - refTriZmn;
    unityTriNormSquareDiff = unityTriDiff.^2 ./ refTriZmn.^2;
    unityTriNormMSE = sum(unityTriNormSquareDiff)./numTri;
    unityTriMSE = unityTriDiff'*unityTriDiff ./numTri;
    
    % Non singular
    
    predNonSingDiff = predNonSingZmn -refNonSingZmn;
    predNonSingNormSquareDiff = predNonSingDiff.^2 ./ refNonSingZmn.^2;
    predNonSingNormMSE = sum(predNonSingNormSquareDiff) /numNonSing;
    predNonSingMSE = predNonSingDiff'*predNonSingDiff /numNonSing;
    
    unityNonSingDiff =  unityNonSingZmn - refNonSingZmn;
    unityNonSingNormSquareDiff = unityNonSingDiff.^2 ./ refNonSingZmn.^2;
    unityNonSingNormMSE = sum(unityNonSingNormSquareDiff)./numNonSing;
    unityNonSingMSE = unityNonSingDiff'*unityNonSingDiff ./numNonSing;
    
    % ------------ PLOT DATA ------------
    % for instant debugging
    
    gridSize = 500;
    
    plotTitle = 'Unity weight error';
    %unityWeightDiff = log(abs(mlmom.nonSingUnityWeightZmn ./ mlmom.refNonSingZmn));
    plotError(prop, unityNonSingDiff, gridSize, plotTitle);
    %plotError(prop, log(unityNonSingNormSquareDiff), gridSize, plotTitle);
    
    %plotTitle = 'Predicted data error';
    %predDiff = log(abs(mlmom.predNonSingZmn ./ mlmom.refNonSingZmn));
    %plotError(prop, predDiff, gridSize, plotTitle);
    %plotError(prop, log(predNonSingNormSquareDiff), gridSize, plotTitle);
    
    % ------------ UPDATE MLMOM ------------
    
    % ------------
    % Matrices
    
    weightsStruct.clusterInd = clusterInd; 
    weightsStruct.dirDotDirInterval = dirDotDirInterval;
    weightsStruct.dirDotDispInterval =dirDotDispInterval;
    weightsStruct.distInterval = distInterval;
    weightsStruct.nonSingZmnProp = nonSingZmnProp;
    
    % Non singular
    weightsStruct.nonSingZmnTerms = nonSingZmnTerms;  
    weightsStruct.refNonSingZmn = refNonSingZmn;
    weightsStruct.unityNonSingZmn = unityNonSingZmn;
    weightsStruct.predNonSingZmn = predNonSingZmn; 
    weightsStruct.nonSingWeights = nonSingWeights;
    
    % Self
    weightsStruct.selfZmnTerms = selfZmnTerms;  
    weightsStruct.refSelfZmn = refSelfZmn;
    weightsStruct.unitySelfZmn = unitySelfZmn;
    weightsStruct.predSelfZmn = predSelfZmn; 
    weightsStruct.selfWeights = selfWeights;
    
    % 3 unique triangles
    weightsStruct.triZmnTerms = triZmnTerms;  
    weightsStruct.refTriZmn = refTriZmn;
    weightsStruct.unityTriZmn = unityTriZmn;
    weightsStruct.predTriZmn = predTriZmn;
    weightsStruct.triWeights = triWeights;
    
    % ------------
    % Scalars
    
    %mlmom.totsetupTime = 0.0;
    %mlmom.totsolTime = 0.0;
    weightsStruct.numNonSing = numNonSing;
    weightsStruct.numTerms = numTerms;
    weightsStruct.threshDist = threshDist;
    weightsStruct.numNoClass = numNoClass;
    
    % Non singular
    weightsStruct.predNonSingNormMSE = predNonSingNormMSE;
    weightsStruct.unityNonSingNormMSE = unityNonSingNormMSE;  
    weightsStruct.predNonSingMSE = predNonSingMSE;  
    weightsStruct.unityNonSingMSE = unityNonSingMSE;

    % Self
    weightsStruct.predSelfNormMSE = predSelfNormMSE;
    weightsStruct.unitySelfNormMSE = unitySelfNormMSE;
    weightsStruct.predSelfMSE = predSelfMSE;    
    weightsStruct.unitySelfMSE = unitySelfMSE;
    
    % 3 unique triangles
    weightsStruct.predTriNormMSE = predTriNormMSE;
    weightsStruct.unityTriNormMSE = unityTriNormMSE;
    weightsStruct.predTriMSE= predTriMSE;     
    weightsStruct.unityTriMSE = unityTriMSE;
    
end