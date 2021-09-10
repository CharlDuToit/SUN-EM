function [weightsModel] = calcWeightsModel(refSelfZmn, refTriZmn, refNonSingZmn,...
    selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp,clusterInd,singDataThresh, addBias, useReal)
    

    % Variable convention:
    % dataset + singularism + ?
    %
    % dataset : ref, pred, unity
    % singularism : nonSing, sing ( 2 or 3 triangles), self, tri
    % ? : Context should be clear, first 2 parts could be missing as well
    
    weightsModel = [];
    emptyWeightsForSingularData = false;
    
 % ------------ SELF TERMS MLR  ------------
    
    [numSelf , ~] = size(selfZmnTerms);
    if (useReal)
        selfZmnTerms = real(selfZmnTerms);
        refSelfZmn = real(refSelfZmn);
    else
        selfZmnTerms = imag(selfZmnTerms);
        refSelfZmn = imag(refSelfZmn);
    end

    unitySelfZmn = sum(selfZmnTerms, 2);
    [selfWeights, predSelfZmn] =MLR(selfZmnTerms, refSelfZmn, singDataThresh, addBias, emptyWeightsForSingularData);
    
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
    [triWeights, predTriZmn] =MLR(triZmnTerms, refTriZmn, singDataThresh, addBias, emptyWeightsForSingularData);

    % ------------ NON SINGULAR TERMS CLUSTERING MLR  ------------
    
    [numNonSing, ~] = size(nonSingZmnProp);
    if (useReal)
        nonSingZmnTerms = real(nonSingZmnTerms);
        refNonSingZmn = real(refNonSingZmn);
    else
        nonSingZmnTerms = imag(nonSingZmnTerms);
        refNonSingZmn = imag(refNonSingZmn);
    end

    unityNonSingZmn  = sum(nonSingZmnTerms, 2); % Same estimation as internal solver   
    [nonSingWeights, predNonSingZmn] = calcClusterWeights(nonSingZmnTerms, refNonSingZmn, clusterInd,singDataThresh, addBias);    
    
    % ------------ ERROR ------------
    
    % self terms
   
    predSelfDiff = predSelfZmn - refSelfZmn;
    %predSelfNormSquareDiff = predSelfDiff.^2 ./ refSelfZmn.^2;
    %predSelfNormMSE = sum(predSelfNormSquareDiff) /numSelf;
    predSelfRelNormPercentError = 100 * sqrt(sum(predSelfDiff.^2) / sum(refSelfZmn.^2));
    predSelfMSE = predSelfDiff'*predSelfDiff /numSelf;
    
    unitySelfDiff =  unitySelfZmn - refSelfZmn;
    %unitySelfNormSquareDiff = unitySelfDiff.^2 ./ refSelfZmn.^2;
    %unitySelfNormMSE = sum(unitySelfNormSquareDiff)./numSelf;
    unitySelfRelNormPercentError = 100 * sqrt(sum(unitySelfDiff.^2) / sum(refSelfZmn.^2));
    unitySelfMSE = unitySelfDiff'*unitySelfDiff ./numSelf;
    
    % 3 unique triangles
   
    predTriDiff = predTriZmn - refTriZmn;
    %predTriNormSquareDiff = predTriDiff.^2 ./ refTriZmn.^2;
    %predTriNormMSE = sum(predTriNormSquareDiff) /numTri;
    predTriRelNormPercentError = 100 * sqrt(sum(predTriDiff.^2) / sum(refTriZmn.^2));
    predTriMSE = predTriDiff'*predTriDiff /numTri;
    
    unityTriDiff =  unityTriZmn - refTriZmn;
    %unityTriNormSquareDiff = unityTriDiff.^2 ./ refTriZmn.^2;
    %unityTriNormMSE = sum(unityTriNormSquareDiff)./numTri;
    unityTriRelNormPercentError = 100 * sqrt(sum(unityTriDiff.^2) / sum(refTriZmn.^2));
    unityTriMSE = unityTriDiff'*unityTriDiff ./numTri;
    
    % Non singular
    
    predNonSingDiff = predNonSingZmn -refNonSingZmn;
    %predNonSingNormSquareDiff = predNonSingDiff.^2 ./ refNonSingZmn.^2;
    %predNonSingNormMSE = sum(predNonSingNormSquareDiff) /numNonSing;
    predNonSingRelNormPercentError = 100 * sqrt(sum(predNonSingDiff.^2) / sum(refNonSingZmn.^2));
    predNonSingMSE = predNonSingDiff'*predNonSingDiff /numNonSing;
    
    unityNonSingDiff =  unityNonSingZmn - refNonSingZmn;
    %unityNonSingNormSquareDiff = unityNonSingDiff.^2 ./ refNonSingZmn.^2;
    %unityNonSingNormMSE = sum(unityNonSingNormSquareDiff)./numNonSing;
    unityNonSingRelNormPercentError = 100 * sqrt(sum(unityNonSingDiff.^2) / sum(refNonSingZmn.^2));
    unityNonSingMSE = unityNonSingDiff'*unityNonSingDiff ./numNonSing;
    
    
    % ------------ PLOT DATA ------------
    % for instant debugging
    
    gridSize = 500;
    
    %---SURFACE PLOTS
    
    if (~useReal)
    
        
        %[sortedRefNonSingZmn, ind] = sort(refNonSingZmn);
        %[sortedRefNonSingZmn, ind] = sort(log(abs(refNonSingZmn)));
        %[sortedRefNonSingZmn, ind] = sort(refNonSingZmn.^2);
        
        %------2D PLOTS----
        
        % NON SINGULAR PLOTS 2D PLOTS
        
        %plot(sortedRefNonSingZmn,'g','LineWidth', 5 )
        hold on
        
        %plot(predNonSingZmn(ind),'b--', 'LineWidth', 1)
        %plot(log(abs(predNonSingZmn(ind))),'b', 'LineWidth', 1)
        %plot(predNonSingZmn(ind).^2,'b', 'LineWidth', 1)
        
        %plot(unityNonSingZmn(ind), 'LineWidth', 1)
        %plot(log(abs(unityNonSingZmn(ind))), 'LineWidth', 1)
        %plot(unityNonSingZmn(ind).^2, 'LineWidth', 1)
        
        hold off;
        hold on;
%        [sortedPredNonSingNormSquareDiff, ind1] = sort(predNonSingNormSquareDiff);
        %[sortedUnityNonSingNormSquareDiff, ind2] = sort(unityNonSingNormSquareDiff);
        %[sortedPredNonSingSquareDiff, ind3] = sort(predNonSingDiff.^2);
        %[sortedUnityNonSingSquareDiff, ind4] = sort(unityNonSingDiff.^2);
        
        %plot(log(sortedPredNonSingNormSquareDiff), 'r.','MarkerSize', 1)
        %plot(log(sortedUnityNonSingNormSquareDiff), 'b.','MarkerSize', 1)
        
        %plot( predNonSingDiff(ind2).^2 , 'r.','MarkerSize', 1)
        %plot(log(sortedUnityNonSingNormSquareDiff), 'b.','MarkerSize', 1)
        
        %plot(log(sortedPredNonSingNormSquareDiff), 'r.','MarkerSize', 1)
        %plot(sortedPredNonSingNormSquareDiff, 'b','LineWidth', 1)
        
        %plot(unityNonSingNormSquareDiff(ind),'r','LineWidth', 1)
        %plot(log(unityNonSingNormSquareDiff(ind)),'b.','MarkerSize', 1)
        %plot(sort(log(unityNonSingNormSquareDiff)),'b.','MarkerSize', 1)
        %ylim([0,10]);
        hold off;
        
        hold on;
        %[sortedPredNonSingSquareDiff, ind] = sort(predNonSingDiff.^2);
        %plot(log(sortedPredNonSingSquareDiff), 'r.','MarkerSize', 1)
        
        %plot(log(unityNonSingDiff(ind).^2), 'b.','MarkerSize', 1)
        %plot(sort(log(unityNonSingDiff(ind).^2)), 'b.','MarkerSize', 1)
        hold off;
        hold on;
              
        %log norm error
        % we want ponits below 0 to show improvement
        %plot(-log(sortedUnityNonSingNormSquareDiff) + log(predNonSingNormSquareDiff(ind2)) , 'r.','MarkerSize', 1)
        
        
        %plot( log(predNonSingNormSquareDiff(ind2)) , 'r.','MarkerSize', 1)
        %plot(log(sortedUnityNonSingNormSquareDiff), 'b.','MarkerSize', 1)
        %numImprov = find(-log(sortedUnityNonSingNormSquareDiff) + log(predNonSingNormSquareDiff(ind2)) < 0);
        %36k
        
        % We want points above 0 to show improvemeent
        %plot(-log(sortedPredNonSingNormSquareDiff) + log(unityNonSingNormSquareDiff(ind1)) , 'r.','MarkerSize', 1)
       % plot(log(sortedPredNonSingNormSquareDiff), 'b.','MarkerSize', 1)
       % mlmom.numImprovedPoints = find(log(sortedPredNonSingNormSquareDiff) - log(unityNonSingNormSquareDiff(ind1)) > 0);
        
        %plot(log(predNonSingNormSquareDiff(ind)), 'r.','MarkerSize', 1)
        
        %plot(sort(log(predNonSingNormSquareDiff)), 'r.','MarkerSize', 1)
        
        hold off;
        hold on;
        %plot(sortedUnityNonSingSquareDiff, 'r.','MarkerSize', 1)
        %a = predNonSingDiff.^2;
        %plot( a(ind4), 'b.','MarkerSize', 1)
        
        hold off;
        
         % SELF TERM PLOTS 2D PLOTS
         
        % [sortedPredSelfNormSquareDiff, ind5] = sort(predSelfNormSquareDiff);
        % [sortedUnitySelfNormSquareDiff, ind6] = sort(unitySelfNormSquareDiff);
         %[sortedPredSelfSquareDiff, ind7] = sort(predSelfDiff.^2);
         %[sortedUnitySelfSquareDiff, ind8] = sort(unitySelfDiff.^2);
         hold on;
        % plot(sortedUnitySelfSquareDiff, 'r.','MarkerSize', 1)
         %a = predSelfDiff.^2;
         %plot(-sortedUnitySelfSquareDiff + a(ind8), 'b.','MarkerSize', 1)
         hold off
         
         % TRI TERM PLOTS 2D PLOTS
         
         %[sortedPredTriNormSquareDiff, ind9] = sort(predTriNormSquareDiff);
        % [sortedUnityTriNormSquareDiff, ind10] = sort(unityTriNormSquareDiff);
         %[sortedPredTriSquareDiff, ind11] = sort(predTriDiff.^2);
         %[sortedUnityTriSquareDiff, ind12] = sort(unityTriDiff.^2);
         hold on;
         %plot(sortedUnityTriSquareDiff, 'r.','MarkerSize', 1)
         %a = predTriDiff.^2;
        % plot(-sortedUnityTriSquareDiff + a(ind12), 'b.','MarkerSize', 1)
         hold off
    end

    % ------------ UPDATE MLMOM ------------
    % ------------
    % Matrices
    %weightsStruct.clusterAvg = clusterAvg;
    %weightsStruct.clusterInd = clusterInd; 
    %weightsStruct.dirDotDirInterval = dirDotDirInterval;
    %weightsStruct.dirDotDispInterval =dirDotDispInterval;
    %weightsStruct.distInterval = distInterval;
    %weightsStruct.nonSingZmnProp = nonSingZmnProp;
    
    % Non singular
    weightsModel.nonSingZmnTerms = nonSingZmnTerms;  
    weightsModel.refNonSingZmn = refNonSingZmn;
    weightsModel.unityNonSingZmn = unityNonSingZmn;
    weightsModel.predNonSingZmn = predNonSingZmn; 
    weightsModel.nonSingWeights = nonSingWeights;
    
    % Self
    weightsModel.selfZmnTerms = selfZmnTerms;  
    weightsModel.refSelfZmn = refSelfZmn;
    weightsModel.unitySelfZmn = unitySelfZmn;
    weightsModel.predSelfZmn = predSelfZmn; 
    weightsModel.selfWeights = selfWeights;
    
    % 3 unique triangles
    weightsModel.triZmnTerms = triZmnTerms;  
    weightsModel.refTriZmn = refTriZmn;
    weightsModel.unityTriZmn = unityTriZmn;
    weightsModel.predTriZmn = predTriZmn;
    weightsModel.triWeights = triWeights;
    
    % ------------
    % Scalars
    
    %mlmom.totsetupTime = 0.0;
    %mlmom.totsolTime = 0.0;
    weightsModel.numNonSing = numNonSing;
    
    % Non singular
    %weightsModel.predNonSingNormMSE = predNonSingNormMSE;
    %weightsModel.unityNonSingNormMSE = unityNonSingNormMSE; 
    weightsModel.predNonSingRelNormPercentError = predNonSingRelNormPercentError;
    weightsModel.unityNonSingRelNormPercentError = unityNonSingRelNormPercentError;
    weightsModel.predNonSingMSE = predNonSingMSE;  
    weightsModel.unityNonSingMSE = unityNonSingMSE;

    % Self
    %weightsModel.predSelfNormMSE = predSelfNormMSE;
    %weightsModel.unitySelfNormMSE = unitySelfNormMSE;
    weightsModel.predSelfRelNormPercentError = predSelfRelNormPercentError;
    weightsModel.unitySelfRelNormPercentError = unitySelfRelNormPercentError;
    weightsModel.predSelfMSE = predSelfMSE;    
    weightsModel.unitySelfMSE = unitySelfMSE;
    
    % 3 unique triangles
    %weightsModel.predTriNormMSE = predTriNormMSE;
    %weightsModel.unityTriNormMSE = unityTriNormMSE;
    weightsModel.predTriRelNormPercentError = predTriRelNormPercentError;
    weightsModel.unityTriRelNormPercentError = unityTriRelNormPercentError;
    weightsModel.predTriMSE= predTriMSE;     
    weightsModel.unityTriMSE = unityTriMSE;
    
end
