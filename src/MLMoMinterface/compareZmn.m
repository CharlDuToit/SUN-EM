%function [comp] = compareZmn(refZmn, predZmn, unityZmn, singInd, freqIndex, useReal)
function [comp] = compareZmn(refZmn, predZmn, unityZmn, singInd)
   %Matrices are all [num_edges, num_edges, 1] ie only 1 frequency
   comp = [];
   [numEdges, ~]  = size(refZmn);
   nonSingCount = 0;
   triCount = 0;
   predTwoUniqueSquaredDiffSum = 0;
   unityTwoUniqueSquaredDiffSum = 0;
   refTwoUniqueSquaredSum = 0;
   predThreeUniqueSquaredDiffSum = 0;
   unityThreeUniqueSquaredDiffSum = 0;
   refThreeUniqueSquaredSum = 0;
   predFourUniqueSquaredDiffSum = 0;
   unityFourUniqueSquaredDiffSum = 0;
   refFourUniqueSquaredSum = 0;
   
   predSquaredSum = 0;
   unitySquaredSum = 0;

   for mm = 1:numEdges
       for nn = 1:numEdges
           predSquaredSum = predSquaredSum + abs(predZmn(mm,nn)).^2;
           unitySquaredSum = unitySquaredSum + abs(unityZmn(mm,nn)).^2;
           if (mm == nn)
               
               predSelfDiff = abs(predZmn(mm,nn) - refZmn(mm,nn));
               predTwoUniqueSquaredDiffSum = predTwoUniqueSquaredDiffSum + predSelfDiff.^2 ;
   
               %predSelfRelNormPercentError = predSelfDiff.^2 / refSelfZmn.^2;
               %predSelfMSE = predSelfDiff'*predSelfDiff /numSelf;
               
               unitySelfDiff =  abs(unityZmn(mm,nn)  - refZmn(mm,nn));
               unityTwoUniqueSquaredDiffSum = unityTwoUniqueSquaredDiffSum + unitySelfDiff.^2;
               
               refTwoUniqueSquaredSum = refTwoUniqueSquaredSum + abs(refZmn(mm,nn)).^2;
               %unitySelfMSE = unitySelfDiff'*unitySelfDiff ./numSelf;
           elseif (singInd(mm,nn))
               triCount = triCount + 1;
               predTriDiff = abs(predZmn(mm,nn) - refZmn(mm,nn));
               predThreeUniqueSquaredDiffSum = predThreeUniqueSquaredDiffSum + predTriDiff.^2 ;
               
               unityTriDiff =  abs(unityZmn(mm,nn)  - refZmn(mm,nn));
               unityThreeUniqueSquaredDiffSum = unityThreeUniqueSquaredDiffSum + unityTriDiff.^2;
               
               refThreeUniqueSquaredSum = refThreeUniqueSquaredSum + abs(refZmn(mm,nn)).^2;
               
           else
               nonSingCount = nonSingCount + 1;
               predNonSingDiff = abs(predZmn(mm,nn) - refZmn(mm,nn));
               predFourUniqueSquaredDiffSum = predFourUniqueSquaredDiffSum + predNonSingDiff.^2 ;
               
               unityNonSingDiff =  abs(unityZmn(mm,nn)  - refZmn(mm,nn));
               unityFourUniqueSquaredDiffSum = unityFourUniqueSquaredDiffSum + unityNonSingDiff.^2;
               
                refFourUniqueSquaredSum = refFourUniqueSquaredSum + abs(refZmn(mm,nn)).^2;
           end
       end % for nn
   end % for mm
   
   predRowRelNormPercentError = zeros(numEdges,1);
   for mm = 1:numEdges
       [~, predRowRelNormPercentError(mm)] = calcError(refZmn(mm,:), predZmn(mm,:));
   end
         
   comp.predFourUniqueRelNormPercentError = 100 * sqrt(predFourUniqueSquaredDiffSum / refFourUniqueSquaredSum);
   comp.unityFourUniqueRelNormPercentError = 100 * sqrt(unityFourUniqueSquaredDiffSum / refFourUniqueSquaredSum);
   comp.predFourUniqueMSE = predFourUniqueSquaredDiffSum/nonSingCount;
   comp.unityFourUniqueMSE = unityFourUniqueSquaredDiffSum/nonSingCount;
   
   comp.predTwoUniqueRelNormPercentError = 100 * sqrt(predTwoUniqueSquaredDiffSum / refTwoUniqueSquaredSum);
   comp.unityTwoUniqueRelNormPercentError = 100 * sqrt(unityTwoUniqueSquaredDiffSum / refTwoUniqueSquaredSum);
   comp.predTwoUniqueMSE = predTwoUniqueSquaredDiffSum/numEdges;
   comp.unityTwoUniqueMSE = unityTwoUniqueSquaredDiffSum/numEdges;
   
   comp.predThreeUniqueRelNormPercentError = 100 * sqrt(predThreeUniqueSquaredDiffSum / refThreeUniqueSquaredSum);
   comp.unityThreeUniqueRelNormPercentError = 100 * sqrt(unityThreeUniqueSquaredDiffSum / refThreeUniqueSquaredSum);
   comp.predThreeUniqueMSE = predThreeUniqueSquaredDiffSum/triCount;
   comp.unityThreeUniqueMSE = unityThreeUniqueSquaredDiffSum/triCount;

   comp.refFrobNorm = sqrt(refFourUniqueSquaredSum + refThreeUniqueSquaredSum + refTwoUniqueSquaredSum);
   comp.predFrobNorm = sqrt(predSquaredSum);
   comp.unityFrobNorm = sqrt(unitySquaredSum);
   
   [~, comp.predRelNormPercentError] = calcError(refZmn, predZmn);
    comp.predRowRelNormPercentError = predRowRelNormPercentError;
    
   [~, comp.unityRelNormPercentError] = calcError(refZmn, unityZmn);
   
end
