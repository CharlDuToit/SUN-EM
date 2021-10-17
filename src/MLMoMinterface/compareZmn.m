%function [comp] = compareZmn(refZmn, predZmn, unityZmn, singInd, freqIndex, useReal)
function [comp] = compareZmn(refZmn, predZmn, unityZmn, singInd)
   %Matrices are all [num_edges, num_edges, 1] ie only 1 frequency
   comp = [];
   [numEdges, ~]  = size(refZmn);
   nonSingCount = 0;
   triCount = 0;
   predSelfSquaredDiffSum = 0;
   unitySelfSquaredDiffSum = 0;
   refSelfSquaredSum = 0;
   predTriSquaredDiffSum = 0;
   unityTriSquaredDiffSum = 0;
   refTriSquaredSum = 0;
   predNonSingSquaredDiffSum = 0;
   unityNonSingSquaredDiffSum = 0;
   refNonSingSquaredSum = 0;
   
   predSquaredSum = 0;
   unitySquaredSum = 0;

   for mm = 1:numEdges
       for nn = 1:numEdges
           predSquaredSum = predSquaredSum + abs(predZmn(mm,nn)).^2;
           unitySquaredSum = unitySquaredSum + abs(unityZmn(mm,nn)).^2;
           if (mm == nn)
               
               predSelfDiff = abs(predZmn(mm,nn) - refZmn(mm,nn));
               predSelfSquaredDiffSum = predSelfSquaredDiffSum + predSelfDiff.^2 ;
   
               %predSelfRelNormPercentError = predSelfDiff.^2 / refSelfZmn.^2;
               %predSelfMSE = predSelfDiff'*predSelfDiff /numSelf;
               
               unitySelfDiff =  abs(unityZmn(mm,nn)  - refZmn(mm,nn));
               unitySelfSquaredDiffSum = unitySelfSquaredDiffSum + unitySelfDiff.^2;
               
               refSelfSquaredSum = refSelfSquaredSum + abs(refZmn(mm,nn)).^2;
               %unitySelfMSE = unitySelfDiff'*unitySelfDiff ./numSelf;
           elseif (singInd(mm,nn))
               triCount = triCount + 1;
               predTriDiff = abs(predZmn(mm,nn) - refZmn(mm,nn));
               predTriSquaredDiffSum = predTriSquaredDiffSum + predTriDiff.^2 ;
               
               unityTriDiff =  abs(unityZmn(mm,nn)  - refZmn(mm,nn));
               unityTriSquaredDiffSum = unityTriSquaredDiffSum + unityTriDiff.^2;
               
               refTriSquaredSum = refTriSquaredSum + abs(refZmn(mm,nn)).^2;
               
           else
               nonSingCount = nonSingCount + 1;
               predNonSingDiff = abs(predZmn(mm,nn) - refZmn(mm,nn));
               predNonSingSquaredDiffSum = predNonSingSquaredDiffSum + predNonSingDiff.^2 ;
               
               unityNonSingDiff =  abs(unityZmn(mm,nn)  - refZmn(mm,nn));
               unityNonSingSquaredDiffSum = unityNonSingSquaredDiffSum + unityNonSingDiff.^2;
               
                refNonSingSquaredSum = refNonSingSquaredSum + abs(refZmn(mm,nn)).^2;
           end
       end % for nn
   end % for mm
   
   predRowRelNormPercentError = zeros(numEdges,1);
   for mm = 1:numEdges
       [~, predRowRelNormPercentError(mm)] = calcError(refZmn(mm,:), predZmn(mm,:));
   end
         
   comp.predNonSingRelNormPercentError = 100 * sqrt(predNonSingSquaredDiffSum / refNonSingSquaredSum);
   comp.unityNonSingRelNormPercentError = 100 * sqrt(unityNonSingSquaredDiffSum / refNonSingSquaredSum);
   comp.predNonSingMSE = predNonSingSquaredDiffSum/nonSingCount;
   comp.unityNonSingMSE = unityNonSingSquaredDiffSum/nonSingCount;
   
   comp.predSelfRelNormPercentError = 100 * sqrt(predSelfSquaredDiffSum / refSelfSquaredSum);
   comp.unitySelfRelNormPercentError = 100 * sqrt(unitySelfSquaredDiffSum / refSelfSquaredSum);
   comp.predSelfMSE = predSelfSquaredDiffSum/numEdges;
   comp.unitySelfMSE = unitySelfSquaredDiffSum/numEdges;
   
   comp.predTriRelNormPercentError = 100 * sqrt(predTriSquaredDiffSum / refTriSquaredSum);
   comp.unityTriRelNormPercentError = 100 * sqrt(unityTriSquaredDiffSum / refTriSquaredSum);
   comp.predTriMSE = predTriSquaredDiffSum/triCount;
   comp.unityTriMSE = unityTriSquaredDiffSum/triCount;

   comp.refFrobNorm = sqrt(refNonSingSquaredSum + refTriSquaredSum + refSelfSquaredSum);
   comp.predFrobNorm = sqrt(predSquaredSum);
   comp.unityFrobNorm = sqrt(unitySquaredSum);
   
   [~, comp.predRelNormPercentError] = calcError(refZmn, predZmn);
    comp.predRowRelNormPercentError = predRowRelNormPercentError;
    
   [~, comp.unityRelNormPercentError] = calcError(refZmn, unityZmn);
   
end
