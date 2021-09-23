function [comp] = compareZmn(refZmn, predZmn, unityZmn, singInd, freqIndex, useReal)
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
   
   if (useReal)
       refZmn = real(refZmn(:,:,freqIndex));
       predZmn = real(predZmn(:,:,freqIndex));
       unityZmn = real(unityZmn(:,:,freqIndex));
   else
       refZmn = imag(refZmn(:,:,freqIndex));
       predZmn = imag(predZmn(:,:,freqIndex));
       unityZmn = imag(unityZmn(:,:,freqIndex));
   end
   
   for mm = 1:numEdges
       for nn = 1:numEdges
           predSquaredSum = predSquaredSum + predZmn(mm,nn).^2;
           unitySquaredSum = unitySquaredSum + unityZmn(mm,nn).^2;
           if (mm == nn)
               
               predSelfDiff = predZmn(mm,nn) - refZmn(mm,nn);
               predSelfSquaredDiffSum = predSelfSquaredDiffSum + predSelfDiff.^2 ;
   
               %predSelfRelNormPercentError = predSelfDiff.^2 / refSelfZmn.^2;
               %predSelfMSE = predSelfDiff'*predSelfDiff /numSelf;
               
               unitySelfDiff =  unityZmn(mm,nn)  - refZmn(mm,nn);
               unitySelfSquaredDiffSum = unitySelfSquaredDiffSum + unitySelfDiff.^2;
               
               refSelfSquaredSum = refSelfSquaredSum + refZmn(mm,nn).^2;
               %unitySelfMSE = unitySelfDiff'*unitySelfDiff ./numSelf;
           elseif (singInd(mm,nn))
               triCount = triCount + 1;
               predTriDiff = predZmn(mm,nn) - refZmn(mm,nn);
               predTriSquaredDiffSum = predTriSquaredDiffSum + predTriDiff.^2 ;
               
               unityTriDiff =  unityZmn(mm,nn)  - refZmn(mm,nn);
               unityTriSquaredDiffSum = unityTriSquaredDiffSum + unityTriDiff.^2;
               
               refTriSquaredSum = refTriSquaredSum + refZmn(mm,nn).^2;
               
           else
               nonSingCount = nonSingCount + 1;
               predNonSingDiff = predZmn(mm,nn) - refZmn(mm,nn);
               predNonSingSquaredDiffSum = predNonSingSquaredDiffSum + predNonSingDiff.^2 ;
               
               unityNonSingDiff =  unityZmn(mm,nn)  - refZmn(mm,nn);
               unityNonSingSquaredDiffSum = unityNonSingSquaredDiffSum + unityNonSingDiff.^2;
               
                refNonSingSquaredSum = refNonSingSquaredSum + refZmn(mm,nn).^2;
           end
       end % for nn
   end % for mm
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
end