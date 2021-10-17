function [predictedSetup] = predictSolverSetup(Const, Solver_setup, mlmom, refZMatrices, returnUnity, includeError)
    % mlmom : trained model
    % Assume frequencies in solver_setup are the same as mlmom
    % includeReal : sometimes only the imaginary part has large errors
    %includeRealCalc = mlmom.includeRealCalc;
    
    numFreq = numel(Solver_setup.frequencies.freq_num);
    %------------ EXTRACT TERMS ------------
    tic;
    [~,~,selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp,nonSingEdgeLabels, singInd] = extractZmnInfo(Const, Solver_setup);
    edgeLengths = Solver_setup.rwg_basis_functions_length_m;
    unityZmnCalcTime = toc;
    
    %------------ PREDICT TERMS ------------
    tic;
    if (includeError)
        returnUnity = 1;
    end
    [predZmn, unityZmn] =...
        predictExtractedTerms( mlmom, selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp, singInd,edgeLengths, returnUnity);
    assignAndMultiplyCalcTime = toc;
    %------------ CALCULATE ERROR ------------
    errorSummary = cell(numFreq, 3);
    refZmn = [];
    if (includeError)
        refZmn = refZMatrices.values;
        for f = 1:numFreq
            [errorSummary{f,1}] = compareZmn(real(refZmn(:, :, f)), real(predZmn(:, :, f)), real(unityZmn(:, :, f)), singInd);
            [errorSummary{f,2}] = compareZmn(imag(refZmn(:, :, f)), imag(predZmn(:, :, f)), imag(unityZmn(:, :, f)), singInd);
            [errorSummary{f,3} ] = compareZmn(refZmn(:, :, f), predZmn(:, :, f),unityZmn(:, :, f), singInd);
        end
    end
    
    %------------ UPDATE STRUCT ------------
     predictedSetup = [];
     
     if (includeError)
         predictedSetup.errorSummary = errorSummary;
         predictedSetup.refZmn = refZmn;
     end
     predictedSetup.predZmn = predZmn;
     predictedSetup.unityZmn = unityZmn;
     predictedSetup.singInd = singInd;
     
     %timing
     predictedSetup.unityZmnCalcTime = unityZmnCalcTime;
     predictedSetup.assignAndMultiplyCalcTime = assignAndMultiplyCalcTime;
     predictedSetup.predictCalcTime = unityZmnCalcTime + assignAndMultiplyCalcTime;
end

