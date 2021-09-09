function [zMatrices] = predictSolverSetup(Const, Solver_setup, mlmom, includeRealCalc)
    % mlmom : trained model
    % Assume frequencies in solver_setup are the same as mlmom
    % includeReal : sometimes only the imaginary part has large errors
    
    [selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp, singInd] = extractZmnInfo(Const, Solver_setup);
    [zMatrices] = predictExtractedTerms( mlmom, selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp, singInd, includeRealCalc);

end

