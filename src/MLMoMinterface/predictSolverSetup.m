function [predZmn, unityZmn, singInd] = predictSolverSetup(Const, Solver_setup, mlmom, includeRealCalc, returnUnity)
    % mlmom : trained model
    % Assume frequencies in solver_setup are the same as mlmom
    % includeReal : sometimes only the imaginary part has large errors
    
    [selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp,nonSingEdgeLabels, singInd] = extractZmnInfo(Const, Solver_setup);
    edgeLengths = Solver_setup.rwg_basis_functions_length_m;
    [predZmn, unityZmn] =...
        predictExtractedTerms( mlmom, selfZmnTerms, triZmnTerms, nonSingZmnTerms, nonSingZmnProp, singInd,edgeLengths, includeRealCalc, returnUnity);

end

