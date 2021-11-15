function [predictedSetup] = predictSolverSetup(Const, Solver_setup, mlmom, refStruct, includeError)
    % mlmom : trained model
    % Assume frequencies in solver_setup are the same as mlmom
    % includeReal : sometimes only the imaginary part has large errors
    %includeRealCalc = mlmom.includeRealCalc;
    
    numFreq = Solver_setup.frequencies.freq_num;
    if(numFreq ~= mlmom.numFreq)
        message_fc(Const, sprintf('  MLMoM and Solver-setup does not have the same number of frequencies'));
        predictedSetup = [];
        return;
    end
    for f = 1:numFreq
        if (Solver_setup.frequencies.samples(f) ~= mlmom.freqSamples(f))
            message_fc(Const, sprintf('  MLMoM and Solver-setup does not have the same frequency samples'));
            predictedSetup = [];
            return;
        end
    end
    %------------ EXTRACT TERMS ------------
    [termStruct, propStruct, indicesStruct, singInd ] = extractZmnInfo(Const, Solver_setup,  mlmom.constMeshSize);
    termAndPropCalcTime = termStruct.calcTermTime + termStruct.assignAndSwapTime + propStruct.centreDistanceTime + propStruct.edgeCentreTime;
    termStruct.refTwoUniqueZmn = [];
    termStruct.refThreeUniqueZmn =[];
    termStruct.refFourUniqueZmn = [];
    if (includeError)
        tic;
        [termStruct.refTwoUniqueZmn, termStruct.refThreeUniqueZmn, termStruct.refFourUniqueZmn] = extractRefZmn(refStruct.zMatrices.values, singInd );
        %[~ , numTerms] = size(termStruct.fourUniqueTerms);
        refZmnExtractTime = toc;
    end
    
    %------------ PREDICT TERMS ------------
    tic;
    numEdges = Solver_setup.num_mom_basis_functions;
    unityZmn = zeros(numEdges,numEdges,numFreq );
    for f = 1:numFreq
        for mm = 1:numEdges
            for nn= 1:numEdges
                unityZmn(mm,nn,f) = sum(termStruct.allTerms(mm,nn,:,f) );
            end
        end
    end
    unityZmnTime = toc;
    %[predZmn, indicesStruct] = loopPrediction_(mlmom.weightModels, mlmom.clusterStruct, propStruct, termStruct,indicesStruct,unityZmn, mlmom.includeRealCalc, numFreq);
    [predZmn, indicesStruct, assignTime, multiplyTime] = loopPrediction(mlmom.weightModels, mlmom.clusterStruct, propStruct, termStruct,indicesStruct,unityZmn, mlmom.includeRealCalc, numFreq);
    
    %assignAndMultiplyCalcTime = toc;
    %------------ CALCULATE ERROR ------------
    errorSummary = cell(numFreq, 3);
    refZmn = [];
    refY = [];
    
    refX = [];
    unityX = [];
    predX = [];
    
    unityXError = zeros(numFreq, 3);
    predXError = zeros(numFreq, 3);
    if (includeError)
        refZmn = refStruct.zMatrices.values;
        refY = refStruct.yVectors.values;
        refX = refStruct.xVectors.Isol;
        
        unityX = zeros(numEdges, numFreq);
        predX = zeros(numEdges, numFreq);
        
        for f = 1:numFreq
            
            unityX(:,f) = unityZmn(:,:,f)\refY(:,f);
            predX(:,f) = predZmn(:,:,f)\refY(:,f);
            
            [~,unityXError(f,1)] = calcError(real(refX(:,f)),real(unityX(:,f)));
            [~,predXError(f,1)] = calcError(real(refX(:,f)),real(predX(:,f)));
            [~,unityXError(f,2)] = calcError(imag(refX(:,f)),imag(unityX(:,f)));
            [~,predXError(f,2)] = calcError(imag(refX(:,f)),imag(predX(:,f)));
            [~,unityXError(f,3)] = calcError(refX(:,f),unityX(:,f));
            [~,predXError(f,3)] = calcError(refX(:,f),predX(:,f));
            
            [errorSummary{f,1}] = compareZmn(real(refZmn(:, :, f)), real(predZmn(:, :, f)), real(unityZmn(:, :, f)), singInd);
            [errorSummary{f,2}] = compareZmn(imag(refZmn(:, :, f)), imag(predZmn(:, :, f)), imag(unityZmn(:, :, f)), singInd);
            [errorSummary{f,3}] = compareZmn(refZmn(:, :, f), predZmn(:, :, f),unityZmn(:, :, f), singInd);
        end
    end
    
    %------------ UPDATE STRUCT ------------
     predictedSetup = [];
     predictedSetup.name = "mlmom";
     predictedSetup.numSols = numFreq;
     predictedSetup.freqSamples = Solver_setup.frequencies.samples;
     
     predictedSetup.indicesStruct = indicesStruct;
     predictedSetup.propStruct = propStruct;
     predictedSetup.termStruct = termStruct;
     
     if (includeError)
         predictedSetup.predXError = predXError;
         predictedSetup.unityXError = unityXError;
         predictedSetup.errorSummary = errorSummary;
         predictedSetup.refStruct = refStruct;
         predictedSetup.Isol = predX;
         predictedSetup.unityX = unityX;
         predictedSetup.refZmnExtractTime = refZmnExtractTime; 
     end
     predictedSetup.predZmn = predZmn;
     predictedSetup.unityZmn = unityZmn;
     predictedSetup.singInd = singInd;
     
     %timing
     predictedSetup.calcTermTime = termStruct.calcTermTime;
     predictedSetup.assignAndSwapTime = termStruct.assignAndSwapTime;
     predictedSetup.centreDistanceTime = propStruct.centreDistanceTime;
     predictedSetup.edgeCentreTime = propStruct.edgeCentreTime;
     
     predictedSetup.termAndPropCalcTime = termAndPropCalcTime;
     predictedSetup.unityZmnTime = unityZmnTime;
     predictedSetup.assignTime = assignTime;
     predictedSetup.multiplyTime = multiplyTime;
     predictedSetup.predictCalcTime = termAndPropCalcTime + unityZmnTime + assignTime + multiplyTime;
     
     % Write the MLMoM solution to a ASCII str file, so that it can be read
     % again by FEKO (for plotting in POSTFEKO) - only if requested (i.e. if the filename is defined)
     if (includeError && ~isempty(Const.SUNEMmlmomstrfilename))
         writeSolToFile(Const, predictedSetup);
     end%if
end

function [predZmn, indicesStruct, assignTimeTot, multiplyTimeTot] = loopPrediction(weightModels, clusterStruct,propStruct, termStruct,indicesStruct,unityZmn ,includeRealCalc, numFreq)
    predZmn = complex(zeros(termStruct.numEdges, termStruct.numEdges, numFreq));
    hasCalcInd = 0;
    assignTimeTot = 0;
    multiplyTimeTot = 0;
    assignTimeReal = 0;
    multiplyTimeReal = 0;
    for f = 1:numFreq
        if (includeRealCalc)
            %[predZmnReal, indicesStruct] = predictTerms( weightModels{f, 1}, clusterStruct,propStruct, extractTermPart(termStruct,1), indicesStruct, hasCalcInd);
            [predZmnReal, indicesStruct, assignTimeReal, multiplyTimeReal] = predictTerms_temp( weightModels{f, 1}, clusterStruct,propStruct, extractTermPart(termStruct,1, f), indicesStruct, hasCalcInd);
            hasCalcInd = 1;
        else
            predZmnReal = real(unityZmn(:,:,f));
        end
        %[predZmnImag, indicesStruct] = predictTerms( weightModels{f, 2}, clusterStruct,propStruct, extractTermPart(termStruct,0), indicesStruct, hasCalcInd);
        [predZmnImag, indicesStruct, assignTimeImag, multiplyTimeImag] = predictTerms_temp( weightModels{f, 2}, clusterStruct,propStruct, extractTermPart(termStruct,0, f), indicesStruct, hasCalcInd);
        hasCalcInd = 1;
        predZmn(:,:,f) = predZmnReal + 1i* predZmnImag;
        
        assignTimeTot = assignTimeTot + assignTimeImag + assignTimeReal;
        multiplyTimeTot = multiplyTimeTot +  multiplyTimeReal + multiplyTimeImag;
        
    end
end


function termStruct = extractTermPart(termStruct, useReal, f)

    if (useReal)
        termStruct.twoUniqueTerms = real(termStruct.twoUniqueTerms(:,:,f) ) ;
        termStruct.threeUniqueTerms = real(termStruct.threeUniqueTerms(:,:,f) );
        termStruct.fourUniqueTerms = real(termStruct.fourUniqueTerms(:,:,f) );
        %termStruct.refTwoUniqueZmn = real(termStruct.refTwoUniqueZmn(:,f) );
        %termStruct.refThreeUniqueZmn = real(termStruct.refThreeUniqueZmn(:,f) );
        %termStruct.refFourUniqueZmn = real(termStruct.refFourUniqueZmn(:,f) );
    else
        termStruct.twoUniqueTerms = imag(termStruct.twoUniqueTerms(:,:,f) ) ;
        termStruct.threeUniqueTerms = imag(termStruct.threeUniqueTerms(:,:,f) );
        termStruct.fourUniqueTerms = imag(termStruct.fourUniqueTerms(:,:,f) );
        %termStruct.refTwoUniqueZmn = imag(termStruct.refTwoUniqueZmn(:,f) );
        %termStruct.refThreeUniqueZmn = imag(termStruct.refThreeUniqueZmn(:,f) );
        %termStruct.refFourUniqueZmn = imag(termStruct.refFourUniqueZmn(:,f) );
    end
    
end

