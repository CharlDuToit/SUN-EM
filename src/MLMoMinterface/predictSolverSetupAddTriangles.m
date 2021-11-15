function [predictedSetup] = predictSolverSetupAddTriangles(Const, Solver_setup, mlmomAddTriangles, includeError)
    % mlmomAddTriangles : trained model
    % Assume frequencies in solver_setup are the same
    

    numFreq = Solver_setup.frequencies.freq_num;
    message_fc(Const, sprintf('Predicting solver setup '));

    % ======================= EXTRACT OLD ======================= 
    % Extract old terms and properties
    %message_fc(Const, sprintf('  Calculating old solver setup matrices'));
    tic;
    oldCentreDistanceProperties = createCentreDistanceProperties(Solver_setup );
    [oldTerms, oldSingInd] = fillZmnTermsByEdge(Const,Solver_setup, oldCentreDistanceProperties);
    oldEdgeCentreProperties = createEdgeCentreProperties(Solver_setup);
    oldRhoProperties = createRhoProperties(Solver_setup);
    oldEdgeLengths = Solver_setup.rwg_basis_functions_length_m;
    %[~ , ~, numTerms] = size(oldTerms);
    oldZmnTime = toc;
    
        % ======================= ADD TRIANGLES ======================= 
    % Add triangles to solver setup
    message_fc(Const, sprintf('  Creating new solver setup '));
    tic;
    [new_solver_setup, newEdgeLinkOldEdge, newEdgeParallelExternalEdgeLinkOldInternalEdge] = addTriangles(Solver_setup);
    newEdgeLengths = new_solver_setup.rwg_basis_functions_length_m;
    numNewEdges = new_solver_setup.num_mom_basis_functions;
    addTrianglesTime = toc;
    
     % ======================= PROJECT ======================= 
     % Project old terms to new terms
     message_fc(Const, sprintf('  Projecting old solver setup matrices '));
     tic;
     newCentreDistanceProperties = createCentreDistanceProperties(new_solver_setup );
     [newTerms, newSingInd, newLinkOld, groupIndices] = ...
         projectOldSolverSetup(new_solver_setup, Solver_setup,newEdgeLinkOldEdge,...
         newEdgeParallelExternalEdgeLinkOldInternalEdge, oldTerms, oldSingInd, newCentreDistanceProperties, oldCentreDistanceProperties );
     newEdgeCentreProperties = createEdgeCentreProperties(new_solver_setup);
     newRhoProperties = createRhoProperties(new_solver_setup);
     projectTime = toc;
     
     % ======================= ASSIGN  ======================= 
     tic
     message_fc(Const, sprintf('  Assigning clusters '));
     groupIndices = assignGroupClusters(mlmomAddTriangles.groupMeans, groupIndices,newLinkOld, newEdgeCentreProperties,...
         oldEdgeCentreProperties,newEdgeLengths, oldEdgeLengths, newRhoProperties, oldRhoProperties);
     assignTime = toc;
     
     %=============== REFERENCE ZMN (if includeError)============ 

     refZmn = [];
     refZmnTime = 0;
     if (includeError)
         tic;
         quadPts = Const.QUAD_PTS;
         Const.QUAD_PTS = 12;
         [refZMatrices] = FillZMatrixByEdge(Const,new_solver_setup) ;
         refZmn = refZMatrices.values;
         Const.QUAD_PTS = quadPts;
         refZmnTime = toc;
     end
     
     %======================= MULTIPLY WEIGHTS  ======================= 
     tic;
     message_fc(Const, sprintf('  Multiplying weights '));
     
     EMag = 1;
     theta_0 = 0;
     phi_0 = 0;
     y = zeros(numNewEdges, numFreq);
     refX = zeros(numNewEdges, numFreq);
     predX = zeros(numNewEdges, numFreq);
     %erros
     predXError = zeros(numFreq, 3);
     
     predZmn = zeros(numNewEdges, numNewEdges,numFreq );
     projZmn = zeros(numNewEdges, numNewEdges,numFreq );
     %Errors
     groupErrors = cell(numFreq, 2);
     predRelNormPercentError = zeros(numFreq, 3);
     errorSummary = cell(numFreq, 2); 
     predRowRelNormPercentError = zeros(numNewEdges, numFreq);
     %Frob norms
     predFrobNorms = zeros(numFreq, 3);
     projFrobNorms = zeros(numFreq, 3);
     refFrobNorms = zeros(numFreq, 3);
     for f = 1:numFreq
         [groupErrorsReal, projZmnReal, predZmnReal] = applyGroupWeights(mlmomAddTriangles.weightModels{f,1}, groupIndices, real(newTerms(:,:,:,f)), real(refZmn(:,:,f)), includeError);
         [groupErrorsImag, projZmnImag, predZmnImag] = applyGroupWeights(mlmomAddTriangles.weightModels{f,2}, groupIndices, imag(newTerms(:,:,:,f)), imag(refZmn(:,:,f)), includeError);
         projZmn(:, :, f) = projZmnReal + 1i * projZmnImag;
         predZmn(:, :, f) = predZmnReal + 1i * predZmnImag;
         
         y(:,f) = FillVVector(Const, new_solver_setup, EMag,theta_0,phi_0);
         predX(:,f) = predZmn(:,:,f)\y(:,f);
         
         if (includeError)
             refX(:,f) = refZmn(:,:,f)\y(:,f);
             [~,predXError(f,1)] = calcError(real(refX(:,f)),real(predX(:,f)));
             [~,predXError(f,2)] = calcError(imag(refX(:,f)),imag(predX(:,f)));
             [~,predXError(f,3)] = calcError(refX(:,f),predX(:,f));
             
             groupErrors{f,1} = groupErrorsReal;
             groupErrors{f,2} = groupErrorsImag;
             [compReal] = compareZmn(real(refZmn(:, :, f)), real(predZmn(:, :, f)), real(projZmn(:, :, f)), newSingInd);
             [compImag] = compareZmn(imag(refZmn(:, :, f)), imag(predZmn(:, :, f)), imag(projZmn(:, :, f)), newSingInd);
             %[compComplex] = compareZmn(refZmn(:, :, f), predZmn(:, :, f),projZmn(:, :, f), newSingInd);
             %[compComplex] = compareZmn(predictedSetup.refZmn, predictedSetup.predZmn,predictedSetup.projZmn, predictedSetup.newSingInd);
             errorSummary{f,1} = compReal;
             errorSummary{f,2} = compImag;
             
             predRelNormPercentError(f,1) = compReal.predRelNormPercentError;
             predRelNormPercentError(f,2)= compImag.predRelNormPercentError;
             [~, predRelNormPercentError(f,3)] = calcError(refZmn(:, :, f), predZmn(:, :, f));
             
             for mm = 1:numNewEdges
                 [~, predRowRelNormPercentError(mm)] = calcError(refZmn(mm, :, f), predZmn(mm, :, f));
             end
             
             predFrobNorms(f,1) = compReal.predFrobNorm;
             predFrobNorms(f,2) = compImag.predFrobNorm;
             predFrobNorms(f,3) = calcFrobNorm(predZmn(:,:,f));
             
             projFrobNorms(f,1) = compReal.unityFrobNorm;
             projFrobNorms(f,2) = compImag.unityFrobNorm;
             projFrobNorms(f,3) = calcFrobNorm(projZmn(:,:,f));
             
             refFrobNorms(f,1) = compReal.refFrobNorm;
             refFrobNorms(f,2) = compImag.refFrobNorm;
             refFrobNorms(f,3) = calcFrobNorm(refZmn(:,:,f));
         end
     end
     multiplyTime = toc;
     
     %======================= UPDATE STRUCTURE   =======================
     predictedSetup = [];
     
     predictedSetup.y = y;
     predictedSetup.refX = refX;
     predictedSetup.Isol = predX;
     %erros
     predictedSetup.predXError = predXError;
     
     predictedSetup.groupIndices = groupIndices;
     predictedSetup.new_solver_setup = new_solver_setup;
     
     if (includeError)
         predictedSetup.errorSummary = errorSummary;
         predictedSetup.predFrobNorms = predFrobNorms;
         predictedSetup.projFrobNorms = projFrobNorms;
         predictedSetup.refFrobNorms = refFrobNorms;
         predictedSetup.groupErrors = groupErrors;
         predictedSetup.predRelNormPercentError = predRelNormPercentError;
         predictedSetup.predRowRelNormPercentError = predRowRelNormPercentError;
         predictedSetup.refZmn = refZmn;
     end
     
     %matrices
     predictedSetup.projZmn = projZmn;
     predictedSetup.predZmn = predZmn;
     predictedSetup.newSingInd = newSingInd;
       
     predictedSetup.newLinkOld = newLinkOld;
     
     predictedSetup.newProperties = newEdgeCentreProperties;
     predictedSetup.newTerms = newTerms;
     predictedSetup.newSingInd = newSingInd;
     
     predictedSetup.oldProperties = oldEdgeCentreProperties;
     predictedSetup.oldTerms = oldTerms;
     predictedSetup.oldSingInd = oldSingInd;
     
     %scalars
     %constants
     predictedSetup.includeError = includeError;
     predictedSetup.numFreq = numFreq;
     predictedSetup.quadPts = Const.QUAD_PTS;
     %timing
     predictedSetup.oldZmnTime = oldZmnTime;
     predictedSetup.addTrianglesTime = addTrianglesTime;
     predictedSetup.projectTime = projectTime;
     predictedSetup.assignTime =assignTime;
     predictedSetup.refZmnTime = refZmnTime;
     predictedSetup.multiplyTime =multiplyTime;
     predictedSetup.predictTime = oldZmnTime + addTrianglesTime + projectTime + assignTime + multiplyTime;
end

