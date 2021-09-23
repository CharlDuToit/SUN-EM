function  plotMLMOM(mlmom, activePlots, gridSize, freqInd, clusterIndex)
    %freqInd : scalar index
    %  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55]
    % [0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
    % [0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0]
    % [0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
    %[0 0 0 0 0 0 0

    
    % 0 or 1 values stored in activePlots(1:NumPlots)
    % 1 Surface Error Real Unity 3 Prop
    % 2 Surface Error Real Pred 3 Prop
    % 3 Surface Error Real Squared Unity 3 Prop
    % 4 Surface Error Real Squared Pred 3 Prop
    % 5 Surface Error Real Log Norm Squared Unity 3 Prop
    % 6 Surface Error Real Log Norm Squared Pred 3 Prop
    
    % 7 Surface Error Imag Unity 3 Prop
    % 8 Surface Error Imag Pred 3 Prop
    % 9 Surface Error Imag Squared Unity 3 Prop
    % 10 Surface Error Imag Squared Pred 3 Prop
    % 11 Surface Error Imag Log Norm Squared Unity 3 Prop
    % 12 Surface Error Imag Log Norm Squared Pred 3 Prop
    
    % 13 Error Squared Sorted NonSing Real Unity, Sorted Pred
    % 14 Error Squared Sorted Tri Real Unity, Sorted Pred
    % 15 Error Squared Sorted Self Real Unity, Sorted Pred
    % 16 Log Norm Squared Sorted NonSing Real Unity, Sorted Pred
    % 17 Log Norm Squared Sorted Tri Real Unity, Sorted Pred
    % 18 Log Norm Squared Sorted Self Real Unity, Sorted Pred
    
    % 19 Error Squared Sorted NonSing Imag Unity, Sorted Pred
    % 20 Error Squared Sorted Tri Imag Unity, Sorted Pred
    % 21 Error Squared Sorted Self Imag Unity, Sorted Pred
    % 22 Log Norm Squared Sorted NonSing Imag Unity, Sorted Pred
    % 23 Log Norm Squared Sorted Tri Imag Unity, Sorted Pred
    % 24 Log Norm Squared Sorted Self Imag Unity, Sorted Pred
       
    % 25 Complex Magnitude Squared Sorted NonSing Ref, Reordered Unity
    % 26 Complex Magnitude Squared Sorted Tri Ref, Reordered Unity
    % 27 Complex Magnitude Squared Sorted Self Ref, Reordered Unity
    % 28 Complex Magnitude Squared Sorted NonSing Ref, Reordered Pred
    % 29 Complex Magnitude Squared Sorted Tri Ref, Reordered Pred
    % 30 Complex Magnitude Squared Sorted Self Ref, Reordered Pred
    
    % 31 Complex Phase Sorted NonSing Ref, Reordered Unity
    % 32 Complex Phase Sorted Tri Ref, Reordered Unity
    % 33 Complex Phase Sorted Self Ref, Reordered Unity
    % 34 Complex Phase Sorted NonSing Ref, Reordered Pred
    % 35 Complex Phase Sorted Tri Ref, Reordered Pred
    % 36 Complex Phase Sorted Self Ref, Reordered Pred
    
    % 37 Real Cluster
    % 38 Imag Cluster
    % 39 Surface Value Real Ref 3 Prop
    % 40 Surface Value Real Pred 3 Prop
    % 41 Surface Value Real Unity 3 Prop
    % 42 Surface Value Imag Ref 3 Prop
    % 43 Surface Value Imag Pred 3 Prop
    % 44 Surface Value Imag Unity 3 Prop
    
    % 45 Surface Error Real Cluster Unity, Pred Length m vs Length n
    % 46 Surface Error Real Cluster Log Norm Squared Unity, Pred Length m vs Length n
    
    % 47 Surface Error Imag Cluster Unity, Pred Length m vs Length n
    % 48 Surface Error Imag Cluster Log Norm Squared Unity, Pred Length m vs Length n
    
    % 49 Surface Value Real Cluster Ref Pred Unity Length m vs Length n
    % 50 Surface Value Imag Cluster Ref Pred Unity Length m vs Length n
    
    % 51 Surface Error Real Cluster Log Norm Squared Pred - Unity Length m vs Length n
    % 52 Surface Error Imag Cluster Log Norm Squared Pred - Unity Length m vs Length n
    prop = mlmom.nonSingZmnProp;
    numPlots = numel(activePlots);
    for k = 1:numPlots
        if (activePlots(k) == 0)
            continue
        end
        switch k
            %real surface
            case 1 % 1 Surface Error Real Unity 3 Prop
                model = mlmom.weightModels{freqInd, 1};
                err = model.unityNonSingZmn - model.refNonSingZmn;
                plotTitle = 'Unity Real Difference Error ';
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 2% 2 Surface Error Real Pred 3 Prop
                model = mlmom.weightModels{freqInd, 1};
                err = model.predNonSingZmn - model.refNonSingZmn;
                plotTitle = 'Prediction Real Difference Error';
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 3% 3 Surface Error Real Squared Unity 3 Prop
                model = mlmom.weightModels{freqInd, 1};
                err = (model.unityNonSingZmn - model.refNonSingZmn).^2;
                plotTitle = 'Unity Real Squared Difference Error';  
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 4 % 4 Surface Error Real Squared Pred 3 Prop
                model = mlmom.weightModels{freqInd, 1};
                err = (model.predNonSingZmn - model.refNonSingZmn).^2;
                plotTitle = 'Prediction Real Squared Difference Error';
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);       
            case 5% 5 Surface Error Real Log Norm Squared Unity 3 Prop
                model = mlmom.weightModels{freqInd, 1};
                err = 0.5*log((model.unityNonSingZmn - model.refNonSingZmn).^2./ model.refNonSingZmn.^2);
                plotTitle = 'Unity Real Log of Sqrt of Normalised Squared Difference Error ';
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 6% 6 Surface Error Real Log Norm Squared Pred 3 Prop
                model = mlmom.weightModels{freqInd, 1};
                err = log((model.predNonSingZmn - model.refNonSingZmn).^2./ model.refNonSingZmn.^2);
                plotTitle = 'Unity Real Log of Sqrt of Normalised Squared Difference Error ';
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            % Imaginary surface
            case 7% 7 Surface Error Imag Unity 3 Prop
                model = mlmom.weightModels{freqInd, 2};
                err = model.unityNonSingZmn - model.refNonSingZmn;
                plotTitle = 'Unity Imag Difference Error ';
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 8% 8 Surface Error Imag Pred 3 Prop
                model = mlmom.weightModels{freqInd, 2};
                err = model.predNonSingZmn - model.refNonSingZmn;
                plotTitle = 'Prediction Imag Difference Error';
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 9% 9 Surface Error Imag Squared Unity 3 Prop
                model = mlmom.weightModels{freqInd, 2};
                err = (model.unityNonSingZmn - model.refNonSingZmn).^2;
                plotTitle = 'Unity Imag Squared Difference Error';
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 10% 10 Surface Error Imag Squared Pred 3 Prop
                model = mlmom.weightModels{freqInd, 2};
                err = (model.predNonSingZmn - model.refNonSingZmn).^2;
                plotTitle = 'Prediction Imag Squared Difference Error';
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 11 % 11 Surface Error Imag Log Norm Squared Unity 3 Prop
                model = mlmom.weightModels{freqInd, 2};
                err = 0.5*log((model.unityNonSingZmn - model.refNonSingZmn).^2 ./ model.refNonSingZmn.^2);
                plotTitle = 'Unity Imag Log of Sqrt of Normalised Squared Difference Error ';
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 12% 12 Surface Error Imag Log Norm Squared Pred 3 Prop
                model = mlmom.weightModels{freqInd, 2};
                err = 0.5*log((model.predNonSingZmn - model.refNonSingZmn).^2 ./ model.refNonSingZmn.^2);
                plotTitle = 'Prediction Imag Log of Sqrt of Normalised Squared Difference Error ';
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 13 %Error Squared Sorted NonSing Real Unity, Sorted Pred
                model = mlmom.weightModels{freqInd, 1};
                predErr = (model.predNonSingZmn -model.refNonSingZmn).^2;
                unityErr = (model.unityNonSingZmn -model.refNonSingZmn).^2;
                plotTitle = 'NonSing Real Squared Difference Error';
                plotBothSortedError( predErr, unityErr, plotTitle);
                
            case 14% 14 Error Squared Sorted Tri Real Unity, Sorted Pred
                model = mlmom.weightModels{freqInd, 1};
                predErr = (model.predTriZmn -model.refTriZmn).^2;
                unityErr = (model.unityTriZmn -model.refTriZmn).^2;
                plotTitle = 'Tri Real Squared Difference Error';
                plotBothSortedError( predErr, unityErr, plotTitle);
            case 15% 15 Error Squared Sorted Self Real Unity, Sorted Pred
                model = mlmom.weightModels{freqInd, 1};
                predErr = (model.predSelfZmn -model.refSelfZmn).^2;
                unityErr = (model.unitySelfZmn -model.refSelfZmn).^2;
                plotTitle = 'Self Real Squared Difference Error';
                plotBothSortedError( predErr, unityErr, plotTitle);
            case 16% 16 Log Norm Squared Sorted NonSing Real Unity, Sorted Pred
                model = mlmom.weightModels{freqInd, 1};
                predErr = 0.5*log((model.predNonSingZmn -model.refNonSingZmn).^2 ./model.refNonSingZmn.^2);
                unityErr = 0.5*log((model.unityNonSingZmn -model.refNonSingZmn).^2 ./model.refNonSingZmn.^2);
                plotTitle = 'NonSing Real Log of Sqrt of Normalised Squared Difference Error';
                plotBothSortedError( predErr, unityErr, plotTitle);
            case 17% 17 Log Norm Squared Sorted Tri Real Unity, Sorted Pred
                model = mlmom.weightModels{freqInd, 1};
                predErr = 0.5*log((model.predTriZmn -model.refTriZmn).^2 ./model.refTriZmn.^2);
                unityErr = 0.5*log((model.unityTriZmn -model.refTriZmn).^2 ./model.refTriZmn.^2);
                plotTitle = 'Tri Real Log of Sqrt of Normalised Squared Difference Error';
                plotBothSortedError( predErr, unityErr, plotTitle);
            case 18% 18 Log Norm Squared Sorted Self Real Unity, Sorted Pred
                model = mlmom.weightModels{freqInd, 1};
                predErr = 0.5*log((model.predSelfZmn -model.refSelfZmn).^2 ./model.refSelfZmn.^2);
                unityErr = 0.5*log((model.unitySelfZmn -model.refSelfZmn).^2 ./model.refSelfZmn.^2);
                plotTitle = 'Self Real Log of Sqrt of Normalised Squared Difference Error';
                plotBothSortedError( predErr, unityErr, plotTitle);
            case 19% 19 Error Squared Sorted NonSing Imag Unity, Sorted Pred
                model = mlmom.weightModels{freqInd, 2};
                predErr = (model.predNonSingZmn -model.refNonSingZmn).^2;
                unityErr = (model.unityNonSingZmn -model.refNonSingZmn).^2;
                plotTitle = 'NonSing Imag Squared Difference Error';
                plotBothSortedError( predErr, unityErr, plotTitle);
            case 20% 20 Error Squared Sorted Tri Imag Unity, Sorted Pred
                model = mlmom.weightModels{freqInd, 2};
                predErr = (model.predTriZmn -model.refTriZmn).^2;
                unityErr = (model.unityTriZmn -model.refTriZmn).^2;
                plotTitle = 'Tri Imag Squared Difference Error';
                plotBothSortedError( predErr, unityErr, plotTitle);
            case 21% 21 Error Squared Sorted Self Imag Unity, Sorted Pred
                model = mlmom.weightModels{freqInd, 2};
                predErr = (model.predSelfZmn -model.refSelfZmn).^2;
                unityErr = (model.unitySelfZmn -model.refSelfZmn).^2;
                plotTitle = 'Self Imag Squared Difference Error';
                plotBothSortedError( predErr, unityErr, plotTitle);
            case 22% 22 Log Norm Squared Sorted NonSing Imag Unity, Sorted Pred
                model = mlmom.weightModels{freqInd, 2};
                predErr = 0.5*log((model.predNonSingZmn -model.refNonSingZmn).^2 ./model.refNonSingZmn.^2);
                unityErr = 0.5*log((model.unityNonSingZmn -model.refNonSingZmn).^2 ./model.refNonSingZmn.^2);
                plotTitle = 'NonSing Imag Log of Sqrt of Normalised Squared Difference Error';
                plotBothSortedError( predErr, unityErr, plotTitle);
            case 23% 23 Log Norm Squared Sorted Tri Imag Unity, Sorted Pred
                 model = mlmom.weightModels{freqInd, 2};
                predErr = 0.5*log((model.predTriZmn -model.refTriZmn).^2 ./model.refTriZmn.^2);
                unityErr = 0.5*log((model.unityTriZmn -model.refTriZmn).^2 ./model.refTriZmn.^2);
                plotTitle = 'Tri Imag Log of Sqrt of Normalised Squared Difference Error';
                plotBothSortedError( predErr, unityErr, plotTitle);
            case 24% 24 Log Norm Squared Sorted Self Imag Unity, Sorted Pred
                model = mlmom.weightModels{freqInd, 2};
                predErr = 0.5*log((model.predSelfZmn -model.refSelfZmn).^2 ./model.refSelfZmn.^2);
                unityErr = 0.5*log((model.unitySelfZmn -model.refSelfZmn).^2 ./model.refSelfZmn.^2);
                plotTitle = 'Self Imag Log of Sqrt of Normalised Squared Difference Error';
                plotBothSortedError( predErr, unityErr, plotTitle);
            case 25% 25 Complex Magnitude Squared Sorted NonSing Ref, Reordered Unity
                refVal = mlmom.weightModels{freqInd, 1}.refNonSingZmn.^2;
                refVal = refVal + mlmom.weightModels{freqInd, 2}.refNonSingZmn.^2;
                unityVal = mlmom.weightModels{freqInd, 1}.unityNonSingZmn.^2;
                unityVal = unityVal + mlmom.weightModels{freqInd, 2}.unityNonSingZmn.^2;
                plotTitle = 'NonSing Magnitude Squared';
                plotSortedRefReorderedUnity(refVal, unityVal, plotTitle);
            case 26% 26 Complex Magnitude Squared Sorted Tri Ref, Reordered Unity
                refVal = mlmom.weightModels{freqInd, 1}.refTriZmn.^2;
                refVal = refVal + mlmom.weightModels{freqInd, 2}.refTriZmn.^2;
                unityVal = mlmom.weightModels{freqInd, 1}.unityTriZmn.^2;
                unityVal = unityVal + mlmom.weightModels{freqInd, 2}.unityTriZmn.^2;
                plotTitle = 'Tri Magnitude Squared';
                plotSortedRefReorderedUnity(refVal, unityVal, plotTitle);
            case 27% 27 Complex Magnitude Squared Sorted Self Ref, Reordered Unity
                refVal = mlmom.weightModels{freqInd, 1}.refSelfZmn.^2;
                refVal = refVal + mlmom.weightModels{freqInd, 2}.refSelfZmn.^2;
                unityVal = mlmom.weightModels{freqInd, 1}.unitySelfZmn.^2;
                unityVal = unityVal + mlmom.weightModels{freqInd, 2}.unitySelfZmn.^2;
                plotTitle = 'Self Magnitude Squared';
                plotSortedRefReorderedUnity(refVal, unityVal, plotTitle);
            case 28% 28 Complex Magnitude Squared Sorted NonSing Ref, Reordered Pred
                refVal = mlmom.weightModels{freqInd, 1}.refNonSingZmn.^2;
                refVal = refVal + mlmom.weightModels{freqInd, 2}.refNonSingZmn.^2;
                predVal = mlmom.weightModels{freqInd, 1}.predNonSingZmn.^2;
                predVal = predVal + mlmom.weightModels{freqInd, 2}.predNonSingZmn.^2;
                plotTitle = 'NonSing Magnitude Squared';
                plotSortedRefReorderedPred( refVal, predVal, plotTitle);
            case 29% 29 Complex Magnitude Squared Sorted Tri Ref, Reordered Pred
                refVal = mlmom.weightModels{freqInd, 1}.refTriZmn.^2;
                refVal = refVal + mlmom.weightModels{freqInd, 2}.refTriZmn.^2;
                predVal = mlmom.weightModels{freqInd, 1}.predTriZmn.^2;
                predVal = predVal + mlmom.weightModels{freqInd, 2}.predTriZmn.^2;
                plotTitle = 'Tri Magnitude Squared';
                plotSortedRefReorderedPred( refVal, predVal, plotTitle);
            case 30% 30 Complex Magnitude Squared Sorted Self Ref, Reordered Pred
                refVal = mlmom.weightModels{freqInd, 1}.refSelfZmn.^2;
                refVal = refVal + mlmom.weightModels{freqInd, 2}.refSelfZmn.^2;
                predVal = mlmom.weightModels{freqInd, 1}.predSelfZmn.^2;
                predVal = predVal + mlmom.weightModels{freqInd, 2}.predSelfZmn.^2;
                plotTitle = 'Self Magnitude Squared';
                plotSortedRefReorderedPred( refVal, predVal, plotTitle);
            case 31 % 31 Complex Phase Sorted NonSing Ref, Reordered Unity
                refVal = mlmom.weightModels{freqInd, 1}.refNonSingZmn;
                refVal = refVal + 1i *mlmom.weightModels{freqInd, 2}.refNonSingZmn;
                refVal = angle(refVal);
                unityVal = mlmom.weightModels{freqInd, 1}.unityNonSingZmn;
                unityVal = unityVal + 1i* mlmom.weightModels{freqInd, 2}.unityNonSingZmn;
                unityVal = angle(unityVal);
                plotTitle = 'NonSing Phase';
                plotSortedRefReorderedUnity(refVal, unityVal, plotTitle);
            case 32% 32 Complex Phase Sorted Tri Ref, Reordered Unity
                refVal = mlmom.weightModels{freqInd, 1}.refTriZmn;
                refVal = refVal + 1i *mlmom.weightModels{freqInd, 2}.refTriZmn;
                refVal = angle(refVal);
                unityVal = mlmom.weightModels{freqInd, 1}.unityTriZmn;
                unityVal = unityVal + 1i* mlmom.weightModels{freqInd, 2}.unityTriZmn;
                unityVal = angle(unityVal);
                plotTitle = 'Tri Phase';
                plotSortedRefReorderedUnity(refVal, unityVal, plotTitle);
            case 33% 33 Complex Phase Sorted Self Ref, Reordered Unity
                refVal = mlmom.weightModels{freqInd, 1}.refSelfZmn;
                refVal = refVal + 1i *mlmom.weightModels{freqInd, 2}.refSelfZmn;
                refVal = angle(refVal);
                unityVal = mlmom.weightModels{freqInd, 1}.unitySelfZmn;
                unityVal = unityVal + 1i* mlmom.weightModels{freqInd, 2}.unitySelfZmn;
                unityVal = angle(unityVal);
                plotTitle = 'Self Phase';
                plotSortedRefReorderedUnity(refVal, unityVal, plotTitle);
            case 34 % 34 Complex Phase Sorted NonSing Ref, Reordered Pred
                refVal = mlmom.weightModels{freqInd, 1}.refNonSingZmn;
                refVal = refVal + 1i *mlmom.weightModels{freqInd, 2}.refNonSingZmn;
                refVal = angle(refVal);
                predVal = mlmom.weightModels{freqInd, 1}.predNonSingZmn;
                predVal = predVal + 1i* mlmom.weightModels{freqInd, 2}.predNonSingZmn;
                predVal = angle(predVal);
                plotTitle = 'NonSing Phase';
                plotSortedRefReorderedPred(refVal, predVal, plotTitle);
            case 35% 35 Complex Phase Sorted Tri Ref, Reordered Pred
                refVal = mlmom.weightModels{freqInd, 1}.refTriZmn;
                refVal = refVal + 1i *mlmom.weightModels{freqInd, 2}.refTriZmn;
                refVal = angle(refVal);
                predVal = mlmom.weightModels{freqInd, 1}.predTriZmn;
                predVal = predVal + 1i* mlmom.weightModels{freqInd, 2}.predTriZmn;
                predVal = angle(predVal);
                plotTitle = 'Tri Phase';
                plotSortedRefReorderedPred(refVal, predVal, plotTitle);
            case 36% 36 Complex Phase Sorted Self Ref, Reordered Pred
                refVal = mlmom.weightModels{freqInd, 1}.refSelfZmn;
                refVal = refVal + 1i *mlmom.weightModels{freqInd, 2}.refSelfZmn;
                refVal = angle(refVal);
                predVal = mlmom.weightModels{freqInd, 1}.predSelfZmn;
                predVal = predVal + 1i* mlmom.weightModels{freqInd, 2}.predSelfZmn;
                predVal = angle(predVal);
                plotTitle = 'Self Phase';
                plotSortedRefReorderedPred(refVal, predVal, plotTitle);
            case 37
                ind = find(mlmom.clusterInd == clusterIndex);
                refVal = mlmom.weightModels{freqInd, 1}.refNonSingZmn(ind);
                predVal = mlmom.weightModels{freqInd, 1}.predNonSingZmn(ind);
                unityVal = mlmom.weightModels{freqInd, 1}.unityNonSingZmn(ind);
                plotTitle = sprintf('Real Cluster %d ', clusterIndex);
                plotCluster(refVal, predVal, unityVal, plotTitle);
            case 38
                ind = find(mlmom.clusterInd == clusterIndex);
                refVal = mlmom.weightModels{freqInd, 2}.refNonSingZmn(ind);
                predVal = mlmom.weightModels{freqInd, 2}.predNonSingZmn(ind);
                unityVal = mlmom.weightModels{freqInd, 2}.unityNonSingZmn(ind);
                prop = mlmom.nonSingZmnProp(ind); 
                [sortedProp, ind] = sort(prop(:,1));
                refVal = refVal(ind);
                predVal = predVal(ind);
                unityVal = unityVal(ind);
                plotTitle = sprintf('Imag Cluster %d ', clusterIndex);
                plotCluster(refVal, predVal, unityVal, plotTitle);
            case 39 % 39 Surface Real Ref 3 Prop
                ref = mlmom.weightModels{1, 1}.refNonSingZmn;
                plotTitle = 'Real Reference Val';
                plotSurfaceError3Prop(prop, ref, gridSize, plotTitle);
            case 40  %40 Surface Real Pred 3 Prop
                ref = mlmom.weightModels{1, 1}.refNonSingZmn;
                plotTitle = 'Real Reference Val';
                plotSurfaceError3Prop(prop, ref, gridSize, plotTitle);
            case 41 %40 Surface Real Unity 3 Prop
                unity = mlmom.weightModels{1, 1}.unityNonSingZmn;
                plotTitle = 'Real Unity Val';
                plotSurfaceError3Prop(prop, unity, gridSize, plotTitle);
            case 42 % 44 Surface Imag Ref 3 Prop
                ref = mlmom.weightModels{1, 2}.refNonSingZmn;
                plotTitle = 'Imag Reference Val';
                plotSurfaceError3Prop(prop, ref, gridSize, plotTitle);
            case 43% 43 Surface Imag Pred 3 Prop
                pred = mlmom.weightModels{1, 2}.predNonSingZmn;
                plotTitle = 'Imag Prediction Val';
                plotSurfaceError3Prop(prop, pred, gridSize, plotTitle);
            case 44% 44 Surface Imag Unity 3 Prop
                unity = mlmom.weightModels{1, 2}.unityNonSingZmn;
                plotTitle = 'Imag Unity Val';
                plotSurfaceError3Prop(prop, unity, gridSize, plotTitle);
            case 45
            case 46
            case 47 % 47 Surface Error Imag Cluster Unity, Pred Length m vs Length n
                model = mlmom.weightModels{freqInd, 2};
                ind = find(mlmom.clusterInd == clusterIndex);
                refVal = model.refNonSingZmn(ind);
                predVal = model.predNonSingZmn(ind);
                unityVal =model.unityNonSingZmn(ind);
                predError = predVal - refVal;
                unityError = unityVal - refVal;
                edgeLabels = mlmom.nonSingEdgeLabels(ind,:);
                [numPoints, ~] = size(edgeLabels);
                lengths = zeros(numPoints, 2);
                lengths(: , 1) = mlmom.edgeLengths(edgeLabels(:, 1));
                lengths(: , 2) = mlmom.edgeLengths(edgeLabels(:, 2));
                plotTitle = sprintf('Imag Cluster Error %d ', clusterIndex);
                %plotSurfaceClusterValue( refVal, predVal, unityVal, lengths ,gridSize, plotTitle);
                plotSurfaceClusterError( predError, unityError, lengths ,3, plotTitle);
            case 48 % 48 Surface Error Imag Cluster Log Norm Squared Unity, Pred Length m vs Length n
                model = mlmom.weightModels{freqInd, 2};
                ind = find(mlmom.clusterInd == clusterIndex);
                refVal = model.refNonSingZmn(ind);
                predVal = model.predNonSingZmn(ind);
                unityVal =model.unityNonSingZmn(ind);
                predError = 0.5*log((predVal - refVal).^2 ./ refVal.^2);
                unityError = 0.5*log((unityVal - refVal).^2 ./ refVal.^2 );
                edgeLabels = mlmom.nonSingEdgeLabels(ind,:);
                [numPoints, ~] = size(edgeLabels);
                lengths = zeros(numPoints, 2);
                lengths(: , 1) = mlmom.edgeLengths(edgeLabels(:, 1));
                lengths(: , 2) = mlmom.edgeLengths(edgeLabels(:, 2));
                
                %[weights, pred ]= MLR(lengths,unityError, 1e-46, 1, 0);
                plotTitle = sprintf('Imag Cluster %d Log of Sqrt of Normalised Squared Difference Error  ', clusterIndex);
                plotSurfaceClusterError( predError, unityError, lengths ,3, plotTitle);
            case 49 %49 Surface Value Real Cluster Ref Pred Unity Length m vs Length n
                model = mlmom.weightModels{freqInd, 1};
                ind = find(mlmom.clusterInd == clusterIndex);
                refVal = model.refNonSingZmn(ind);
                predVal = model.predNonSingZmn(ind);
                unityVal =model.unityNonSingZmn(ind);
                edgeLabels = mlmom.nonSingEdgeLabels(ind,:);
                [numPoints, ~] = size(edgeLabels);
                lengths = zeros(numPoints, 2);
                lengths(: , 1) = mlmom.edgeLengths(edgeLabels(:, 1));
                lengths(: , 2) = mlmom.edgeLengths(edgeLabels(:, 2));
                plotTitle = sprintf('Real Cluster %d ', clusterIndex);
                %plotSurfaceClusterValue( refVal, predVal, unityVal, lengths ,gridSize, plotTitle);
                plotSurfaceClusterValue( refVal, predVal, unityVal, lengths ,3, plotTitle);
            case 50 %50 Surface Value Imag Cluster Ref Pred Unity Length m vs Length n
                model = mlmom.weightModels{freqInd, 2};
                ind = find(mlmom.clusterInd == clusterIndex);
                refVal = model.refNonSingZmn(ind);
                predVal = model.predNonSingZmn(ind);
                unityVal =model.unityNonSingZmn(ind);
                edgeLabels = mlmom.nonSingEdgeLabels(ind,:);
                [numPoints, ~] = size(edgeLabels);
                lengths = zeros(numPoints, 2);
                lengths(: , 1) = mlmom.edgeLengths(edgeLabels(:, 1));
                lengths(: , 2) = mlmom.edgeLengths(edgeLabels(:, 2));
                
                
                plotTitle = sprintf('Imag Cluster %d ', clusterIndex);
                %plotSurfaceClusterValue( refVal, predVal, unityVal, lengths ,gridSize, plotTitle);
                plotSurfaceClusterValue( refVal, predVal, unityVal, lengths ,3, plotTitle);
            case 51  %51 Surface Error Real Cluster Log Norm Squared Pred - Unity Length m vs Length n
                model = mlmom.weightModels{freqInd, 1};
                ind = find(mlmom.clusterInd == clusterIndex);
                refVal = model.refNonSingZmn(ind);
                predVal = model.predNonSingZmn(ind);
                unityVal =model.unityNonSingZmn(ind);
                predError = 0.5*log((predVal - refVal).^2 ./ refVal.^2);
                unityError = 0.5*log((unityVal - refVal).^2 ./ refVal.^2 );
                edgeLabels = mlmom.nonSingEdgeLabels(ind,:);
                predMinusUnityError = predError - unityError; %we want these values to be negative
                [numPoints, ~] = size(edgeLabels);
                lengths = zeros(numPoints, 2);
                lengths(: , 1) = mlmom.edgeLengths(edgeLabels(:, 1));
                lengths(: , 2) = mlmom.edgeLengths(edgeLabels(:, 2));
                plotTitle = sprintf('Imag Cluster %d Log of Sqrt of Normalised Squared Difference Error  ', clusterIndex);
                plotSurfaceClusterErrorDifference( predMinusUnityError, unityError, lengths ,3, plotTitle);
            case 52 %52 Surface Error Imag Cluster Log Norm Squared Pred - Unity Length m vs Length n
                model = mlmom.weightModels{freqInd, 2};
                ind = find(mlmom.clusterInd == clusterIndex);
                refVal = model.refNonSingZmn(ind);
                predVal = model.predNonSingZmn(ind);
                unityVal =model.unityNonSingZmn(ind);
                predError = 0.5*log((predVal - refVal).^2 ./ refVal.^2);
                unityError = 0.5*log((unityVal - refVal).^2 ./ refVal.^2 );
                edgeLabels = mlmom.nonSingEdgeLabels(ind,:);
                predMinusUnityError = predError - unityError; %we want these values to be negative
                [numPoints, ~] = size(edgeLabels);
                lengths = zeros(numPoints, 2);
                lengths(: , 1) = mlmom.edgeLengths(edgeLabels(:, 1));
                lengths(: , 2) = mlmom.edgeLengths(edgeLabels(:, 2));

                plotTitle = sprintf('Imag Cluster %d Log of Sqrt of Normalised Squared Difference Error  ', clusterIndex);
                plotSurfaceClusterErrorDifference( predMinusUnityError, unityError, lengths ,3, plotTitle);
                
                %[weights, predUnityError ]= MLR(lengths,unityError, 1e-46, 1, 0);
                %plotSurfaceClusterErrorDifference( unityError, predUnityError, lengths ,3, plotTitle);
                
                %[weights, predpredMinusUnityError ]= MLR(lengths,predMinusUnityError, 1e-46, 1, 0);
                %plotSurfaceClusterErrorDifference( predMinusUnityError, predpredMinusUnityError, lengths ,3, plotTitle);
                
            case 53
                
            case 54
            otherwise
                
        end % switch
    end

end



function  plotSurfaceError3Prop(prop, err, gridSize, plotTitle)

    distLin = linspace(min(prop(:,1)) , max(prop(:,1)), gridSize);
    dirDotDirLin = linspace(min(prop(:,2)) , max(prop(:,2)), gridSize);
    dirDotDispLin = linspace(min(prop(:,3)) , max(prop(:,3)), gridSize);
    markerSize = 330000/numel(prop(:,1));
    subplot(2,2,1);
    % x = dist, y = dir dot dir   
    [X , Y] = meshgrid(distLin, dirDotDirLin);
    Z = griddata(prop(:,1) , prop(:, 2) , err, X , Y, 'cubic');  
    mesh(X,Y,Z);   
    axis tight; hold on;
    plot3(prop(:,1) , prop(:, 2) , err , '.', 'MarkerSize', markerSize);
    title(plotTitle);
    xlabel('Distance');
    ylabel('Edge direction dot edge direction');
    %zlim([-0.01 0.01]);
    
    subplot(2,2,2);
    % x = dist, y = dir dot disp
    [X , Y] = meshgrid(distLin, dirDotDispLin);
    Z = griddata(prop(:,1) , prop(:, 3) , err, X , Y, 'cubic');  
    mesh(X,Y,Z);   
    axis tight;
    hold on;
    plot3(prop(:,1) , prop(:, 3) , err , '.', 'MarkerSize', markerSize);
    title(plotTitle);
    xlabel('Distance');
    ylabel('Edge direction dot displacement');
    
      subplot(2,2,3);
    % x = dir dot dir, y = dir dot disp
    [X , Y] = meshgrid(distLin, dirDotDispLin);
    Z = griddata(prop(:,2) , prop(:, 3) , err, X , Y, 'cubic');  
    mesh(X,Y,Z);   
    axis tight;
    hold on;
    plot3(prop(:,2) , prop(:, 3) , err , '.', 'MarkerSize', markerSize);
    title(plotTitle);
    xlabel('Edge direction dot edge direction');
    ylabel('Edge direction dot displacement');
end

function  plotSurfaceClusterValue( refVal, predVal, unityVal, lengths ,gridSize, plotTitle)

    lengthMLin = linspace(min(lengths(:,1)) , max(lengths(:,1)), gridSize);
    lengthNLin = linspace(min(lengths(:,2)) , max(lengths(:,2)), gridSize);
    %markerSize = 33000/numel(refVal);
    markerSize = 15;
    subplot(1,1,1);
    % x = dist, y = dir dot dir   
    [X , Y] = meshgrid(lengthMLin, lengthNLin);
    %ref
     refZ = griddata(lengths(:,1) , lengths(:,2) , refVal, X , Y, 'cubic');  
    mesh(X,Y,refZ);   
    axis tight; hold on;
    plot3(lengths(:,1) , lengths(:,2) , refVal , 'g.', 'MarkerSize', markerSize);
    %pred
     predZ = griddata(lengths(:,1) , lengths(:,2) , predVal, X , Y, 'cubic');  
     mesh(X,Y,predZ);   
     plot3(lengths(:,1) , lengths(:,2) , predVal , 'r.', 'MarkerSize', markerSize);
%     %unity
     unityZ = griddata(lengths(:,1) , lengths(:,2) , unityVal, X , Y, 'cubic');  
     mesh(X,Y,unityZ);   
     plot3(lengths(:,1) , lengths(:,2) , unityVal , 'b.', 'MarkerSize', markerSize);
    
    title(plotTitle);
    xlabel('Edge m length');
    ylabel('Edge n length');
    zlabel('Value');
    legend('Reference','Prediction', 'Unity');
    %zlim([-0.01 0.01]);
    
end

function  plotSurfaceClusterError( predError, unityError, lengths ,gridSize, plotTitle)
    lengthMLin = linspace(min(lengths(:,1)) , max(lengths(:,1)), gridSize);
    lengthNLin = linspace(min(lengths(:,2)) , max(lengths(:,2)), gridSize);
    %markerSize = 33000/numel(refVal);
    markerSize = 15;
    subplot(1,2,1);
    [X , Y] = meshgrid(lengthMLin, lengthNLin);
    %pred
    pred = griddata(lengths(:,1) , lengths(:,2) , predError, X , Y, 'cubic');  
    mesh(X,Y,pred);   
    axis tight; hold on;
    plot3(lengths(:,1) , lengths(:,2) , predError , 'r.', 'MarkerSize', markerSize);
    title(strcat('Prediction',plotTitle));
    xlabel('Edge m length');
    ylabel('Edge n length');
    zlabel('Error');
    subplot(1,2,2);
    % unity
    unity = griddata(lengths(:,1) , lengths(:,2) , predError, X , Y, 'cubic');  
    mesh(X,Y,unity);   
    axis tight; hold on;
    plot3(lengths(:,1) , lengths(:,2) , unityError , 'b.', 'MarkerSize', markerSize);
    title(strcat('Unity',plotTitle));
    xlabel('Edge m length');
    ylabel('Edge n length');
    zlabel('Error');
end

function  plotSurfaceClusterErrorDifference( predMinusUnityError, unityError, lengths ,gridSize, plotTitle)
    lengthMLin = linspace(min(lengths(:,1)) , max(lengths(:,1)), gridSize);
    lengthNLin = linspace(min(lengths(:,2)) , max(lengths(:,2)), gridSize);
    %markerSize = 33000/numel(refVal);
    markerSize = 15;
    subplot(1,2,1);
    [X , Y] = meshgrid(lengthMLin, lengthNLin);
    %pred - Unit
    diff = griddata(lengths(:,1) , lengths(:,2) , predMinusUnityError, X , Y, 'cubic');  
    mesh(X,Y,diff);   
    axis tight; hold on;
    plot3(lengths(:,1) , lengths(:,2) , predMinusUnityError , 'k.', 'MarkerSize', markerSize);
    title(strcat('Prediction - Unity,',plotTitle));
    xlabel('Edge m length');
    ylabel('Edge n length');
    zlabel('Error Difference');
    subplot(1,2,2);
    % unity
    unity = griddata(lengths(:,1) , lengths(:,2) , predMinusUnityError, X , Y, 'cubic');  
    mesh(X,Y,unity);   
    axis tight; hold on;
    plot3(lengths(:,1) , lengths(:,2) , unityError , 'b.', 'MarkerSize', markerSize);
    title(strcat('Unity',plotTitle));
    xlabel('Edge m length');
    ylabel('Edge n length');
    zlabel('Unity Error');
end

function plotCluster(refVal, predVal, unityVal, plotTitle)
     figure;
     markerSize = 800/numel(refVal);
     plot(refVal, 'g.','MarkerSize', markerSize);
     hold on
     plot(predVal, 'r.','MarkerSize', markerSize);
     hold on
     plot(unityVal, 'b.','MarkerSize', markerSize);
     title(plotTitle);
     xlabel('Matrix element');
     ylabel('Value');
     legend('Reference','Prediction', 'Unity');
     hold off;
end

function  plotBothSortedError( predErr, unityErr, plotTitle)
    [sortedPredErr, ~] = sort(predErr);
    [sortedUnityErr, ~] = sort(unityErr);
     figure;
     markerSize = 5000/numel(unityErr);
     plot(sortedPredErr, 'r.','MarkerSize', markerSize);
     hold on
     plot(sortedUnityErr, 'b.','MarkerSize', markerSize);
     title(plotTitle);
     xlabel('Matrix element');
     legend('Sorted Prediction', 'Sorted Unity');
     hold off;
end

function  plotSortedRefReorderedUnity( refVal, unityVal, plotTitle)
    [sortedRefVal, ind] = sort(refVal);
    reorderedUnityVal = unityVal(ind);
    %[sortedUnityErr, ~] = sort(unityVal);
     figure;
     markerSize = 5000/numel(unityVal);
     plot(sortedRefVal, 'g.','MarkerSize', markerSize);
     hold on
     plot(reorderedUnityVal, 'b.','MarkerSize', markerSize);
     title(plotTitle);
     xlabel('Matrix element');
     legend('Sorted Reference', 'Reordered Unity');
     hold off;
end

function  plotSortedRefReorderedPred( refVal, predVal, plotTitle)
    [sortedRefVal, ind] = sort(refVal);
    reorderedPredVal = predVal(ind);
    %[sortedUnityErr, ~] = sort(unityVal);
     figure;
     markerSize = 5000/numel(predVal);
     plot(sortedRefVal, 'g.','MarkerSize', markerSize);
     hold on
     plot(reorderedPredVal, 'r.','MarkerSize', markerSize);
     title(plotTitle);
     xlabel('Matrix element');
     legend('Sorted Reference', 'Reordered Prediction');
     hold off;
end