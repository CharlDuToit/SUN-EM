function  plotMLMOM(mlmom, activePlots, gridSize, freqInd)
    %freqInd : scalar index
    %  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21]
    % [0 0 0 0 0 0 0 1 0 0  0  0  0  0  0  0  0  0  0  0  0]
    % [0 0 0 0 0 0 1 0 0 0  0  0  0  0  0  0  0  0  0  0  0]
    
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
    prop = mlmom.nonSingZmnProp;
    numPlots = numel(activePlots);
    for k = 1:numPlots
        if (activePlots(k) == 0)
            continue
        end
        switch k
            %real surface
            case 1 % 1 Surface Error Real Unity 3 Prop
                err = model.unityNonSingZmn - model.refNonSingZmn;
                plotTitle = ('Unity Real Difference Error ');
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 2% 2 Surface Error Real Pred 3 Prop
                model = mlmom.weightModels{freqInd, 1};
                err = model.predNonSingZmn - model.refNonSingZmn;
                plotTitle = ('Prediction Real Difference Error');
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 3% 3 Surface Error Real Squared Unity 3 Prop
                model = mlmom.weightModels{freqInd, 1};
                err = (model.unityNonSingZmn - model.refNonSingZmn).^2;
                plotTitle = ('Unity Real Squared Difference Error');
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 4 % 4 Surface Error Real Squared Pred 3 Prop
                model = mlmom.weightModels{freqInd, 1};
                err = (model.predNonSingZmn - model.refNonSingZmn).^2;
                plotTitle = ('Prediction Real Squared Difference Error');
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);       
            case 5% 5 Surface Error Real Log Norm Squared Unity 3 Prop
                model = mlmom.weightModels{freqInd, 1};
                err = log((model.unityNonSingZmn - model.refNonSingZmn).^2./ model.refNonSingZmn.^2);
                plotTitle = ('Unity Real Log of Normalised Squared Difference Error ');
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 6% 6 Surface Error Real Log Norm Squared Pred 3 Prop
                model = mlmom.weightModels{freqInd, 1};
                err = log((model.predNonSingZmn - model.refNonSingZmn).^2./ model.refNonSingZmn.^2);
                plotTitle = ('Unity Real Log of Normalised Squared Difference Error ');
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            % Imaginary surface
            case 7% 7 Surface Error Imag Unity 3 Prop
                model = mlmom.weightModels{freqInd, 2};
                err = model.unityNonSingZmn - model.refNonSingZmn;
                plotTitle = ('Unity Imag Difference Error ');
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 8% 8 Surface Error Imag Pred 3 Prop
                model = mlmom.weightModels{freqInd, 2};
                err = model.predNonSingZmn - model.refNonSingZmn;
                plotTitle = ('Prediction Imag Difference Error');
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 9% 9 Surface Error Imag Squared Unity 3 Prop
                model = mlmom.weightModels{freqInd, 2};
                err = (model.unityNonSingZmn - model.refNonSingZmn).^2;
                plotTitle = ('Unity Imag Squared Difference Error');
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 10% 10 Surface Error Imag Squared Pred 3 Prop
                model = mlmom.weightModels{freqInd, 2};
                err = (model.predNonSingZmn - model.refNonSingZmn).^2;
                plotTitle = ('Prediction Imag Squared Difference Error');
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 11 % 11 Surface Error Imag Log Norm Squared Unity 3 Prop
                model = mlmom.weightModels{freqInd, 2};
                err = log((model.unityNonSingZmn - model.refNonSingZmn).^2 ./ model.refNonSingZmn.^2);
                plotTitle = ('Unity Imag Log of Normalised Squared Difference Error ');
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 12% 12 Surface Error Imag Log Norm Squared Pred 3 Prop
                model = mlmom.weightModels{freqInd, 2};
                err = log((model.predNonSingZmn - model.refNonSingZmn).^2 ./ model.refNonSingZmn.^2);
                plotTitle = ('Prediction Imag Log of Normalised Squared Difference Error ');
                plotSurfaceError3Prop(prop, err, gridSize, plotTitle);
            case 13
            case 14
            case 15
            case 16
            otherwise
                
        end % switch
    end

end

function  plotSurfaceError3Prop(prop, err, gridSize, plotTitle)

    distLin = linspace(min(prop(:,1)) , max(prop(:,1)), gridSize);
    dirDotDirLin = linspace(min(prop(:,2)) , max(prop(:,2)), gridSize);
    dirDotDispLin = linspace(min(prop(:,3)) , max(prop(:,3)), gridSize);
    
    subplot(2,2,1);
    % x = dist, y = dir dot dir   
    [X , Y] = meshgrid(distLin, dirDotDirLin);
    Z = griddata(prop(:,1) , prop(:, 2) , err, X , Y, 'cubic');  
    mesh(X,Y,Z);   
    axis tight; hold on;
    plot3(prop(:,1) , prop(:, 2) , err , '.', 'MarkerSize', 5);
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
    plot3(prop(:,1) , prop(:, 3) , err , '.', 'MarkerSize', 5);
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
    plot3(prop(:,2) , prop(:, 3) , err , '.', 'MarkerSize', 5);
    title(plotTitle);
    xlabel('Edge direction dot edge direction');
    ylabel('Edge direction dot displacement');
end