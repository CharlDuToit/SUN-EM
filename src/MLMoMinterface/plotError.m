function  plotError(prop, err, gridSize, plotTitle)
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