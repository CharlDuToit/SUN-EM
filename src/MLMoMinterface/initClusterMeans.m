function [clusterMeans, numClass, numNoClass] = initClusterMeans(prop,clusterSizes, threshDist, minClusterSize, minPointsDynamicRegion)
  % Assume planar, similar edge lengths
  % prop 1: distance
  % prop 2: edge direction dot edge direction
  % prop 3: edge direction dot displacement
  
   [~, numProp] = size(prop);
  
   %Initialise cluster size of each property
    numDirDotDirClusters = clusterSizes(1); 
    numDirDotDispClusters = clusterSizes(2);  
    numPreThreshDistClusters = clusterSizes(3);
    numPostThreshDistClusters = clusterSizes(4);
    
    maxDist = max(prop(:,1));
    minDist = min(prop(:,1));
    
    % intervals for properties
    %distInterval = hypspace(minDist, threshDist, numPreThreshDistClusters);
    distInterval = linspace(minDist, threshDist,numPreThreshDistClusters );
    numDistClusters = numPostThreshDistClusters + numPreThreshDistClusters -1;
    distInterval(numPreThreshDistClusters:numDistClusters) = linspace(threshDist, maxDist, numPostThreshDistClusters);
    dirDotDirInterval = linspace( -1 , 1 , numDirDotDirClusters);
    dirDotDispInterval = linspace( -1 , 1 , numDirDotDispClusters);  

    %initialise
    
    numNoClass = 0;
    numClass = 0;
    
    distOffset = 0;
    dirDotDirOffset = 0;
    dirDotDispOffset = 0;
    
    maxIncrease = 0.5;
    numSteps = 10;
    totalClusters = 0;
    maxClusters = (numDistClusters- 1) * (numDirDotDirClusters -1) * (numDirDotDispClusters-1);
    clusterMeans = zeros(numel(maxClusters), numProp);
  
    %Absolote value of dot products (might want negative values for non
    %planar structures)
   % prop(:,2) = abs(prop(:,2));
    
    for i = 1:(numDistClusters-1)
        low_dist = distInterval(i) + distOffset;
        high_dist = distInterval(i+1) +distOffset;
        for j = 1:(numDirDotDirClusters-1)
            if (j == 1)
                dirDotDirOffset = 0;
            end
            low_dir_dot_dir = dirDotDirInterval(j)+ dirDotDirOffset;
            high_dir_dot_dir = dirDotDirInterval(j+1) + dirDotDirOffset;
            for k = 1:(numDirDotDispClusters-1)
                if (k == 1)
                    dirDotDispOffset = 0;
                end
                low_dir_dot_disp = dirDotDispInterval(k) + dirDotDispOffset;
                high_dir_dot_disp = dirDotDispInterval(k+1) + dirDotDispOffset;
                     
                if (low_dist > maxDist || low_dir_dot_dir > 1 || low_dir_dot_disp > 1) 
                    break 
                end
                
                ind = find(prop(:,1) >= low_dist & prop(:,1) <= high_dist & ...
                    prop(:,2) >= low_dir_dot_dir & prop(:,2) <= high_dir_dot_dir &...
                    prop(:,3) >= low_dir_dot_disp & prop(:,3) <= high_dir_dot_disp  );              
                [N, ~]= size(ind);
                
                if (N >= minClusterSize)                   
                    
                    numClass = numClass + N;
                    totalClusters = totalClusters + 1;
                    clusterMeans(totalClusters, :) = sum(prop(ind,:),1)/N;

                % dynmaic region size changes    
                else 
                    if (N >= minPointsDynamicRegion) % dont increase a region if there are few points in it
                        %distStep = (maxIncrease / numSteps) * (distInterval(i+1) - distInterval(i));
                        %dirDotDirStep = (maxIncrease / numSteps) * (dirDotDirInterval(j+1) - dirDotDirInterval(j));
                        dirDotDispStep = (maxIncrease / numSteps) * (dirDotDispInterval(k+1) - dirDotDispInterval(k));
                        
                        step = 0;
                        while (N < minClusterSize && step < numSteps)
                            step = step + 1;
                            %distOffset = distOffset + distStep;
                            %dirDotDirOffset = dirDotDirOffset + dirDotDirStep;
                            dirDotDispOffset = dirDotDispOffset + dirDotDispStep;
                            
                            %high_dist = high_dist + distStep;
                            %high_dir_dot_dir = high_dir_dot_dir + dirDotDirStep;                       
                            high_dir_dot_disp = high_dir_dot_disp + dirDotDispStep;
                            
                            ind = find(prop(:,1) >= low_dist & prop(:,1) <= high_dist & ...
                                prop(:,2) >= low_dir_dot_dir & prop(:,2) <= high_dir_dot_dir &...
                                prop(:,3) >= low_dir_dot_disp & prop(:,3) <= high_dir_dot_disp  );
                            
                            [N, ~]= size(ind);
                            
                        end % while 
                    end % N > 0
                    if (N >= minClusterSize)   
                        numClass = numClass + N;
                        totalClusters = totalClusters + 1;
                        clusterMeans(totalClusters, :) = sum(prop(ind,:),1)/N;      
                    else
                        numNoClass = numNoClass + N;

                    end % if (N >= numTerms + 1)  
                end % if (N >= numTerms + 1)                
            end % for k
        end % for j
    end % for i
    clusterMeans = clusterMeans(1:totalClusters , :);
    
end

function [interval] = hypspace(low, high, numPoints)
    interval = zeros(1,numPoints);
    if (numPoints == 2)
        interval(1) = low;
        interval(2) = high;
    end
    delta = (1/low - 1/high)/(numPoints-1);
    for k = 1:numPoints
        interval(k) = 1/(1/low - (k -1)*delta );
    end


end