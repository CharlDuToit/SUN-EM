function [clusterMeans, numClass, numNoClass] = initClusterMeansDynamic(properties,clusterIntervals,  minSize)
   
    [~, numProp] = size(properties);
    
    numNoClass = 0;
    numClass= 0;
    intervalLengths = zeros(numProp, 1);
    
    maxNumClusters = 1;
    for p = 1:numProp
        intervalLengths(p) = numel(clusterIntervals{p}(:));
        maxNumClusters = maxNumClusters .* (intervalLengths(p)-1);
    end
    clusterMeans = zeros(maxNumClusters, numProp);
    
    %=========================
    %lowHigh = zeros(numProp, 2);
    intervalIndices = ones(numProp, 1);
    lowestPropChanged = 1;
    totalClusters = 0;
    tempProp = cell(numProp-1, 1);
    for c = 1:maxNumClusters
        
        %N = 0;
        for p = lowestPropChanged:numProp
            %prop = properties(;
            %lowHigh(p, : ) = [ clusterIntervals{p}(intervalIndices(p)), clusterIntervals(p,intervalIndices(p)+1) ];
            low = clusterIntervals{p}(intervalIndices(p));
            high = clusterIntervals{p}(intervalIndices(p)+1 );
            if (p == 1) 
                ind = find(properties(:,p) >= low & properties(:,p) <= high );
                prop = properties(ind, :);
                if (numProp > 1)
                    tempProp{1} = prop;
                end
            elseif (p >= lowestPropChanged) % p >= 2
                ind = find(tempProp{p-1}(:,p) >= low & tempProp{p-1}(:,p)  <= high );
                prop = tempProp{p-1}(ind, :);
                if (p < numProp)
                    tempProp{p} = prop;
                end
            end
%             [N, ~]= size(ind);
%             if (N < minSize)
%                 break
%             end
            
        end %  for p = 1:numProp
        
        %--------------------
        % Create cluster
        [N, ~]= size(ind);
        if (N < minSize)
            numNoClass = numNoClass+ N;
        else
            numClass = numClass + N;
            totalClusters = totalClusters + 1;
            clusterMeans(totalClusters, :) = sum(prop,1)/N;
        end
        
        %--------------------
        % adjust intervalIndices
        lowestPropChanged = numProp;
        for p = numProp:-1:1
            if ( p == numProp )
                if(intervalIndices(p) == intervalLengths(p) -1 )
                    intervalIndices(p) = 1;
                    if( p > 1 )
                        % more than 1 property
                        intervalIndices(p-1) = intervalIndices(p-1) + 1;
                        lowestPropChanged = p -1;
                    else
                        %should never execute. test this
                        break;
                    end
                else
                    intervalIndices(p) = intervalIndices(p) + 1;
                end
            elseif (intervalIndices(p) == intervalLengths(p) )
                intervalIndices(p) = 1;
                if( p > 1 )
                    intervalIndices(p-1) = intervalIndices(p-1) +1;
                    lowestPropChanged = p -1;
                else
                    % LAST LOOP C= MAX+1
                    % should never execute. test this
                    break;
                end
                
            else
                % lower properties wont change
                break;
            end
            
        end %for p = numProp:1
   
    end%for c = 1:maxNumClusters
    clusterMeans = clusterMeans(1:totalClusters, :);

end
