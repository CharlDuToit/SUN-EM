function [err, ind] = calcMinClusterError(clusterMeans, prop, propScale)
    %propScale: the scalars by which every property contribution to the
    %error should be multiplied

    %err: error of mean which has the smallest error
    %ind: index of mean which has the smallest error
    
%    [err, ind] =min(propScale(1) * abs( log(clusterMeans(:,1)/prop(1) ) ) +...
%    propScale(2)*abs(clusterMeans(:,2) - prop(2)) + ...
%    propScale(3)*abs(clusterMeans(:,3) - prop(3)));

%    if (  prop(2) ~= pi/2 )
%    %if ((clusterMeans(:,2) ~= 0) && (prop(2) ~= 0 ))
%        [err, ind] =min(propScale(1) * abs( log(clusterMeans(:,1) / prop(1) ) ) +...
%            propScale(2)*abs( log( abs(cos(clusterMeans(:,2)) / cos(prop(2))) ) ) + ...
%            propScale(3)*abs(clusterMeans(:,3) - prop(3)));
%    else
%        [err, ind] =min(propScale(1) * abs( log(clusterMeans(:,1)/prop(1) ) ) +...
%            propScale(2)*abs(clusterMeans(:,2) - prop(2)) + ...
%            propScale(3)*abs(clusterMeans(:,3) - prop(3)));
%    end

   [err, ind] =min(propScale(1) * abs( log(clusterMeans(:,1)/prop(1) ) ) +...
   propScale(2)*abs(cos(clusterMeans(:,2)) - cos(prop(2))) + ...
   propScale(3)*abs(clusterMeans(:,3) - prop(3)));

end