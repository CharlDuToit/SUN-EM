function [err, ind] = calcMinClusterError(clusterMeans, prop, propScale, errorCode)
    %propScale: the scalars by which every property contribution to the
    %error should be multiplied
    %errorcode: the type of error function

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
    
    switch errorCode
        case 1 % 3 properties, dist, orientation, radial
            [err, ind] =min(propScale(1) * abs( log(clusterMeans(:,1)/prop(1) ) ) +...
                propScale(2)*abs(cos(clusterMeans(:,2)) - cos(prop(2))) + ...
                propScale(3)*abs(clusterMeans(:,3) - prop(3)));
        case 2 % 2 new properties, 2 old properties, no orientation
            [err, ind] =min( propScale(1) * abs( log(clusterMeans(:,1)/prop(1) ) ) +...
                propScale(2)*abs(clusterMeans(:,2) - prop(2))  +...
                propScale(3) * abs( log(clusterMeans(:,3)/prop(3) ) )+... 
                 propScale(4)*abs(clusterMeans(:,4) - prop(4))    );
        otherwise
            err = 0;
            ind = 0;
    end


end