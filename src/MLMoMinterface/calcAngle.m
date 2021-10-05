function angleRad = calcAngle(a,b)
    %angle between -pi and pi
    angleRad = atan2(norm(cross(a,b)), dot(a,b));
    if (angleRad > pi)
        angleRad = angleRad - 2*pi;
    end
end