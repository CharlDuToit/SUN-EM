function angleRad = calcAngle(a,b)
    %angle between -pi and pi
    %angleRad = atan2(norm(cross(a,b)), dot(a,b));
    angleRad = atan2(a(1)*b(2)-b(1)*a(2),a(1)*b(1)+a(2)*b(2));
    %angleRad2 = 2 * atan(norm(a*norm(b) - norm(a)*b) / norm(a * norm(b) + norm(a) * b));
    %angleRad3 = subspace(a,b);
    %angleRad4 = atan2(a(1)*b(2)-b(1)*a(2),a(1)*b(1)+a(2)*b(2));
    %angleRad4 = atan2(a(1)*b(2)-b(1)*a(2),a(1)*a(2)+b(1)*b(2));
    %angleRad5 = atan2(norm(cross(b,a)), dot(b,a));
    %if (angleRad <0 || angleRad2 < 0)
        %b= 5;
    %end
    %if (angleRad > pi)
        %angleRad = angleRad - 2*pi;
    %end
end