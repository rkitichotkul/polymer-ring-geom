function output = helixBase(s, a, b)
% Return (x, y, z) of a helix corresponding to the arclength s
% and the parameters a (radius) and b (pitch)
output = [
    a*cos(s/sqrt(a^2 + b^2)); 
    -a*sin(s/sqrt(a^2 + b^2)); 
    b*s/sqrt(a^2 + b^2);
    ];
end

