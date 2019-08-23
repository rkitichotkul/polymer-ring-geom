function output = circleXYBase(s, r)
% Return (x, y, z) of a circle on xy-plane according to the arclength s
% and the radius r
sizeS = size(s);
sizeS = sizeS(2);
output = [
    r*cos(s/r);
    r*sin(s/r);
    zeros([1, sizeS])
    ];
end

