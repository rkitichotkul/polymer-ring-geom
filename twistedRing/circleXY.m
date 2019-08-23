function output = circleXY(r, points)
% Return a circle of radius r with appropriate arclength parameter s
s = 0:(2*pi*r)/points:(2*pi*r);
output = circleXYBase(s, r);
end

