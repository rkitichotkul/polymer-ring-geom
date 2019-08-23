function [n, b, u, r] = extractConfig(c)
% Extract orientation and position from configuration matrix in homogeneous
% coordinate (4D)
% Last coordinate is discarded

sizeC = size(c);
sizeC = sizeC(2);

points = sizeC / 4;
n = zeros(3, points);
b = zeros(3, points);
u = zeros(3, points);
r = zeros(3, points);

for ii = 1:points
    n(1:3, ii) = c(1:3, ii*4-3);
    b(1:3, ii) = c(1:3, ii*4-2);
    u(1:3, ii) = c(1:3, ii*4-1);
    r(1:3, ii) = c(1:3, ii*4);
end

% Discard last point
n = n(1:3, 1:points-1);
b = b(1:3, 1:points-1);
u = u(1:3, 1:points-1);
r = r(1:3, 1:points-1);

n = transpose(n);
b = transpose(b);
u = transpose(u);
r = transpose(r);

end

