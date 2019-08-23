function output = helixScrew(r0, oHelix, L, T, h, points)
% Helix to fit points from rigid body operations
% r0: position of the first point
% oHelix: orientation of the helix with t3 as the axis of rotation, and t1
% pointing toward the first point (in lab frame)

% Parameters
b = h/(2*pi);
a = sqrt((L/T)^2 - b^2);
m = sqrt(a^2 + b^2);
s = 0:(T*m/points):T*m;
sSize = size(s);
sSize = sSize(2);

% Canonical helix along z-axis, with 1 padded as the 4th coordinate
H_canon = [helixBase(s, a, b); ones(1, sSize)];

% Transport the helix onto the frame defined by r0 and oHelix
% rHelix is the vector from the first point in the canonical helix 
% (a, 0, 0) to the first point in the target helix.
% oHelix, which transforms the orientation of the canonical helix to the
% one in the target helix, is the orientation of the target helix since the
% canonical one is simply the identity matrix.
rHelix = r0 - oHelix * [a; 0; 0];
H_canon_helix = [
    [oHelix rHelix];
    [zeros(1,3) 1]];

% Subsequently apply the frame transformation to obtain the helix
% in the lab frame
output = [];
for ii = 1:sSize
    thisH_canon = H_canon(:, ii);
    thisH_helix = H_canon_helix * thisH_canon;
    output = [output thisH_helix];
end
end

