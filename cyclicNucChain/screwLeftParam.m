function [oScrew, LScrew, TScrew, hScrew] = screwLeftParam(operator, steps)
% Determine screw parameters from a given transformation by using Chasles
% theorem (Screw-axis theorem). This function returns parameters of
% fit-left-handed helix.
% Note that the return values are different from those of the associated
% right-handed helix only for LScrew and TScrew. The orientations are the
% same since we draw left-handed helix by reversing a sign in the
% parametric equation.
% operator: transformation matrix in homogeneous coordinate
% steps: number of steps in which the operator is applied

% Rotation and displacement in rigid displacement operator
R = operator(1:3, 1:3);
d = operator(1:3, 4);

% Sweep angle
phi = acos(0.5 * (trace(R) - 1));
% Turns
% Here, 2*pi - phi is the modification of Chasles theorem to fit
% left-handed helix
TScrew = (2*pi - phi) * steps;

% Axis of rotation
% Since the initial orientation is identity matrix, the reference
% orientation here is the lab frame.
uhat = (1 / (2 * sin(phi))) * (R - transpose(R));
u = [uhat(3, 2); uhat(1, 3); uhat(2, 1)];
u = u / norm(u);

% Displacement along axis of rotation
k = dot(u, d);
%  Screw pitch
% Here, 2*pi - phi is the modification of Chasles theorem to fit
% left-handed helix
hScrew = (2 * pi * k) / (2*pi - phi);

% Find the vector from the initial point to a point on the axis of rotation
% Projection of dOp onto plane perpendicular to u
dp = d - k * u;
% Vectors orthogonal to u
v = dp / norm(dp);
w = cross(u, v);
% Vector to a point on the axis of rotation
c = 0.5 * norm(dp) * (v + (sin(phi) / (1 - cos(phi))) * w);
% Screw radius
rScrew = norm(c - dot(c, u) * u);

% Arc length
LScrew = TScrew * sqrt(rScrew^2 + (hScrew/(2*pi))^2);

% Change of orientation: from that in the canonical helix (identity matrix)
% to the one corresponding to the screw operator
t1Screw = -c/norm(c);
t3Screw = u;
t2Screw = cross(t3Screw, t1Screw);
oScrew = [t1Screw t2Screw t3Screw];

end

