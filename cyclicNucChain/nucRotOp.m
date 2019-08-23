function output = nucRotOp(L, T, h, tauDNA, lL)
% Rotation matrix of nucleosome operator

% Parameters
b = h/(2*pi);
a = sqrt((L/T)^2 - b^2);
m = sqrt(a^2 + b^2);

% Inverse of standard entry frame
omegaEntryinv = [
    [-1 0 0];
    [0 -b/m -a/m];
    [0 -a/m b/m]];

% Exit frame
omegaExit = [
    [-cos(L/m) -(b/m)*sin(L/m) -(a/m)*sin(L/m)];
    [sin(L/m) -(b/m)*cos(L/m) -(a/m)*cos(L/m)];
    [0 -(a/m) (b/m)]];

% Rotation along z to fix twist
theta_helix = (b*L)/(m^2);         % Detwist inherent torsion of helix 
theta_DNA = lL * tauDNA;            % Twist to the amount natural in DNA
theta = theta_helix + theta_DNA;
rotZ = [
    [cos(theta) -sin(theta) 0];
    [sin(theta) cos(theta) 0];
    [0 0 1]];

% Rotation operator
output = omegaEntryinv * omegaExit * rotZ;

% Put into homogeneous coordinate
output = [
    [output zeros(3, 1)];
    [zeros(1, 3) 1]];

end

