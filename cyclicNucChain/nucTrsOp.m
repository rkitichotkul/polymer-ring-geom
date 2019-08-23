function output = nucTrsOp(L, T, h)
% Operator of nucleosome translation
% This is a relative operator and should be multiplied to the right of a
% triad matrix

% Parameters
b = h / (2 * pi);
a = sqrt((L/T)^2 - b^2);
m = sqrt(a^2 + b^2);

% Translation operator in helix frame
d_helix_trs = [a*(cos(T) - 1); -a*sin(T); b*T];

omega_helix_trs = [
    [eye(3) d_helix_trs];
    [zeros(1,3) 1]];
% Transform helix frame to the frame of the entry triad
R_entry_helix = [
    [-1 0 0];
    [0 -b/m -a/m];
    [0 -a/m b/m]];

d_entry_helix = [a; 0; 0];

H_entry_helix = [
    [R_entry_helix d_entry_helix];
    [zeros(1,3) 1]];
% The nucleosome translation operator in the entry frame
output = H_entry_helix * omega_helix_trs * inv(H_entry_helix);


end

