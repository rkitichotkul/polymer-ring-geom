% Transport position and orientation of points obtained from successive
% applications of step operator onto a circle whose circumference matches
% the length along z-axis of the helix lined up along the z-axis (the
% canonical one)

% Parameters for step operator
% https://github.com/brunobeltran/nuc_chain/blob/master/nuc_chain/__init__.py
% Paramters from Richmond and Davey (2003)
% Note that 1 bp = 0.34 nm
L = 147 * 0.34;     % wrapping (nm)
T = 1.67 * 2 * pi;  % turns (rad)
h = 2.59;           % pitch (nm) (25.9 A)
lL = 35 * 0.34;     % linker length (nm)
tauDNA = 2*pi/(10.43*0.34);         % DNA intrinsic twist (rad/nm) (10.43 bp/twist)
steps = 200;         % number of successive step operations applied (number of beads)

% Helical parameters for reference
b = h / (2*pi);             % pitch = 2*pi*b (bp)
a = sqrt((L/T)^2 - b^2);    % radius (bp)
m = sqrt(a^2 + b^2);

% Initial configuration
c0 = [
    [1 0 0 0];
    [0 1 0 0];
    [0 0 1 0];
    [0 0 0 1]];

% Twist of nucleosome chain (rad)
Tw = 0;

% Generate points
step = stepOp(L, T, h, tauDNA, lL);     % step operator
c = zeros(4, (steps + 1) * 4);  % config of points, 4 columns/point
c(:, 1:4) = c0;                 % initial config

for ii = 1:steps
    c(:, ii*4+1:ii*4+4) = c(:, ii*4-3:ii*4) * step;
end

% Transport points onto z-axis of fit-helix
[oScrew, LScrew, TScrew, hScrew] = screwLeftParam(step, steps);     % fit-helix parameters
bScrew = hScrew / (2*pi);
aScrew = sqrt((LScrew/TScrew)^2 - bScrew^2);

o_lab_fit = c0(1:3, 1:3) * oScrew;
r_lab_fit = c0(1:3, 4) - o_lab_fit * [aScrew; 0; 0];
T_lab_fit = [
    [o_lab_fit r_lab_fit];
    [zeros(1, 3) 1]];
        
c_canon = T_lab_fit \ c;        % Now this helix is represented in the canonical way (along z-axis)

% Transport points from fit-helix frame (in which helix is "canonical")
% into end-to-end vector frame
d_canon_ee = [aScrew; 0; 0];                % Translation
u_canon_ee = c_canon(:, end) - c_canon(:, 4);   % z-axis in end-to-end frame
u_canon_ee = u_canon_ee(1:3);
u_canon_ee = u_canon_ee / norm(u_canon_ee);
u_canon = [0; 0; 1];
theta = acos(dot(u_canon_ee, u_canon));     % note that norms are 1
k_canon_ee = cross(u_canon, u_canon_ee);    % axis of rotation
k_canon_ee = k_canon_ee / norm(k_canon_ee);
K_canon_ee = [
    [0 -k_canon_ee(3) k_canon_ee(2)];
    [k_canon_ee(3) 0 -k_canon_ee(1)];
    [-k_canon_ee(2) k_canon_ee(1) 0]];
R_canon_ee = eye(3) + sin(theta) * K_canon_ee + (1 - cos(theta)) * K_canon_ee * K_canon_ee;
T_canon_ee = [
    [R_canon_ee d_canon_ee];
    [zeros(1, 3) 1]];
c_ee = T_canon_ee \ c_canon;

% Transport position and orientation onto circle
rCir = bScrew * TScrew;         % Radius of circle matching length along helix's z axis
sCir = 0:(2*pi*rCir/steps):2*pi*rCir;
tau = Tw / (2*pi*rCir);

c_ee_z0 = c_ee;
c_ee_z0(3, 4:4:end) = 0;     % Set z-coordinate of position to 0

c_circle = zeros(4, (steps + 1) * 4);

for ii = 1:steps+1
    s = sCir(ii);
    omegaP = [
        [-cos(s/rCir).*cos(tau*s) -cos(s/rCir).*sin(tau*s) -sin(s/rCir)];
        [-sin(s/rCir).*cos(tau*s) -sin(s/rCir).*sin(tau*s) cos(s/rCir)];
        [-sin(tau*s) cos(tau*s) 0]];
    rP = [rCir*cos(s/rCir); rCir*sin(s/rCir); 0];
    TP = [                          % Parallel tranport onto circle matrix
        [omegaP rP];
        [zeros(1, 3) 1]];
    c_circle(:, ii*4-3:ii*4) = TP * c_ee_z0(:, ii*4-3:ii*4);
end

    
% Reference circle
refCircle = circleXY(rCir, 200);

%% Plot

r_circle = zeros(3, steps+1);
for ii = 1:steps+1
    r_circle(:, ii) = c_circle(1:3, ii*4);
end

figure
hold on
plot3(refCircle(1,:), refCircle(2,:), refCircle(3,:), 'DisplayName', 'Reference circle')
plot3(r_circle(1, :), r_circle(2, :), r_circle(3, :), 'DisplayName', 'Transported helix')
for ii = 1:steps+1
    n_circleThis = 10*c_circle(1:3, ii*4-3);
    b_circleThis = 10*c_circle(1:3, ii*4-2);
    u_circleThis = 10*c_circle(1:3, ii*4-1);
    r_circleThis = c_circle(1:3, ii*4);
    scatter3(r_circleThis(1), r_circleThis(2), r_circleThis(3))
    quiver3(r_circleThis(1), r_circleThis(2), r_circleThis(3), n_circleThis(1), n_circleThis(2), n_circleThis(3));
    quiver3(r_circleThis(1), r_circleThis(2), r_circleThis(3), b_circleThis(1), b_circleThis(2), b_circleThis(3));
    quiver3(r_circleThis(1), r_circleThis(2), r_circleThis(3), u_circleThis(1), u_circleThis(2), u_circleThis(3));
end
xlabel('x')
ylabel('y')
zlabel('z')
title('Points from successive applications of step operator, transported onto matching-size (z-axis and circumference) circle')

%% For export
[n_circle, b_circle, u_circle, r_circle] = extractConfig(c_circle);
dlmwrite('r0', r_circle, 'delimiter', '\t', 'precision', '%.12f')
dlmwrite('u0', u_circle, 'delimiter', '\t', 'precision', '%.12f')
dlmwrite('v0', n_circle, 'delimiter', '\t', 'precision', '%.12f')

