%% Parameters
L0 = 10;           % inter-bead distance
points = 120;      % number of beads
Tw = 24*pi;     % twist ‰(rad)

L = L0 * points;       % circumference
r = L/(2*pi);         % radius

tau = Tw/(2*pi*r);      % torsion to be introduced

%% Construct circle on xy plane
% Note that the last extra point is the same as the first point

s = 0:(2*pi*r)/points:2*pi*r;       % arclength parameter
omegaTwC = zeros(4, 4*(points+1));            % config, circle with twist
for ii = 1:points+1
    thisS = s(ii);
    
    rC = [r*cos(thisS/r); r*sin(thisS/r); 0];
    oC = [
        [-cos(thisS/r) 0 -sin(thisS/r)];
        [-sin(thisS/r) 0 cos(thisS/r)];
        [0 1 0]];
    omegaC = [                            % config, circle without twist
        [oC rC];
        [zeros(1, 3) 1]];
        
    TrotZ = makehgtform('zrotate', tau*thisS);
    
    omegaTwC(:, 4*ii-3:4*ii) = omegaC * TrotZ;
end

%% Plot

refCircle = circleXY(r, 100);

figure
hold on
for ii = 1:points
    nTwC = 10*omegaTwC(1:3, ii*4-3);
    bTwC = 10*omegaTwC(1:3, ii*4-2);
    uTwC = 10*omegaTwC(1:3, ii*4-1);
    rTwC = omegaTwC(1:3, ii*4);
    scatter3(rTwC(1), rTwC(2), rTwC(3))
    quiver3(rTwC(1), rTwC(2), rTwC(3), nTwC(1), nTwC(2), nTwC(3));
    quiver3(rTwC(1), rTwC(2), rTwC(3), bTwC(1), bTwC(2), bTwC(3));
    quiver3(rTwC(1), rTwC(2), rTwC(3), uTwC(1), uTwC(2), uTwC(3));
end
plot3(refCircle(1,:), refCircle(2,:), refCircle(3,:), 'DisplayName', 'reference circle')
xlabel('x')
ylabel('y')
zlabel('z')
title('Twisted circle')

%% Export

[n_out, b_out, u_out, r_out] = extractConfig(omegaTwC);
dlmwrite('r0', r_out, 'delimiter', '\t', 'precision', '%.12f')
dlmwrite('u0', u_out, 'delimiter', '\t', 'precision', '%.12f')
dlmwrite('v0', b_out, 'delimiter', '\t', 'precision', '%.12f')

