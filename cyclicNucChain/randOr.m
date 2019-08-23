function output = randOr()
% Random orientation (3 vectors concatenated along row)
u = randUnitVec();
n = randUnitVec();
b = cross(u, n);
b = b / norm(b);
n = cross(b, u);
output = [n b u];
end

