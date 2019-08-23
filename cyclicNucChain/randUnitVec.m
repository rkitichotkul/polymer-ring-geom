function output = randUnitVec()
% Generate a random unit vector based on uniform distribution
output = rand(3, 1);

% Normalize
output = output / norm(output);
end

