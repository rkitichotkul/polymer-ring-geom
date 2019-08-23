function output = linkTrsOp(lL)
% Linker translation operator
% This is a relative operator and should be multiplied to the right of a
% triad matrix
d = [0; 0; lL];
output = [
    [eye(3) d];
    [zeros(1, 3) 1]];
end

