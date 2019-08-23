function output = stepOp(L, T, h, tauDNA, lL)
% Operator equivalent to applying nucleosome translation, nucleosome
% rotation, and linker translation to a triad. This is an relative
% operator, and hence should be multiplied to the right of a triad matrix.

omega_trs = nucTrsOp(L, T, h);
omega_rot = nucRotOp(L, T, h, tauDNA, lL);      % Include twist from linker also in nucRot (this is cheating)
omega_link = linkTrsOp(lL);

output = omega_trs * omega_rot * omega_link;
end

