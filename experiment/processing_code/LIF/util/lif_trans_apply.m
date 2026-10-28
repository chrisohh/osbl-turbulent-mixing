function J = lif_trans_apply(T, I)
%LIF_TRANS_APPLY Resample a transverse frame onto the laser-sheet plane.
%
%   J = LIF_TRANS_APPLY(T, I)
%
% T from lif_trans_rectify.  Outside the camera's footprint J is 0 (the LIF
% frames are black where there is no dye, so 0 is the natural fill and keeps
% NaNs out of the surface / concentration code).
    J = interp2(I, T.U, T.V, 'linear', 0);
    J(isnan(J)) = 0;
end
