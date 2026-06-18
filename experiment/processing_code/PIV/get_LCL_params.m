function p = get_LCL_params(exp_name)
% get_LCL_params  Backward-compatible wrapper around get_piv_params.
%   Kept so existing longitudinal scripts keep working.
%   New code should call get_piv_params directly (supports LCL/LCTA/LCTB).
%
%   p = get_LCL_params('ExpLCL_2_01')
p = get_piv_params(exp_name);
end
