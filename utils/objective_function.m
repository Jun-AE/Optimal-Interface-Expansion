function cor = objective_function(x, frequency_range, y_n, y_e, c, e, ...
                                  truncation, reduction, trust_func_flag, ...
                                  interface_dofs, method, dim, trust_func_val, gen_excit)
% OBJECTIVE_FUNCTION  GCCM fitness for SEMM expansion with sensor/excitation selection.
%
%   cor = objective_function(x, frequency_range, y_n, y_e, c, e, ...
%                            truncation, reduction, trust_func_flag, ...
%                            interface_dofs, method, dim, trust_func_val)
%
%   dim = 1 : x picks a row of c; sensors = excitations (drive-point).
%   dim = 2 : x = [i, j]; sensors = c(i,:), excitations = e(j,:).
%   dim = 3 : x = [i, j]; sensors = c(i,:), excitations = e(j,:,i)
%             (or excitationCombs(j,:) when gen_excit is supplied).
%
%   Reference:
%     Junaid et al. (2026), Journal of Sound and Vibration.

if nargin < 14, gen_excit = []; end
x = round(x);

if trust_func_flag
    tf.freq_Hz   = trust_func_val;
    tf.steepness = 15;
else
    tf = [];
end

switch dim
    case 1
        cdofs = c(x, :);
        edofs = cdofs;
    case 2
        cdofs = c(x(1), :);
        edofs = e(x(2), :);
    case 3
        cdofs = c(x(1), :);
        if isempty(gen_excit)
            edofs = e(x(2), :, x(1));
        else
            exc_combs = gen_excit(cdofs);
            edofs     = exc_combs(x(2), :);
        end
end

ys = semm(cdofs, edofs, y_n, y_e, frequency_range, ...
          truncation, reduction, trust_func_flag, tf);

ys_sub = ys(interface_dofs, interface_dofs, :);
ye_sub = y_e(interface_dofs, interface_dofs, :);

switch upper(method)
    case 'COH'
        [~, ~, cor] = func_coh(ys_sub, ye_sub);
    case 'LAC'
        [~, ~, cor] = func_lac(ys_sub, ye_sub);
end

end
