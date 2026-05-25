function [frequency_range_rad, frequency_range_Hz] = frequency_generation(freqstart, freqend, stepsize)
% FREQUENCY_GENERATION  Build a uniform frequency axis in Hz and rad/s.
%
%   [W, F] = FREQUENCY_GENERATION(FREQSTART, FREQEND, STEPSIZE) returns a
%   uniform frequency axis F = FREQSTART:STEPSIZE:FREQEND (in Hz) and its
%   angular counterpart W = 2*pi*F (in rad/s).
%
%   Inputs:
%     FREQSTART - (double, scalar) [Hz] Lower bound of the axis.
%     FREQEND   - (double, scalar) [Hz] Upper bound of the axis.
%     STEPSIZE  - (double, scalar) [Hz] Uniform step.
%
%   Outputs:
%     FREQUENCY_RANGE_RAD - (1×N double) [rad/s] Angular frequency axis.
%     FREQUENCY_RANGE_HZ  - (1×N double) [Hz]    Frequency axis. The two
%                           outputs satisfy W = 2*pi*F element-wise.
%
%   Example:
%     [w, f] = frequency_generation(1, 700, 1);
%
%   See also: compute_frf, semm.

frequency_range_Hz  = freqstart : stepsize : freqend;
frequency_range_rad = 2 * pi * frequency_range_Hz;

end
