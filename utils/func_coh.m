function [COH, COH2DMAT, GCCM] = func_coh(Y_1, Y_2)
    COH = ( (Y_1 + Y_2) .* conj(Y_1 + Y_2) ) ./ ( 2 * ( Y_1 .* conj(Y_1) + Y_2 .* conj(Y_2) ) );
    GCCM = mean(COH,'all');
    % Compute the 2D average along the 3rd dimension.
    COH2DMAT = mean(COH, 3);
end
