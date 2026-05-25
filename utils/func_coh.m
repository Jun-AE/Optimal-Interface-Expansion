% function [COH, COH2DMAT, OverallCOH] = func_coh(Y_1, Y_2)
% % Preserve original calculation structure with MATLAB-optimized operations
% sum_Y = Y_1 + Y_2;
% COH = (sum_Y .* conj(sum_Y)) ./ (2 * (Y_1.*conj(Y_1) + Y_2.*conj(Y_2)));
% 
% OverallCOH = mean(COH, 'all');
% COH2DMAT = mean(COH, 3);
% end
% 

function [COH, COH2DMAT, OverallCOH] = func_coh(Y_1, Y_2)
   
    COH = ( (Y_1 + Y_2) .* conj(Y_1 + Y_2) ) ./ ( 2 * ( Y_1 .* conj(Y_1) + Y_2 .* conj(Y_2) ) );
    
    OverallCOH = mean(COH(:));
  
    % Compute the 2D average along the 3rd dimension.
    COH2DMAT = mean(COH, 3);
end