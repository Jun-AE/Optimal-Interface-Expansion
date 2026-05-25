function [LAC,LAC2DMAT,OverallLAC] = func_lac(Y_1,Y_2)

numerator = 2*(abs((conj(Y_1).*Y_2)));
denumerator = (Y_1 .* conj(Y_1)) + (Y_2 .* conj(Y_2));
LAC = numerator ./ denumerator;

OverallLAC=mean(LAC,'all');
LAC2DMAT=mean(LAC,3);


