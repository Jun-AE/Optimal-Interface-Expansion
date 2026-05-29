function FRF_TSVD = svd_truncation(FRF, num_drop)
% svd_truncation  Page-wise TSVD: discard trailing singular values per frequency.
%
%   FRF_TSVD = svd_truncation(FRF, NUM_DROP)
%
%   For each frequency page of FRF, computes an SVD and reconstructs using only
%   the leading (size(s,1) - NUM_DROP) singular values. Equivalently, NUM_DROP
%   trailing singular values are discarded (largest rank reduction at the tail).
%
%   Inputs:
%     FRF      - (n x m x nFreq) complex FRF array
%     NUM_DROP - (non-negative integer) number of singular values to DROP from
%                each page. Example: NUM_DROP = 15
%                removes the 15 smallest singular values when the page is
%                square with at least 15 singular values.
%
%   Output:
%     FRF_TSVD - same size as FRF, rank-reduced reconstruction
%
%   See also: semm, pagesvd.

if nargin < 2 || isempty(num_drop)
    num_drop = 0;
end

num_page = size(FRF, 3);

[U, s, V] = pagesvd(FRF, 'vector');

VH = pagectranspose(V);

kk = size(s, 1) - num_drop;
if kk < 1
    error('svd_truncation:TooManyDropped', ...
        'num_drop=%d leaves no singular values (page has %d).', num_drop, size(s, 1));
end

Uk = U(:, 1:kk, :);
Vk = VH(1:kk, :, :);

Sk = zeros(kk, kk, num_page);

for i = 1:num_page
    Sk(:, :, i) = diag(s(1:kk, :, i));
end

US       = pagemtimes(Uk, Sk);
FRF_TSVD = pagemtimes(US, Vk);

end
