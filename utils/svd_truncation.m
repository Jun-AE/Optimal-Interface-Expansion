function [FRF_TSVD] = svd_truncation(FRF,reduction)

num_page=size(FRF,3);

[U,s,V]=pagesvd(FRF,'vector');

VH=pagectranspose(V);

kk= size(s,1)-reduction;

Uk=U(:,1:kk,:);
Vk= VH(1:kk,:,:);

Sk=zeros(kk,kk,num_page);

for i=1:num_page
    Sk(:,:,i) = diag(s(1:kk, :, i));
end

US = pagemtimes(Uk, Sk);
FRF_TSVD = pagemtimes(US, Vk);

end