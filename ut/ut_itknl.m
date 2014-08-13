function kernels=ut_itknl(siz)
% kernels=ut_itknl([5, 5])
% Creates a 3D kernel for block processing of metrics

m=siz(1);
n=siz(2);

kernels=zeros(m, n, m*n);

for i=1:m*n
    k=kernels(:, :, i);
    k(i)=1;
    kernels(:, :, i)=k;
end
