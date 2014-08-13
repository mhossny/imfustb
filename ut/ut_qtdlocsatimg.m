function infty=ut_qtdlocsatimg(a, md)
% This function calculates saturated subimages based on the quadtree map.
% If you are using it please refer to our work in [1-3].
% 
% [1] Hossny, M.; Nahavandi, S.; Creighton, D.; Bhatti, A, "Image fusion performance
% metric based on mutual information and entropy driven quadtree decomposition," 
% Electronics Letters , vol.46, no.18, pp.1266,1268, September 2010
% doi: 10.1049/el.2010.1778
% 
% [2] Hossny, M.; Nahavandi, S.; Creighton, D., "Comments on 'Information measure 
% for performance of image fusion'," Electronics Letters , vol.44, no.18, 
% pp.1066,1067, August 28 2008 (doi: 10.1049/el:20081754)
% 
% [3] Hossny, M.; Nahavandi, S., "Measuring the capacity of image fusion," Image 
% Processing Theory, Tools and Applications (IPTA), IEEE International
% Conference on , pp.415-420, 2012. (doi:10.1109/IPTA.2012.6469548)

%ent_a=ut_entropy(a)

%md=2^(log2(md)+1)

ent_a=mean2(bk_entropy(a, zeros(md, md)))
ent_b=ent_a;

%[q_a, b_a, d_a]=ut_qtd(a, ent_a, 10);


entropythreshold=min(ent_a, ent_b);%*2/3;
entropythreshold=(ent_a+ ent_b)/2
entropythreshold=ent_a

%[q_a, b_a, d_a]=ut_qtd(a, entropythreshold, 30);

mindim=16;%2;
maxdim=size(a, 1)/2;%min(size(a, 1))/2;

im=a;

q=qtdecomp(im, @ut_qtentropy, entropythreshold, [mindim, maxdim]);

numblks=length(find(q));
dims=2.^[log2(maxdim):-1:log2(mindim)];
%dims=[512, 256, 128, 64, 32, 16, 8, 4, 2, 1];
numdims=length(dims);
blkid=1;
vrtid=2;
dd=512;
counter=0;

totalweight=0;
totalerr=0;

infty=uint8(zeros(size(a)));

for dim=dims
    clear blks;
    [ablks, r, c]=qtgetblk(a, q, dim);
    siz=size(ablks);
    numsubblks=size(ablks, 3);
    k=0;

    if ~numsubblks; continue; end;
    s=size(ablks, 1);
    step=1;
    eblks=[];
    for i=1:step:numsubblks
        tmpu=0:s*s-1;
        eblks(:, :, i)=reshape(tmpu(randperm(numel(tmpu))), [s,s]);
        eblks(:, :, i)=histeq(uint8(eblks(:, :, i)), 256);
    end
    
    eblks=reshape(eblks, siz);
    infty=uint8(qtsetblk(uint8(infty), q, dim+k, eblks));
end
