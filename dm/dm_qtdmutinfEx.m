function [err, errmap, err_a, errmap_a, err_b, errmap_b]=dm_qtdmutinfEx(a, b, md)
% Localised quadtree mutual information metric as modified by Hossny et al. in
% 
% Hossny, M.; Nahavandi, S.; Creighton, D.; Bhatti, A, "Image fusion performance
% metric based on mutual information and entropy driven quadtree decomposition," 
% Electronics Letters , vol.46, no.18, pp.1266,1268, September 2010
% doi: 10.1049/el.2010.1778
% 
% The MI is normalised by max entropy.
% 
% The MI is normalised by avergaing marginal entropies as explained by Hossny et al. in 
% 
% Hossny, M.; Nahavandi, S.; Creighton, D., "Comments on 'Information measure 
% for performance of image fusion'," Electronics Letters , vol.44, no.18, 
% pp.1066,1067, August 28 2008 (doi: 10.1049/el:20081754)


% init 
% ent_a=ut_entropy(a);
% ent_b=ut_entropy(b);
% 
% entropythreshold=min(ent_a, ent_b)*2/3;

ent_a=mean2(bk_entropy(a, zeros(md, md)));
ent_b=mean2(bk_entropy(b, zeros(md, md)));

entropythreshold=min(ent_a, ent_b);

%[qa, ba, da]=ut_qtd(a, entropythreshold, 0);
%[qb, bb, db]=ut_qtd(b, entropythreshold, 0);

% decomposiing a

entropythreshold=ent_a;
mindim=md;%2;
maxdim=min(size(a, 1))/2;

% if ent_a<ent_b
%     im=a;
% else
%     im=b;
% end

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

errmap=zeros(size(a));

for dim=dims
    clear blks;
    [ablks, r, c]=qtgetblk(a, q, dim);
    [bblks, r, c]=qtgetblk(b, q, dim);
    siz=size(ablks);
    numsubblks=size(ablks, 3);
    k=0;

    if ~numsubblks; continue; end;
    
    w=1;%/dim;
    totalweight=totalweight+numsubblks*w;

    ablks=ut_prepblks4bulkproc(ablks);
    bblks=ut_prepblks4bulkproc(bblks);
    step=500;
    eblks=[];
    for i=1:step:numsubblks

        si=i;
        ei=min(i+step-1, numsubblks);
        
        G=256;

        ha=hist(ablks(:, si:ei), 0:G-1); ha=ha/dim/dim;
        hb=hist(bblks(:, si:ei), 0:G-1); hb=hb/dim/dim;
        hab=ut_jhist(ablks(:, si:ei), bblks(:, si:ei));  hab=hab/dim/dim;

        if numsubblks==1; ha=ha'; hb=hb'; end;

        ent_as=-sum(log2(ha.^ha), 1);
        ent_bs=-sum(log2(hb.^hb), 1);
        ent_abs=-sum(log2(hab.^hab), 1);
        e=ent_as+ent_bs-ent_abs;
        nume=2*e;
        deno=ent_as+ent_bs;

        e=round(nume./deno*1000)/1000;

        err=e*w;
        eblks=[eblks, repmat(e, numel(ablks(:, 1)), 1)];
        totalerr=totalerr+sum(err);
    end
    
    eblks=reshape(eblks, siz);
    errmap=qtsetblk(errmap, q, dim+k, eblks);
end

err=totalerr/totalweight;

errmap_a=errmap;
err_a=err;


% decomposing b

entropythreshold=ent_b;
im=b;
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

errmap=zeros(size(a));

for dim=dims
    clear blks;
    [ablks, r, c]=qtgetblk(a, q, dim);
    [bblks, r, c]=qtgetblk(b, q, dim);
    siz=size(ablks);
    numsubblks=size(ablks, 3);
    k=0;

    if ~numsubblks; continue; end;
    
    w=1;%/dim;
    totalweight=totalweight+numsubblks*w;

    ablks=ut_prepblks4bulkproc(ablks);
    bblks=ut_prepblks4bulkproc(bblks);
    step=500;
    eblks=[];
    for i=1:step:numsubblks

        si=i;
        ei=min(i+step-1, numsubblks);
        
        G=256;

        ha=hist(ablks(:, si:ei), 0:G-1); ha=ha/dim/dim;
        hb=hist(bblks(:, si:ei), 0:G-1); hb=hb/dim/dim;
        hab=ut_jhist(ablks(:, si:ei), bblks(:, si:ei));  hab=hab/dim/dim;

        if numsubblks==1; ha=ha'; hb=hb'; end;

        ent_as=-sum(log2(ha.^ha), 1);
        ent_bs=-sum(log2(hb.^hb), 1);
        ent_abs=-sum(log2(hab.^hab), 1);
        e=ent_as+ent_bs-ent_abs;
        nume=2*e;
        deno=ent_as+ent_bs;

        e=round(nume./deno*1000)/1000;

        err=e*w;
        eblks=[eblks, repmat(e, numel(ablks(:, 1)), 1)];
        totalerr=totalerr+sum(err);
    end
    
    eblks=reshape(eblks, siz);
    errmap=qtsetblk(errmap, q, dim+k, eblks);
end

err=totalerr/totalweight;


errmap_b=errmap;
err_b=err;


errmap=errmap_a+errmap_b;
errmap=errmap/2;

err=err_a+err_b;
err=err/2;


