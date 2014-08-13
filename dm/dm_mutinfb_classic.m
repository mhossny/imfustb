function [b_map, d_map]=dm_mutinfb_classic(a, b, w)
% Localised mutual information metric as modified by Hossny et al. in
% 
% Hossny, M.; Nahavandi, S.; Creighton, D.; Bhatti, A, "Image fusion performance
% metric based on mutual information and entropy driven quadtree decomposition," 
% Electronics Letters , vol.46, no.18, pp.1266,1268, September 2010
% doi: 10.1049/el.2010.1778
% 
% The MI is un-normalised and biased towords the source image with the
% highest entropy as explained by Hossny et al. in 
% 
% Hossny, M.; Nahavandi, S.; Creighton, D., "Comments on 'Information measure 
% for performance of image fusion'," Electronics Letters , vol.44, no.18, 
% pp.1066,1067, August 28 2008 (doi: 10.1049/el:20081754)
% 
% Get rid of the progress bar for quick batch processing. Also maximise
% step to benefit from block processing. You may run out of memory
% depending on the specs of your machinand the size of the images. 

ws=ut_itknl(w);

[wm, wn, wd]=size(ws);
[am, an, ad]=size(a);
[bm, bn, bd]=size(b);


a_sampling_maps=(ut_cnvit(a, ws));
a_sampling_maps=reshape(a_sampling_maps, [prod(size(a_sampling_maps))/wd, wd])';
b_sampling_maps=(ut_cnvit(b, ws));
b_sampling_maps=reshape(b_sampling_maps, [prod(size(b_sampling_maps))/wd, wd])';
%a_sampling_maps=hist(a_sampling_maps, 255)/wd;
%a_sampling_maps=a_sampling_maps+.000000001;%(find(a_histograms>0));

G=256;

step=15000;
es=repmat(0, size(a_sampling_maps, 2), 1);

hprg=waitbar(0, 'calculating entropy map...');
for i=1:step:size(a_sampling_maps, 2)
    if i+step>size(a_sampling_maps, 2)
        step=size(a_sampling_maps, 2)-i;
    end
	ha=hist(double(a_sampling_maps(:, i:i+step)), G)/wd;
	ha=ha;%(find(a_histograms>0));
	hb=hist(double(b_sampling_maps(:, i:i+step)), G)/wd;
	hb=hb;%(find(a_histograms>0));
	%hab=hist(double([a_sampling_maps(:, i:i+step);b_sampling_maps(:, i:i+step)]), 256)/(2*wd);
	[hab]=ut_jhist(a_sampling_maps(:, i:i+step), b_sampling_maps(:, i:i+step));
	%[hab]=hist(double([a_sampling_maps(:, i:i+step), b_sampling_maps(:, i:i+step)]), 256);
    hab=hab/(1*wd);

    ent_a=-sum(log2(ha.^ha), 1);
	ent_b=-sum(log2(hb.^hb), 1);
    ent_ab=-sum(log2(hab.^hab), 1);

    e=ent_a+ent_b-ent_ab;
    

    es(i:i+step)=e;
    waitbar(i/size(a_sampling_maps, 2), hprg);
	clear h;
end
close(hprg);

b_map=reshape(es, [am, an]);
[b_map, d_map]=ut_prpmp(b_map, hot(256), size(a));
% a_sampling_maps1=a_sampling_maps(:, 1:size(a_sampling_maps, 2)/4);
% a_sampling_maps1=-a_sampling_maps1'*log2(a_sampling_maps1);
% 
% a_sampling_maps2=a_sampling_maps(:, size(a_sampling_maps, 2)/2+1:end);
% a_sampling_maps2=-a_sampling_maps2'*log2(a_sampling_maps2);
% 
% a_sampling_maps=[a_sampling_maps1, a_sampling_maps2];
% 
%a_sampling_maps1=-a_sampling_maps(:, 1:25000)'*log2(a_sampling_maps(:, 1:25000));

% a_sampling_maps=conv_iterator(a, ws);
% a_sampling_vectors=reshape(a_sampling_maps, [prod(size(a_sampling_maps))/wd, wd])';
% a_histograms=hist(a_sampling_vectors, 0:255)/wd;
% a_histograms=a_histograms+.000000001;%(find(a_histograms>0));
% a_entropy_map=-a_histograms*log2(a_histograms');
