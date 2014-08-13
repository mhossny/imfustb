function e=dm_mutinfTE(a, b, alpha)

a=double(a);
b=double(b);

% G=256;
wd=numel(a);

[hab, hahb]=ut_jhist(a(:), b(:));

hab=hab/(1*wd);
hahb=hahb/wd; sum(hahb)

ent_a=ut_entropyTE(a, alpha);%-sum(log2(ha.^ha), 1);
ent_b=ut_entropyTE(b, alpha);%-sum(log2(hb.^hb), 1);

ent_ab=(sum(hab.^alpha)-1)/(1-alpha);

e=ent_a+ent_b-ent_ab;


nume=2*e;
deno=ent_a+ent_b;

if nume~=deno
    e=nume/deno;
else
    e=1;
end

e=1./(1-alpha)*(1-sum(hab.^alpha./hahb.^(alpha-1), 1))