function y=oned_optcquation(c)
global theta r delta alpha gamma w phi beta 
global s_grid h_grid e_grid
global is ih it

at1=(1+r).*s_grid(is);
ht1=(1-delta).*(h_grid(ih)+gamma/alpha.*e_grid(it,is,ih).^alpha);
y=c-(beta.*(1+r).*(1-phi/(1+ht1)).*hyinter_avalue(at1,ht1)).^(-1/theta);

end