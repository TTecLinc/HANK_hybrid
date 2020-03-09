function y=oned_optequation(e)
global theta r delta alpha gamma w phi beta 
global s_grid h_grid
global is ih 

at1=(1+r).*s_grid(is);
ht1=(1-delta).*(h_grid(ih)+gamma/alpha.*e.^alpha);
V_t1a=hyinter_avalue(at1,ht1);
V_t1h=hyinter_hvalue(at1,ht1);
V_t1=hyinter_value(at1,ht1);
mul1=phi/(1+ht1-phi)/(1+ht1);
y=e-1/gamma.*((1+r)/(1-delta).*V_t1a/(mul1.*V_t1+V_t1h)).^(-1/(1-alpha));
end