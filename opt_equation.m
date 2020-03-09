function [x,y]=opt_equation(c,e)
global V_grid V_agird V_hgird na nh 
global theta rho r delta alpha gamma w phi beta 
global a_grid h_grid
global ia ih it

at1=(1+r).*(a_grid(ia)+w.*h_grid(ih)-c-e);
ht1=(1-delta).*(h_grid(ih)+gamma/alpha.*e.^alpha);
V_t1a=inter_avalue(at1,ht1);
V_t1h=inter_hvalue(at1,ht1);
V_t1=inter_value(at1,ht1);

x=beta.*(1+r).*(1-phi/(1+(1-delta).*(h_grid(ih)+gamma/alpha.*(e).^alpha))).*V_t1a-c.^(-theta);
y=(1+r)/(1-delta).*V_t1a/(phi/(1+ht1-phi)/(1+ht1).*V_t1+V_t1h);
end