function y = hyfoc_c( c )
global beta g_A tr r pen delta_h Psi Xi w tau tw
global it ih ia is 
global na nh phi ex_n
global a_grid h_grid s_grid V_agrid e_grid
if it>tw;
    %a1=1/(1+g_A).*((a_grid(ia)+tr).*(1+r)+pen-c);
    a1=s_grid(is).*(1+g_A);
    h1=h_grid(ih).*(1-delta_h);
    V=squeeze(V_agrid(it+1,:,:));
    %'v4' 'natural'
    y=uc(c,0,0)-beta/(1+g_A).*griddata(a_grid,h_grid,V,a1,h1,'v4');
else
    a1=s_grid(is).*(1+g_A);
    h1=(h_grid(ia,ih).*e_grid(it,is+ex_n,ih)).^(Psi).*Xi+h_grid(ia,ih).*(1-delta_h);
    l=1-(1-phi)/phi/h_grid(ia,ih)/w/(1-tau).*c-e_grid(it,is+ex_n,ih);
    V=squeeze(V_agrid(it+1,:,:));
    y=uc(c,l,e_grid(it,is+ex_n,ih))-beta/(1+g_A).*griddata(a_grid,h_grid,V,a1,h1,'v4');
end
