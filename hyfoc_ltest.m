function y = hyfoc_ltest( l )
global beta g_A tr r pen delta_h Psi Xi w tau tw
global it ih ia is 
global na nh phi
global a_grid h_grid s_grid V_agrid hgr

    %a1=1/(1+g_A).*((a_grid(ia)+tr).*(1+r)+pen-c);
    a1=s_grid(is).*(1+g_A);
    h1=h_grid(ih).*(1-delta_h);
    
    V=squeeze(V_agrid(it+1,:,:));
    
        c=(1-(1-phi)/phi/w/(1-tau)/h_grid(ih)/(1-l)).^(-1);
        y=uc(c,l,0)-beta/(1+g_A).*interp2(a_grid,hgr,V',a1,h1,'Spline');
    
end