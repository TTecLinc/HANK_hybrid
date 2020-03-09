function y = foc_c( c,l,e )
global beta g_A tr r pen delta_h Psi Xi w tau tw
global it ih ia V_agrid phi
global a_grid h_grid 
if it>tw;
    a1=1/(1+g_A).*((a_grid(ia)+tr).*(1+r)+pen.*h_grid(ih)-c);
    
    h1=h_grid(ih).*(1-delta_h);
    %y=uc(c,0,0)-beta/(1+g_A).*inter_avalue(a1,h1);
    V=squeeze(V_agrid(it+1,:,:));
    y=uc(c,0,0)-beta/(1+g_A).*griddata(a_grid,h_grid,V,a1,h1,'v4');
else
    a1=1/(1+g_A).*((a_grid(ia,ih)+tr).*(1+r)+l.*w.*h_grid(ia,ih).*(1-tau)-c);
    
    h1=(h_grid(ia,ih).*e).^(Psi).*Xi+h_grid(ia,ih).*(1-delta_h);
    V=squeeze(V_agrid(it+1,:,:));
    y=uc(c,l,e)-beta/(1+g_A).*griddata(a_grid,h_grid,V,a1,h1,'v4');
end
end


