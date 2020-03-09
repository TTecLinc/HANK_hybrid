function y = foc_l( c,l,e )
global beta g_A tr r delta_h Psi Xi w tau 
global ih ia V_agrid it
global a_grid h_grid 
    a1=1/(1+g_A).*((a_grid(ia,ih)+tr).*(1+r)+l.*w.*h_grid(ia,ih).*(1-tau)-c);
    
    h1=(h_grid(ia,ih).*e).^(Psi).*Xi+h_grid(ia,ih).*(1-delta_h);
    V=squeeze(V_agrid(it+1,:,:));
    y=u_le(c,l,e)-beta.*h_grid(ia,ih).*w.*(1-tau)/(1+g_A)*griddata(a_grid,h_grid,V,a1,h1,'v4');
end