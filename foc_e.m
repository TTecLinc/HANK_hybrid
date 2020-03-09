function y = foc_e( c,l,e )
global beta g_A tr r delta_h Psi Xi w tau 
global ih ia it
global a_grid h_grid V_hgrid
    a1=1/(1+g_A).*((a_grid(ia,ih)+tr).*(1+r)+l.*w.*h_grid(ia,ih).*(1-tau)-c);
    
    h1=(h_grid(ia,ih).*e).^(Psi).*Xi+h_grid(ia,ih).*(1-delta_h);
    g_e=Xi.*Psi.*(e.*h_grid(ia,ih)).^(Psi-1).*h_grid(ia,ih);
    V=squeeze(V_hgrid(it+1,:,:));
    y=u_le(c,l,e)-beta.*g_e.*griddata(a_grid,h_grid,V,a1,h1,'v4');
end