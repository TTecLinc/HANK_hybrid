global sigma phi 
global beta g_A tr r pen delta_h Psi Xi w tau tre tw
global it ih ia
global a_grid h_grid
global V_agrid na nh V_hgrid
global max_a max_h min_a min_h
sigma=2;
%beta=1.021;
phi=0.370;
%Human Capital
Xi=0.176;
Psi=0.576;
delta_h=1.4/100;
h_0=1;
%Production
alpha=0.33;
delta=3.5/100;
tau=0.5;
omega=0.7;
r=0.07;
w=0.1;
g_A=1.9/100;
beta=0.985;
beta=beta.*(1+g_A).^(phi.*(1-sigma));
% The number of Grid
tw=7;
tre=2;
nT=tw+tre;
min_a=0.1;
min_h=1;
max_a=6;
max_h=3;
na=15;
nh=15;
%Exogenous Social Security
pen=0.15;
tr=0.16;

%h_grid=linspace(min_h,max_h,nh);
%a_grid=linspace(min_a,max_a,na);
h_grid=linspace(min_h,max_h,nh);
h_grid=repmat(h_grid,na,1);
a_grid=linspace(min_a,max_a,na);
a_grid=repmat(a_grid,nh,1);
a_grid=a_grid';

c_profile=zeros(1,nT);
l_profile=zeros(1,tw);
e_profile=zeros(1,tw);
h_profile=zeros(1,nT);
a_profile=zeros(1,nT);

V_agrid=zeros(nT,na,nh);
V_hgrid=zeros(nT,na,nh);
V_grid=zeros(nT,na,nh);

c_grid=zeros(nT,na,nh);
e_grid=zeros(nT,na,nh);
l_grid=zeros(nT,na,nh);
a1_grid=zeros(nT,na,nh);
h1_grid=zeros(nT,na,nh);


% The last Period
for ia=1:na;
    for ih=1:nh;
        c_grid(nT,ia,ih)=(a_grid(ia,ih)+tr).*(1+r)+pen.*h_grid(ia,ih);
        V_grid(nT,ia,ih)=u(c_grid(nT,ia,ih),0,0);
        V_agrid(nT,ia,ih)=uc(c_grid(nT,ia,ih),0,0).*(1+r);
        V_hgrid(nT,ia,ih)=uc(c_grid(nT,ia,ih),0,0).*pen;
    end
end 
%From nT-1,...tw+1
for it=nT-1:-1:tw+1;
    for ia=1:na;
        for ih=1:nh;
            
            options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','MaxFunEvals',1000);
            solution=fsolve(@(x)foc_c(x,0,0),c_grid(it+1,ia,ih),options);
            
            c_grid(it,ia,ih)=solution;
            a1_grid(it+1,ia,ih)=(a_grid(ia,ih)+tr).*(1+r)+pen.*h_grid(ia,ih);
            h1_grid(it+1,ia,ih)=(h_grid(ia,ih)).*(1-delta_h);
            %Iteration the grid
            V=squeeze(V_hgrid(it+1,:,:));
            V_agrid(it,ia,ih)=(1+r).*uc(c_grid(it,ia,ih),0,0);
            V_hgrid(it,ia,ih)=beta.*(1-delta_h).*griddata(a_grid,h_grid,V,a1_grid(it,ia,ih),h1_grid(it,ia,ih),'v4')...
                +h_grid(ia,ih).*uc(c_grid(it,ia,ih),0,0);
        end
    end
end
for it=tw:-1:1;
    for ia=1:na;
        for ih=1:nh;
            options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','MaxFunEvals',2000);
            solution=fsolve(@(x)opt_equs(x(1),x(2),(3)),[c_grid(it+1,ia,ih),l_grid(it+1,ia,ih),e_grid(it+1,ia,ih)],options);
            %f=@(x)[foc_c(x(1),x(2),x(3));foc_l(x(1),x(2),x(3));foc_e(x(1),x(2),x(3))];
            %solution=broyden(f,[0.5;0.6;0.5]);
            c_grid(it,ia,ih)=solution(1);
            l_grid(it,ia,ih)=solution(2);
            e_grid(it,ia,ih)=solution(3);
            pause;
            a1_grid(it+1,ia,ih)=1/(1+g_A).*((a_grid(ia)+tr).*(1+r)+l_grid(it,ia,ih).*w.*h_grid(ih).*(1-tau)-c_grid(it,ia,ih));
            h1_grid(it+1,ia,ih)=(h_grid(ih).* e_grid(it,ia,ih)).^(Psi).*Xi+h_grid(ih).*(1-delta_h);
            V_agrid(it,ia,ih)=(1+r).*uc(c_grid(it,ia,ih),l_grid(it,ia,ih),e_grid(it,ia,ih));
            g_h=(1-delta_h)+Xi.*Psi.*(e_grid(it,ia,ih).*h_grid(ih)).^(Psi-1).*e_grid(it,ia,ih);
            g_e=Xi.*Psi.*h_grid(ih).*(e_grid(it,ia,ih).*h_grid(ih)).^(Psi-1);
            V_hgrid(it,ia,ih)=beta.*(l_grid(it,ia,ih).*w.*(1-tau).*uc(c_grid(it,ia,ih),l_grid(it,ia,ih),e_grid(it,ia,ih))/beta+...                
                g_h.*u_le(c_grid(it,ia,ih),l_grid(it,ia,ih),e_grid(it,ia,ih))/beta/g_e);
        end
    end
end
% 
% for i=1:nT;
%     c_profile(i)=mean(mean(c_grid(i,:,:)));
%     e_profile(i)=mean(mean(e_grid(i,:,:)));
%     l_profile(i)=mean(mean(l_grid(i,:,:)));
%     a_profile(i)=mean(mean(a1_grid(i,:,:)));
%     h_profile(i)=mean(mean(h1_grid(i,:,:)));
% end
% 
% plot(c_profile);
% title('Consumption Profile');
% 
% figure();
% plot(e_profile);
% title('Education Profile');
% 
% figure();
% plot(l_profile);
% title('Labor Supply Profile');
% 
% figure();
% plot(h_profile);
% title('Human Capital Profile');
% 
% figure();
% plot(a_profile);
% title('Asset Profile');
