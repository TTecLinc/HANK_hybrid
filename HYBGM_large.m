clear;
global sigma phi 
global beta g_A tr r pen delta_h Psi Xi w tau tre tw
global it ih ia is
global a_grid h_grid s_grid e_grid
global V_agrid V_hgrid na nh ns
global max_a max_h min_a min_h ex_n
sigma=2;
%beta=1.021;
phi=0.370;

%Human Capital
Xi=0.176;
%Psi=0.576;
Psi=0.076;
delta_h=1.4/100;
h_0=1;
%Production
alpha=0.33;
delta=3.5/100;
tau=0.5;

r=0.07;
w=1.5;
g_A=1.9/100;

beta=0.985;
beta=beta.*(1+g_A).^(phi.*(1-sigma));
% The number of Grid
tw=20;
tre=10;
nT=tw+tre;
min_a=0;
min_h=0.5;
max_a=5;
max_h=2.5;

na=25;
nh=na;
ns=na;
ex_n=na-ns;
%Exogenous Social Security
pen=0.15;
tr=0.16;

h_grid=linspace(min_h,max_h,nh);
h_grid=repmat(h_grid,na,1);
a_grid=linspace(min_a,max_a,na);
a_grid=repmat(a_grid,nh,1);
a_grid=a_grid';
s_grid=linspace(min_a+0.2,max_a-1.5,na);

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


%% The last Period
for ia=1:na;
    for ih=1:nh;
        c_grid(nT,ia,ih)=(a_grid(ia,ih)+tr).*(1+r)+pen.*h_grid(ia,ih);
        %Iteration the Grid
        V_grid(nT,ia,ih)=u(c_grid(nT,ia,ih),0,0);
        V_agrid(nT,ia,ih)=uc(c_grid(nT,ia,ih),0,0).*(1+r);
        V_hgrid(nT,ia,ih)=uc(c_grid(nT,ia,ih),0,0).*pen;
    end
end 
%%
% nT-1,...,tw+1

for it=nT-1:-1:tw+1;
    for is=1:ns;
        for ih=1:nh;
            c_grid(it,is+ex_n,ih)=fsolve(@(x)hyfoc_c(x),c_grid(it+1,is+ex_n,ih));
            
            %
            a1_grid(it,is+ex_n,ih)=(s_grid(is).*(1+g_A)-pen.*h_grid(is+ex_n,ih)+c_grid(it,is+ex_n,ih))/(1+r)-tr;
            h1_grid(it,is+ex_n,ih)=h_grid(is+ex_n,ih).*(1-delta_h);
            V_agrid(it,is+ex_n,ih)=(1+r).*uc(c_grid(it,is+ex_n,ih),0,0);
            V=squeeze(V_hgrid(it+1,:,:));
            
            V_hgrid(it,is+ex_n,ih)=beta.*(1-delta_h).*griddata(a_grid,h_grid,V,a1_grid(it,is+ex_n,ih),h1_grid(it,is+ex_n,ih),'v4')...
                +h_grid(is+ex_n,ih).*uc(c_grid(it,is+ex_n,ih),0,0);
        end
    end
    
    % The binding positions
    a_grid=a1_grid(it,:,:);
    a_grid=squeeze(a_grid);
    
%     if a_grid(ex_n+1,1)~=0;
%         for i_aux=1:nh;
%             a_grid(1:ex_n+1,i_aux)=linspace(0,a_grid(ex_n+2,i_aux).*0.6,ex_n+1);
%         end
% 
%         for ia=1:ex_n+1;
%             for ih=1:nh;
%                 c_grid(it,ia,ih)=(a_grid(ia,ih)+tr).*(1+r)+pen.*h_grid(ia,ih);
%                 a1_grid(it,ia,ih)=a_grid(ia,ih);
%                 h1_grid(it,ia,ih)=h_grid(ia,ih).*(1-delta_h);
%                 V_agrid(it,ia,ih)=(1+r).*uc(c_grid(it,ia,ih),0,0);
%                 V_hgrid(it,ia,ih)=beta.*(1-delta_h).*griddata(a_grid,h_grid,V,a1_grid(it,ia,ih),h1_grid(it,ia,ih),'v4');
%             end
%         end
%     end
end

%%
for it=tw:-1:1;
    for is=1:ns;
        for ih=1:nh;
            if hyfoc_e(1e-4)>0;
                e_grid(it,is+ex_n,ih)=0;
            elseif hyfoc_e(1)<0;
                e_grid(it,is+ex_n,ih)=0.8;
            else
                e_grid(it,is+ex_n,ih)=fsolve(@(x)hyfoc_e(x),0.001);
            end
            
            
          
            c_grid(it,is+ex_n,ih)=fsolve(@(x)hyfoc_c(x),0.05);
            l_grid(it,is+ex_n,ih)=1-(1-phi)/phi/h_grid(ia,ih)/w/(1-tau).*c_grid(it,is+ex_n,ih)-e_grid(it,is+ex_n,ih);
            if l_grid(it,is+ex_n,ih)<0;
                l_grid(it,is+ex_n,ih)=0;
                c_grid(it,is+ex_n,ih)=(1-e_grid(it,is+ex_n,ih))/((1-phi)/phi/h_grid(ia,ih)/w/(1-tau));
            end
                %pause;
            a1_grid(it,is+ex_n,ih)=((1+g_A).*s_grid(is)+c_grid(it,is+ex_n,ih)-l_grid(it,is+ex_n,ih).*h_grid(ih).*w.*(1-tau))/(1+r)-tr;
            h1_grid(it,is+ex_n,ih)=h_grid(ia,ih).*(1-delta_h)+(h_grid(ia,ih).*e_grid(it,is+ex_n,ih)).^Psi.*Xi;
            V_agrid(it,is+ex_n,ih)=(1+r).*uc(c_grid(it,is+ex_n,ih),l_grid(it,is+ex_n,ih),e_grid(it,is+ex_n,ih));
            g_h=(1-delta_h)+Xi.*Psi.*(e_grid(it,is+ex_n,ih).*h_grid(is+ex_n,ih)).^(Psi-1).*e_grid(it,is+ex_n,ih);
            g_e=Xi.*Psi.*h_grid(ia,ih).*(e_grid(it,is+ex_n,ih).*h_grid(is+ex_n,ih)).^(Psi-1);
            V_hgrid(it,is+ex_n,ih)=abs(beta.*(l_grid(it,is+ex_n,ih).*w.*(1-tau).*uc(c_grid(it,is+ex_n,ih),l_grid(it,is+ex_n,ih),e_grid(it,is+ex_n,ih))/beta+...                
               g_h.*u_le(c_grid(it,is+ex_n,ih),l_grid(it,is+ex_n,ih),e_grid(it,is+ex_n,ih))/beta/g_e));
           if isnan(V_hgrid(it,is+ex_n,ih));
               V_hgrid(it,is+ex_n,ih)=500;
           end
        end
    end
    
    % The binding Positions
    a_grid=a1_grid(it,:,:);  
    a_grid=squeeze(a_grid);
%     if a_grid(ex_n+1,1)~=0;
%         for i_aux=1:nh;
%             a_grid(1:ex_n+1,i_aux)=linspace(0,a_grid(ex_n+2,i_aux).*0.2,ex_n+1);
%         end
% 
%         for ia=1:ex_n+1;
%             for ih=1:nh;
%                 e_grid(it,is+ex_n,ih)=e_grid(it,ex_n+ia,ih).*1.1;
%                 l_grid(it,is+ex_n,ih)=1-e_grid(it,is+ex_n,ih);
%                 c_grid(it,is+ex_n,ih)=(a_grid(ia,ih)+tr).*(1+r)+l_grid(it,is+ex_n,ih).*h_grid(ia,ih).*w.*(1-tau);
%                 a1_grid(it,ia,ih)=a_grid(ia,ih);
%                 h1_grid(it,ia,ih)=h_grid(ia,ih).*(1-delta_h);
%                 V_agrid(it,ia,ih)=(1+r).*uc(c_grid(it,ia,ih),l_grid(it,is+ex_n,ih),e_grid(it,is+ex_n,ih));
%                 V_hgrid(it,ia,ih)=beta.*(1-delta_h).*griddata(a_grid,h_grid,V,a1_grid(it,ia,ih),h1_grid(it,ia,ih),'v4');
%             end
%         end
%     end
end

for i=1:nT;
    c_profile(i)=mean(mean(c_grid(i,:,:)));
    e_profile(i)=mean(mean(e_grid(i,:,:)));
    l_profile(i)=mean(mean(l_grid(i,:,:)));
    a_profile(i)=mean(mean(a1_grid(i,:,:)));
    h_profile(i)=mean(mean(h1_grid(i,:,:)));
end

plot(c_profile);
title('Consumption Profile');

figure();
plot(e_profile);
title('Education Profile');

figure();
plot(l_profile);
title('Labor Supply Profile');

figure();
plot(h_profile);
title('Human Capital Profile');

figure();
plot(a_profile);
title('Asset Profile');
