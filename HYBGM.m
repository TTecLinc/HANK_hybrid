%Parameter
clear;

global V_grid V_agird V_hgird ns nh na
global theta rho r delta alpha gamma w phi beta 
global s_grid h_grid a_grid 
global is ih it
global c_grid e_grid 

theta=0.5;
rho=0.04;
r=0.05;
delta=0.05;
alpha=0.35;
gamma=1;
w=0.1;
phi=0.5;
beta=0.95;
nh=30;
ns=30;
na=30;
nT=5; 

%Step1: Discrete Space
h_grid=linspace(0.01,3,nh);
s_grid=linspace(0.01,10,ns);
a_grid=linspace(0.01,10,na);
a_grid=repmat(a_grid,nh,1);
a_grid=a_grid';
c_grid=zeros(nT,ns,nh);
e_grid=zeros(nT,ns,nh);
h1_grid=zeros(nT,ns,nh);
a1_grid=zeros(nT,ns,nh);
V_grid=zeros(nT,ns,nh);
V_agird=zeros(nT,ns,nh);
V_hgird=zeros(nT,ns,nh);


c_profile=zeros(nT,1);
e_profile=zeros(nT,1);
a_profile=zeros(nT,1);
h_profile=zeros(nT,1);
%Step2: The last Period
for ia=1:na;
    for ih=1:nh;
        c_grid(nT,ia,ih)=a_grid(ia,ih)+w.*h_grid(ih);
        
        V_grid(nT,ia,ih)=c_grid(nT,ia,ih).^(1-theta)/(1-theta);
        V_agird(nT,ia,ih)=c_grid(nT,ia,ih).^(-theta);
        V_hgird(nT,ia,ih)=w.*c_grid(nT,ia,ih).^(-theta);
    end
end 
%Step3: Iterate Backwards
for it=nT-1:-1:1;
    for is=1:ns;
        for ih=1:nh;   
            at1=(1+r).*s_grid(is);
            solutione=fsolve(@(x)oned_optequation(x),0.5);
            e_grid(it,is,ih)=solutione;
            %(c) compute c
            h1_grid(it,is,ih)=(1-delta).*(h_grid(ih)+gamma/alpha.*e_grid(it,is,ih).^alpha);
            c_grid(it,is,ih)=(beta.*(1+r).*(1-phi/(1+h1_grid(it,is,ih))).*hyinter_avalue(at1,h1_grid(it,is,ih))).^(-1/theta);   
            a1_grid(it,is,ih)=s_grid(is)-w.*h_grid(ih)+c_grid(it,is,ih)+e_grid(it,is,ih);
            V_grid(it,is,ih)=1/(1-theta).*c_grid(it,is,ih).^(1-theta)+beta.*(1-phi/(1+h1_grid(it,is,ih))).*hyinter_value(at1,h1_grid(it,is,ih));
            V_agird(it,is,ih)=c_grid(it,is,ih).^(-theta);
            V_hgird(it,is,ih)=(w+1/gamma.*e_grid(it,is,ih).^(1-alpha)).*c_grid(it,is,ih).^(-theta);            
        end
    end
    a_grid=a1_grid(it,:,:);  
    a_grid=squeeze(a_grid);
end

for i=1:nT;
    c_profile(i)=mean(mean(c_grid(i,:,:)));
    e_profile(i)=mean(mean(e_grid(i,:,:)));
    a_profile(i)=mean(mean(a1_grid(i,:,:)));
    h_profile(i)=mean(mean(h1_grid(i,:,:)));
end

plot(c_profile);
title('Consumption Profile');

figure();
plot(e_profile);
title('Education Profile');

figure();
plot(a_profile);
title('Asset Profile');

figure();
plot(h_profile);
title('Human Capital Profile');
