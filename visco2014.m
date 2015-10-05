clear all
close all
addpath('../FEM_beta')

eta=0.5;
eps=0;
load data/mesh %charge le mesh

m.disp

%parametres physiques
%sigma_real=m.tensor('5+exp(-2*(x.^2+ 2*y.^2) ).*( 0.3*sin(10*x.*y)+ x.^2+x.*y.^2 +5*exp(-10*((x-3/5).^2+(y+0.02).^2)))',0,0);
sigma_real=m.tensor('5+( 1/4* x.^2+x.*y.^2 +5*exp(-10*((x-3/5).^2+(y+0.02).^2)))',0,0);
%sigma_real=m.tensor('1+4*(((x-1/3).^2+(y-0.1).^2)<0.01)+8*(((x-2/3).^2+(y+0.1).^2)<0.01)',0,0);
%sigma_real=m.tensor('4+2*exp(-80*((x-3/5).^2+(y+0.02).^2))+0.5*sin(25*x.*y.^2+x.^2-18*(x+y))',0,0); %la vrai conductivité
%sigma_real=m.tensor('3+4*exp(-25*((x-3/5).^2+(y+0.02).^2))',0,0);
sigma_real_t=m.nod2tri(sigma_real);
m.surf(sigma_real)
%%
Ax=m.tensor('1/100 *(1/2*y+1)',0,0);
Ay=m.tensor('1/100 *(-1/2*x+1)',0,0);
A=[Ax;Ay]; % le potentiel A
Ax_t=m.tensor('1/100 *(1/2*y+1)',0,1);
Ay_t=m.tensor('1/100 *(-1/2*x+1)',0,1);
A_t=[Ax_t;Ay_t];


bc=BoundaryCondition([1 2],1,0,'1/100*( -1/4* x.*y - x/2   +  x.*y -  2*y)./sqrt(1/4*x.^2 + 2*y.^2)'); %expression correspondant a A scalaire nu
source_r=[sigma_real.*A(1,:);sigma_real.*A(2,:)];
div_source_r=m.div(source_r);
V_r=m.solvepde(bc,sigma_real,0,div_source_r);
g_V_r=m.grad(V_r);


J=[sigma_real_t.*(g_V_r(1,:)+A_t(1,:));sigma_real_t.*(g_V_r(2,:)+A_t(2,:))];
F=[-J(2,:);J(1,:)];

M_f=[F(1,:).^2; F(1,:).*F(2,:);F(1,:).*F(2,:); F(2,:).^2];
M_f_approx=[F(1,:).^2 + eta ; F(1,:).*F(2,:);F(1,:).*F(2,:); F(2,:).^2+eta];



FFA_t=[M_f(1,:).*A_t(1,:)+ M_f(2,:).*A_t(2,:);M_f(3,:).*A_t(1,:) + M_f(4,:).*A_t(2,:)];
FFA=m.tri2nod(FFA_t);

div_FFA=m.div(FFA);



V_approx=m.solvepde(bc,M_f_approx,0,-div_FFA);

g_V_approx=m.grad(V_approx);

r_t=sqrt( (g_V_approx(1,:)+ A_t(1,:)).^2 + (g_V_approx(2,:)+ A_t(2,:)).^2) ./ sqrt( J(1,:).^2 + J(2,:).^2+eps);

sigma_try_t=  sqrt( J(1,:).^2 + J(2,:).^2+eps) ./  sqrt((g_V_approx(1,:)+ A_t(1,:)).^2 + (g_V_approx(2,:)+ A_t(2,:)).^2);
sigma_try=m.tri2nod(sigma_try_t);


sigma_t=1./r_t;

sigma=m.tri2nod(sigma_t);

figure
subplot(211)
m.surf(V_r)
subplot(212)
m.surf(V_approx)

figure
subplot(311)
m.surf(sigma_real)
subplot(312)
m.surf(sigma)
subplot(313)
m.surf(sigma_try)
