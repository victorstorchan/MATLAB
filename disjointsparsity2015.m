%% 1. Parameters
%addpath(genpath('/import/maths/alberti/Dropbox/univ/matlab/2014-sparse/Toolboxes'))
clear
set(0,'DefaultFigureWindowStyle','docked')
n=128;
itermax=1500;

B={'x','y','x+y','x-y'};

%Ww= createwv(n,'Haar',4)
if n==64
    load('Ww.mat','Ww')
end
if n==128
    load('Wwbig.mat','Ww')
end
Wf=createfourier(n,n/8);


n2=n^2;
X=0:1/(n-1):1;
Y=X;

sigma=ones(n);


for i=1:n
    for j=1:n
       
        if (i-n/2)^2+(j-n/2)^2<n*2
            sigma(i,j)=2;
        end
        %if abs(i-n/3*2)+abs(j-n/3)<n/8
        %    mureal(i,j)=1.6;
        %end
        %if abs(i-n/1.4)<n/6 && abs(j-n/1.3)<n/18
        %    mureal(i,j)=1.5;
        %end
    end
end
suitereg = zeros(n);
for i=1:n
    for j=1:n
       
        if (i-n/2)^2+(j-n/2)^2<n*2
            suitereg(i,j)=exp(-1/(n*2-(i-n/2)^2-(j-n/2)^2));
        end
        
    end
end

%figure
%subplot(1,2,1)
%imagesc(sigma)
%title('gamma')
%axis equal
%colorbar
%caxis([0 2])
%subplot(1,2,2)
%imagesc(2*log(sigma))
%title('2*log(gamma)')
%axis equal
%colorbar



%% 2. Mesh

nodes=[0 1 1 0;0 0 1 1];
g=Geometry(nodes);
m=Mesh(g,0.012);
zone='(x-n/2).^2+(y-n/2).^2<(n/2)*1.4';
m=m.refine(g,zone);
m.plot

p=size(m.nodes,2);
t=size(m.triangles,2);
bc0=BoundaryCondition([1 2 3 4],0,1,0);

sigmat=grid2tri(m,X,Y,sigma);
sigmap=grid2nod(m,X,Y,sigma);

suiteregt=grid2tri(m,X,Y,suitereg);
suiteregp=grid2nod(m,X,Y,suitereg);


p=size(m.nodes,2);
t=size(m.triangles,2);


%% 3. Costruction of the data

NB=length(B);
UrealXY=B;

bc=B;
Ureal=B;



for l=1:NB
    bc{l}=BoundaryCondition([1 2 3 4],0,1,B{l});
    Ureal{l}=m.solvepde(bc{l},sigmat,0,0);
    UrealXY{l}=mesh2grid(m,X,Y,Ureal{l}');
    gu{l}=grad(m,Ureal{l});
    H{l}(1,:)=sigmat.*gu{l}(1,:);
    H{l}(2,:)=sigmat.*gu{l}(2,:);
    normeH{l}=H{l}(1,:).^2+H{l}(2,:).^2;
    
end

%% figures
figure
subplot(1,2,1)
m.surf(normeH{1})
colorbar
set(gca,'YDir','normal')
title('Real normeH')
%premier tau
[x,y] = meshgrid(0:1/(n-1):1,0:1/(n-1):1);
u = (.5-y)/.125+10*((x-.5).^2+(y-.5).^2-0.016).^10;
v = (x-.5)/.125+10*((x-.5).^2+(y-.5).^2-0.016).^10;
%u=(.5-y)/.18;
%v = (x-.5)/.18;


%deuxieme tau
[z,w] = meshgrid(0:1/(n-1):1,0:1/(n-1):1);
r = (.5-w)/0.125+0.0001*((z-.5).^2+(w-.5).^2-0.016).^8;
s =(z-.5)/0.125+0.0001*((z-.5).^2+(w-.5).^2-0.016).^8;

figure
subplot(1,2,1)
quiver(x(1:8:n,1:8:n),y(1:8:n,1:8:n),u(1:8:n,1:8:n),v(1:8:n,1:8:n))
title('values of \tau 1')
axis equal
subplot(1,2,2)
quiver(z(1:8:n,1:8:n),w(1:8:n,1:8:n),r(1:8:n,1:8:n),s(1:8:n,1:8:n))
title('values of \tau 2')
axis equal



%% smooth part

for l=1:NB
    gux{l}=tri2nod(m,gu{l}(1,:));
    guy{l}=tri2nod(m,gu{l}(2,:));
    guxXY{l}=mesh2grid(m,X,Y,gux{l}');
    guyXY{l}=mesh2grid(m,X,Y,guy{l}');
    suiteregx=tri2nod(m,suiteregt(1,:));
    %suiteregy=tri2nod(m,suiteregt(2,:));
    suiteregxXY=mesh2grid(m,X,Y,suiteregx');
    %suiteregyXY=mesh2grid(m,X,Y,suiteregy);
    Z{l}=guxXY{l}.*u+guyXY{l}.*v;
    Z2{l}=abs(Z{l}+awgn2(Z{l},11,'measured')).^2;
    Z3{l}=guxXY{l}.*r+guyXY{l}.*s;
    Z4{l}=abs(Z3{l}+awgn2(Z3{l},11,'measured')).^2;
    Z5{l}= ((guxXY{l}+awgn2(guxXY{l},11,'measured')).^2+(guyXY{l}+awgn2(guyXY{l},11,'measured')).^2).*suiteregxXY;
end
figure
imagesc(suiteregxXY)

figure
subplot(1,3,1)
imagesc(guxXY{l}.^2+guyXY{l}.^2)
colorbar
colormap(jet(300))
title('norm gradu')

subplot(1,3,2)
imagesc(Z2{1}+Z2{3}+Z5{2}.*1.5+Z4{2})

colorbar
title('sumreg1')
colormap(jet(300))

subplot(1,3,3)
imagesc(Z2{3}+Z2{2}+Z5{1}.*1.5+Z4{1})
colorbar
title('sumreg2')
colormap(jet(300))


%% separation

itermax=800;
NB=5;
somme{1}=2*log(sigma)+log(Z2{1}+Z2{3}+Z5{2}+Z4{4});
somme{2}=2*log(sigma)+log(Z2{3}+Z2{2}+Z5{1}+Z4{1});
somme{3}=2*log(sigma)+log(Z2{1}+Z2{4}+Z5{3}+Z4{2});
somme{4}=2*log(sigma)+log(Z2{1}+Z2{2}+Z5{4}+Z4{3});
somme{5}=2*log(sigma)+log(Z2{2}+Z2{3}+Z5{2}+Z4{1});

[irr,smoo]=separate(Wf,Ww,NB,somme,itermax,n);
figure
subplot(1,2,1)
imagesc(irr./2)
axis equal
colorbar
caxis([0 1])
subplot(1,2,2)
imagesc(log(sigma))
axis equal
colorbar
caxis([0 1])

errornorm(sigma, exp(irr./2))
