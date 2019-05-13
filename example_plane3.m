% This example reconstructs the strain field and shows the results on
% cutting plane 2 as per the associated paper

clear all
clc

if isunix
    addpath('./inpolyhedron')
    addpath('./Beltrami_stress')
    addpath('./data')
else
    addpath('.\inpolyhedron')
    addpath('.\Beltrami_stress')
    addpath('.\data')
end

%% load sample info
sample = load('sample_coarse');
sample.V = sample.V*1e-3;       % convert to m from mm
sample.C = sample.C*1e-3;

%% define mesh

kow_pts = load('CP_points');

[X, Y] = meshgrid(linspace(-8.499e-3,8.499e-3,50),linspace(-8.5e-3,8.5e-3,50));
Z = 2.5e-3*ones(size(X));

IN = inpolyhedron(sample.F,sample.V,[X(:),Y(:),Z(:)],'FaceNormals',sample.N);


figure(1)
clf
p = patch('Faces',sample.F,'Vertices',sample.V);
p.FaceAlpha =0.1;
axis equal
hold on
h1 = plot3(X(IN),Y(IN),Z(IN),'bo');
hold off
xlabel('x')
ylabel('y')
zlabel('z')
view([-19.5000, 58])
title('Cutting plane 3')
legend([h1],'Cutting plane 3')


%% load strain image measurement data
load('strain_image_data_set')

num_m = length(y);

disp('Loading strain image data')
disp(['Number of measurements: ' num2str(num_m)])
fprintf('Average measurement standard deviation: %e \n', mean(y_std))


fprintf('\n')
disp('Reconstructing Strain field')

%% set up GP model parameters
num_basis = 15;
load('optimisation_results.mat')    % load hyperparemeters 
thetaopt = optim_theta;

A.lx = thetaopt(2);
A.ly = thetaopt(3);
A.lz = thetaopt(4);
A.sig_f = thetaopt(5);

B.lx = thetaopt(6);
B.ly = thetaopt(7);
B.lz = thetaopt(8);
B.sig_f = thetaopt(9);

C.lx = thetaopt(10);
C.ly = thetaopt(11);
C.lz = thetaopt(12);
C.sig_f = thetaopt(13);

D.lx = thetaopt(14);
D.ly = thetaopt(15);
D.lz = thetaopt(16);
D.sig_f = thetaopt(17);

E.lx = thetaopt(18);
E.ly = thetaopt(19);
E.lz = thetaopt(20);
E.sig_f = thetaopt(21);

F.lx = thetaopt(22);
F.ly = thetaopt(23);
F.lz = thetaopt(24);
F.sig_f = thetaopt(25);


%% form basis functions
nu = 0.28;      % poissons ratio
load('traction_points_20')  % load traction measurement point locations

Xcom = [X(IN);Xt(:)];
Ycom = [Y(IN);Yt(:)];
Zcom = [Z(IN);Zt(:)];

num_t = length(X(IN));
num_b = size(Xt,1); 

[Phi, SLambda,~,Phi_T] = beltrami_approx(num_basis,A,B,C,D,E,F,nu,entry,exit,Xcom,Ycom,Zcom,nsegs);

Phi_bound = Phi_T(num_t*6+1:num_t*6+num_b*6,:);       % extract boundary phis
Phi_T = Phi_T(1:num_t*6,:);

%% work out the traction basis functions
% hooke's law in 3d (scaled by 1/E)
HookesS = [1-nu, nu, nu, 0, 0, 0;
           nu, 1-nu, nu, 0, 0, 0;
           nu, nu, 1-nu, 0, 0, 0;
           0, 0, 0, 1-2*nu, 0, 0;
           0, 0, 0, 0, 1-2*nu, 0;
           0, 0, 0, 0, 0, 1-2*nu]/(1+nu)/(1-2*nu);
       
HCell = repmat({HookesS}, 1, length(Xt));
HH = blkdiag(HCell{:}); 
Phi_bstress = HH*Phi_bound;

clear nbCell;
for i = 1:length(Xt)
    nbCell{i} = [Nt(1,i), 0, 0, Nt(2,i), Nt(3,i), 0;
                0, Nt(2,i), 0, Nt(1,i), 0, Nt(3,i);
                0, 0, Nt(3,i), 0, Nt(1,i), Nt(2,i)];
end
NbNb = blkdiag(nbCell{:});
Phi_trac =NbNb*Phi_bstress;
y_trac = zeros(length(Xt)*3,1);
%% solve with tractions 
sig_t = 1e-8;       % traction variance

Phi_com = [Phi;Phi_trac];
y_com = [y;y_trac];

[n,m]= size(Phi);
[nt,~] = size(Phi_trac);


Gamma = [Phi_com./[y_std;ones(nt,1)*sig_t];diag(1./sqrt(SLambda))];
R = triu(qr(Gamma));
CZ = R(1:m,1:m);

optsT.TRANSA = true; optsT.UT = true; opts.TRANSA = false; opts.UT = true;
f_approx = Phi_T*(linsolve(CZ,linsolve(CZ,Phi_com'*(y_com./[y_std.^2;sig_t^2*ones(nt,1)]),optsT),opts));


epsxx_pred = NaN(size(X));
epsyy_pred = NaN(size(X));
epszz_pred = NaN(size(X));
epsxy_pred = NaN(size(X));
epsxz_pred = NaN(size(X));
epsyz_pred = NaN(size(X));

epsxx_pred(IN) = f_approx(1:6:end);
epsyy_pred(IN) = f_approx(2:6:end);
epszz_pred(IN) = f_approx(3:6:end);
epsxy_pred(IN) = f_approx(4:6:end);
epsxz_pred(IN) = f_approx(5:6:end);
epsyz_pred(IN) = f_approx(6:6:end);


absval_xx = max(abs([epsxx_pred(IN)]));
absval_yy = max(abs([epsyy_pred(IN)]));
absval_zz = max(abs([epszz_pred(IN)]));
absval_xy = max(abs([epsxy_pred(IN)]));
absval_xz = max(abs([epsxz_pred(IN)]));
absval_yz = max(abs([epsyz_pred(IN)]));


%% load the fea results
% coordinate system fix
% fea y is -y of ours
% and so fea exy is - of our exy
% same with eyz
xx_FEA =importdata('exx.txt','\t');
fxx = scatteredInterpolant(xx_FEA.data(:,2)*1e-3,-xx_FEA.data(:,3)*1e-3,xx_FEA.data(:,4)*1e-3,xx_FEA.data(:,5));

yy_FEA =importdata('eyy.txt','\t');
fyy = scatteredInterpolant(yy_FEA.data(:,2)*1e-3,-yy_FEA.data(:,3)*1e-3,yy_FEA.data(:,4)*1e-3,yy_FEA.data(:,5));

zz_FEA =importdata('ezz.txt','\t');
fzz = scatteredInterpolant(zz_FEA.data(:,2)*1e-3,-zz_FEA.data(:,3)*1e-3,zz_FEA.data(:,4)*1e-3,zz_FEA.data(:,5));

xy_FEA =importdata('gxy.txt','\t');
fxy = scatteredInterpolant(xy_FEA.data(:,2)*1e-3,-xy_FEA.data(:,3)*1e-3,xy_FEA.data(:,4)*1e-3,-xy_FEA.data(:,5)/2);

xz_FEA =importdata('gxz.txt','\t');
fxz = scatteredInterpolant(xz_FEA.data(:,2)*1e-3,-xz_FEA.data(:,3)*1e-3,xz_FEA.data(:,4)*1e-3,xz_FEA.data(:,5)/2);

yz_FEA =importdata('gyz.txt','\t');
fyz = scatteredInterpolant(yz_FEA.data(:,2)*1e-3,-yz_FEA.data(:,3)*1e-3,yz_FEA.data(:,4)*1e-3,-yz_FEA.data(:,5)/2);



eps_xx_fea = nan(size(X));
eps_yy_fea = nan(size(X));
eps_zz_fea = nan(size(X));
eps_xy_fea = nan(size(X));
eps_xz_fea = nan(size(X));
eps_yz_fea = nan(size(X));

eps_xx_fea(IN) = fxx(X(IN),Y(IN),Z(IN));
eps_yy_fea(IN) = fyy(X(IN),Y(IN),Z(IN));
eps_zz_fea(IN) = fzz(X(IN),Y(IN),Z(IN));
eps_xy_fea(IN) = fxy(X(IN),Y(IN),Z(IN));
eps_xz_fea(IN) = fxz(X(IN),Y(IN),Z(IN));
eps_yz_fea(IN) = fyz(X(IN),Y(IN),Z(IN));

%%
tmp = load('RdYlBu_up.mat');
cmap = tmp.RdYlBu;

FigHandle = figure(1);
clf

FigHandle.Position = [150   240   920   654];
colormap(cmap)

axis_width = 0.25;
axis_height = 0.25;
top = 0.65;
left = 0.1;
v_space = 0.035;
h_space = -0.04;
% gap = 0.05
fontsize = 25;
tmp = load('RdYlBu_up.mat');
cmap = tmp.RdYlBu;



% reconstruction results
subplot 331
tmp = flipud(fliplr(epsyy_pred));       % because the kowari and other arent lined up the same
pcolor(X,Y,tmp)
title('$\epsilon_{yy}$','Interpreter','latex','FontSize',fontsize*1.5)
shading interp
axis equal
caxis([-absval_yy absval_yy])
axis off
h1 = gca;
% h1.Position
h1.Position = [left  top    axis_width    axis_height];
t1 = text(-0.011,-0.00,'RADEN','Interpreter','latex','FontSize',fontsize);
set(t1,'Rotation',90,'HorizontalAlignment','center')


subplot 332
tmp = flipud(fliplr(epsxy_pred)); 
pcolor(X,Y,tmp)
shading interp
title('$\epsilon_{xy}$','Interpreter','latex','FontSize',fontsize*1.5)
axis equal
caxis([-absval_xy absval_xy])
axis off
h1 = gca;
h1.Position = [left+(axis_width+h_space)  top    axis_width    axis_height];

subplot(3,3,3)
tmp = flipud(fliplr(epsyz_pred)); 
pcolor(X,Y,tmp)
title('$\epsilon_{yz}$','Interpreter','latex','FontSize',fontsize*1.5)
shading interp
axis equal
caxis([-absval_yz absval_yz])
axis off
h1 = gca;
h1.Position = [left+2*(axis_width+h_space)  top    axis_width    axis_height];


% % FEA
subplot(3,3,4)
tmp = flipud(fliplr(eps_yy_fea)); 
pcolor(X,Y,tmp)
shading interp
axis equal
caxis([-absval_yy absval_yy])
axis off
colorbar('SouthOutside')
xlim([-8.5e-3 8.5e-3])
h1 = gca;
h1.Position = [left  top-1*(axis_height+v_space)   axis_width    axis_height];
t1 = text(-0.011,-0.00,'FEA','Interpreter','latex','FontSize',fontsize);
set(t1,'Rotation',90,'HorizontalAlignment','center')
% 
% 
subplot(3,3,5)
tmp = flipud(fliplr(eps_xy_fea)); 
pcolor(X,Y,tmp)
shading interp
axis equal
caxis([-absval_xy absval_xy])
axis off
colorbar('SouthOutside')
xlim([-8.5e-3 8.5e-3])
h1 = gca;
h1.Position = [left+1*(axis_width+h_space)  top-1*(axis_height+v_space)   axis_width    axis_height];

subplot(3,3,6)
tmp = flipud(fliplr(eps_yz_fea)); 
pcolor(X,Y,tmp)
shading interp
axis equal
caxis([-absval_yz absval_yz])
axis off
colorbar('SouthOutside')
xlim([-8.5e-3 8.5e-3])
h1 = gca;
h1.Position = [left+2*(axis_width+h_space)  top-1*(axis_height+v_space)   axis_width    axis_height];
%%
fprintf('\n')
disp('Computed difference betwee FEA results and reconstruction')
error = [epsxx_pred(IN);epsyy_pred(IN);epszz_pred(IN);epsxy_pred(IN);epsxz_pred(IN);epsyz_pred(IN)]-...
    [eps_xx_fea(IN);eps_yy_fea(IN);eps_zz_fea(IN);eps_xy_fea(IN);eps_xz_fea(IN);eps_yz_fea(IN)];
disp(['mean difference = ' num2str(mean(error))])
disp(['mean difference magnitude = ' num2str(mean(abs(error)))])
rel_error = mean(abs(error))/max(abs([eps_xx_fea(IN);eps_yy_fea(IN);eps_zz_fea(IN);eps_xy_fea(IN);eps_xz_fea(IN);eps_yz_fea(IN)]));
disp(['mean relative difference = ' num2str(rel_error)])

