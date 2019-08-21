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

[Y_before_rot, Z] = meshgrid(linspace(-12.0208e-3,12.0208e-3,50*round(sqrt(2))),linspace(-8.5e-3,8.5e-3,50));
X = Y_before_rot*sind(45);
Y = Y_before_rot*cosd(45);

IN = inpolyhedron(sample.F,sample.V,[X(:),Y(:),Z(:)],'FaceNormals',sample.N);

% determine mesh points inside kowari gauge volumes
y = repmat(kow_pts.p2_y'*1e3,5,1);
z = repmat(kow_pts.p2_z'*1e3,5,1);
corners_y = [0.5;0.5;-0.5;-0.5;0.5];
corners_z = [0.5;-0.5;-0.5;0.5;0.5];

y2 = [y+corners_y;nan(1,length(y))]; y2_save = y2*1e-3;
z2 = [z+corners_z;nan(1,length(y))]; z2_save = z2*1e-3;
 y2 = y2(:)*1e-3;
z2 = z2(:)*1e-3;
IN_kow = inpolygon(Y_before_rot,Z,y2,z2);

figure(1)
clf
p = patch('Faces',sample.F,'Vertices',sample.V);
p.FaceAlpha =0.1;
axis equal
hold on
h1 = plot3(X(IN),Y(IN),Z(IN),'bo');
h3 = plot3(X(IN_kow),Y(IN_kow),Z(IN_kow),'go');
hold off
xlabel('x')
ylabel('y')
zlabel('z')
view([76.1, 39.6])
title('Cutting plane 2')
legend([h1,h3],'Cutting plane 2','KOWARI region')


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

% load hypterparameters
load('optimisation_results.mat')
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
nu = 0.28;      % a good assumption

load('traction_points_20') % load traction point locations

Xcom = [X(IN);Xt(:)];
Ycom = [Y(IN);Yt(:)];
Zcom = [Z(IN);Zt(:)];

num_t = length(X(IN));
num_b = size(Xt,1); 

[Phi, SLambda,~,Phi_T] = beltrami_approx(num_basis,A,B,C,D,E,F,nu,entry,exit,Xcom,Ycom,Zcom,nsegs);

Phi_bound = Phi_T(num_t*6+1:num_t*6+num_b*6,:);       % extract boundary phi s
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

Var_f =  Phi_T*(linsolve(CZ,linsolve(CZ,Phi_T.',optsT),opts));
av_std_recon = mean(sqrt(diag(Var_f)));
disp(['Average reconstruction standard deviation = ' num2str(av_std_recon)])

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


%% load plane 2 KOWARI data
kow_eps = load('CP_kowstrains');
fxx_kow = scatteredInterpolant(kow_pts.p2_y,kow_pts.p2_z,kow_eps.p2_eps(1,:)','nearest','nearest');
fyy_kow = scatteredInterpolant(kow_pts.p2_y,kow_pts.p2_z,kow_eps.p2_eps(2,:)','nearest','nearest');
fzz_kow = scatteredInterpolant(kow_pts.p2_y,kow_pts.p2_z,kow_eps.p2_eps(3,:)','nearest','nearest');
fxy_kow = scatteredInterpolant(kow_pts.p2_y,kow_pts.p2_z,kow_eps.p2_eps(4,:)','nearest','nearest');
fxz_kow = scatteredInterpolant(kow_pts.p2_y,kow_pts.p2_z,kow_eps.p2_eps(5,:)','nearest','nearest');
fyz_kow = scatteredInterpolant(kow_pts.p2_y,kow_pts.p2_z,kow_eps.p2_eps(6,:)','nearest','nearest');

eps_xx_kow = NaN(size(X));
eps_yy_kow = NaN(size(X));
eps_zz_kow = NaN(size(X));
eps_xy_kow = NaN(size(X));
eps_xz_kow = NaN(size(X));
eps_yz_kow = NaN(size(X));
eps_xx_kow(IN_kow) = fxx_kow(Y_before_rot(IN_kow),Z(IN_kow));
eps_yy_kow(IN_kow) = fyy_kow(Y_before_rot(IN_kow),Z(IN_kow));
eps_zz_kow(IN_kow) = fzz_kow(Y_before_rot(IN_kow),Z(IN_kow));
eps_xy_kow(IN_kow) = fxy_kow(Y_before_rot(IN_kow),Z(IN_kow));
eps_xz_kow(IN_kow) = fxz_kow(Y_before_rot(IN_kow),Z(IN_kow));
eps_yz_kow(IN_kow) = fyz_kow(Y_before_rot(IN_kow),Z(IN_kow));

absval_xx = max(abs([eps_xx_kow(:);epsxx_pred(IN_kow)]));
absval_yy = max(abs([eps_yy_kow(:);epsyy_pred(IN_kow)]));
absval_zz = max(abs([eps_zz_kow(:);epszz_pred(IN_kow)]));
absval_xy = max(abs([eps_xy_kow(:);epsxy_pred(IN_kow)]));
absval_xz = max(abs([eps_xz_kow(:);epsxz_pred(IN_kow)]));
absval_yz = max(abs([eps_yz_kow(:);epsyz_pred(IN_kow)]));

xlims = [min(y2) max(y2)];
ylims = [min(z2) max(z2)];

%% load the fea results
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

FigHandle = figure(2);
clf

FigHandle.Position = [150   240   920   654];
colormap(cmap)

axis_width = 0.13;
axis_height = 0.15;
top = 0.8;
left = 0.1;
v_space = 0.01;
h_space = 0.01;
gap = 0.04;
fontsize = 20;
tmp = load('RdYlBu_up.mat');
cmap = tmp.RdYlBu;


subplot 461
pcolor(Y_before_rot,Z,eps_xx_kow)
shading flat
title('$\epsilon_{xx}$','Interpreter','latex','FontSize',fontsize*1.25)
axis equal
caxis([-absval_xx absval_xx])
xlim(xlims)
ylim(ylims)
axis off
ylabel('KOWARI')
h1 = gca;
h1.Position = [left  top    axis_width    axis_height];
t1 = text(+0.002,-0.0045,'KOWARI','Interpreter','latex','FontSize',fontsize);
set(t1,'Rotation',90,'HorizontalAlignment','center')

subplot 462
pcolor(Y_before_rot,Z,eps_yy_kow)
shading flat
title('$\epsilon_{yy}$','Interpreter','latex','FontSize',fontsize*1.25)
axis equal
caxis([-absval_yy absval_yy])
xlim(xlims)
ylim(ylims)
axis off
h1 = gca;
h1.Position = [left+axis_width+h_space  top    axis_width    axis_height];


subplot 463
pcolor(Y_before_rot,Z,eps_zz_kow)
shading flat
title('$\epsilon_{zz}$','Interpreter','latex','FontSize',26)
axis equal
caxis([-absval_zz absval_zz])
xlim(xlims)
ylim(ylims)
axis off
h1 = gca;
h1.Position = [left+2*(axis_width+h_space)  top    axis_width    axis_height];

subplot 464
pcolor(Y_before_rot,Z,eps_xy_kow)
shading flat
title('$\epsilon_{xy}$','Interpreter','latex','FontSize',fontsize*1.25)
axis equal
caxis([-absval_xy absval_xy])
xlim(xlims)
ylim(ylims)
axis off
h1 = gca;
h1.Position = [left+3*(axis_width+h_space)  top    axis_width    axis_height];

subplot 465
pcolor(Y_before_rot,Z,eps_xz_kow)
shading flat
title('$\epsilon_{xz}$','Interpreter','latex','FontSize',26)
axis equal
caxis([-absval_xz absval_xz])
xlim(xlims)
ylim(ylims)
axis off
h1 = gca;
h1.Position = [left+4*(axis_width+h_space)  top    axis_width    axis_height];

subplot 466
pcolor(Y_before_rot,Z,eps_yz_kow)
shading flat
title('$\epsilon_{yz}$','Interpreter','latex','FontSize',fontsize*1.25)
axis equal
caxis([-absval_yz absval_yz])
xlim(xlims)
ylim(ylims)
axis off
h1 = gca;
h1.Position = [left+5*(axis_width+h_space)  top    axis_width    axis_height];



% reconstruction results
subplot 467
tmp = flipud(fliplr(epsxx_pred));       % because the kowari and other arent lined up the same
tmp(~IN_kow) = NaN;
pcolor(Y_before_rot,Z,tmp)
shading interp
axis equal
caxis([-absval_xx absval_xx])
xlim(xlims)
ylim(ylims)
axis off
% colorbar('SouthOutside')
h1 = gca;
h1.Position = [left  top-axis_height-v_space    axis_width    axis_height];
t1 = text(+0.002,-0.005,'RADEN','Interpreter','latex','FontSize',fontsize);
set(t1,'Rotation',90,'HorizontalAlignment','center')

subplot 468
tmp = flipud(fliplr(epsyy_pred)); 
tmp(~IN_kow) = NaN;
pcolor(Y_before_rot,Z,tmp)
shading interp
axis equal
caxis([-absval_yy absval_yy])
xlim(xlims)
ylim(ylims)
axis off
% colorbar('SouthOutside')
h1 = gca;
h1.Position = [left+axis_width+h_space  top-axis_height-v_space    axis_width    axis_height];

subplot 469
tmp = flipud(fliplr(epszz_pred)); 
tmp(~IN_kow) = NaN;
pcolor(Y_before_rot,Z,tmp)
shading interp
axis equal
caxis([-absval_zz absval_zz])
xlim(xlims)
ylim(ylims)
axis off
% colorbar('SouthOutside')
h1 = gca;
h1.Position = [left+2*(axis_width+h_space)  top-axis_height-v_space    axis_width    axis_height];

subplot(4,6,10)
tmp = flipud(fliplr(epsxy_pred)); 
tmp(~IN_kow) = NaN;
pcolor(Y_before_rot,Z,tmp)
shading interp
axis equal
caxis([-absval_xy absval_xy])
xlim(xlims)
ylim(ylims)
axis off
% colorbar('SouthOutside')
h1 = gca;
h1.Position = [left+3*(axis_width+h_space)  top-axis_height-v_space    axis_width    axis_height];


subplot(4,6,11)
tmp = flipud(fliplr(epsxz_pred)); 
tmp(~IN_kow) = NaN;
pcolor(Y_before_rot,Z,tmp)
shading interp
axis equal
caxis([-absval_xz absval_xz])
xlim(xlims)
ylim(ylims)
axis off
% colorbar('SouthOutside')
h1 = gca;
h1.Position = [left+4*(axis_width+h_space)  top-axis_height-v_space    axis_width    axis_height];


subplot(4,6,12)
tmp = flipud(fliplr(epsyz_pred)); 
tmp(~IN_kow) = NaN;
pcolor(Y_before_rot,Z,tmp)
shading interp
axis equal
caxis([-absval_yz absval_yz])
xlim(xlims)
ylim(ylims)
axis off
% colorbar('SouthOutside')
h1 = gca;
h1.Position = [left+5*(axis_width+h_space)  top-axis_height-v_space    axis_width    axis_height];


% FEA
subplot(4,6,13)
tmp = eps_xx_fea; 
tmp(~IN_kow) = NaN;
pcolor(Y_before_rot,Z,tmp)
shading interp
axis equal
caxis([-absval_xx absval_xx])
xlim(xlims)
ylim(ylims)
axis off
colorbar('SouthOutside')
h1 = gca;
h1.Position = [left  top-2*(axis_height-v_space)-gap   axis_width    axis_height];
t1 = text(+0.002,-0.005,'FEA','Interpreter','latex','FontSize',fontsize);
set(t1,'Rotation',90,'HorizontalAlignment','center')
% 
subplot(4,6,14)
tmp = flipud(fliplr(eps_yy_fea)); 
tmp(~IN_kow) = NaN;
pcolor(Y_before_rot,Z,tmp)
shading interp
axis equal
caxis([-absval_yy absval_yy])
xlim(xlims)
ylim(ylims)
axis off
colorbar('SouthOutside')
h1 = gca;
h1.Position = [left+axis_width+h_space  top-2*(axis_height-v_space)-gap   axis_width    axis_height];
% 
subplot(4,6,15)
tmp = flipud(fliplr(eps_zz_fea)); 
tmp(~IN_kow) = NaN;
pcolor(Y_before_rot,Z,tmp)
shading interp
axis equal
caxis([-absval_zz absval_zz])
xlim(xlims)
ylim(ylims)
axis off
colorbar('SouthOutside')
h1 = gca;
h1.Position = [left+2*(axis_width+h_space)  top-2*(axis_height-v_space)-gap   axis_width    axis_height];
% 
subplot(4,6,16)
tmp = flipud(fliplr(eps_xy_fea)); 
tmp(~IN_kow) = NaN;
pcolor(Y_before_rot,Z,tmp)
shading interp
axis equal
caxis([-absval_xy absval_xy])
xlim(xlims)
ylim(ylims)
axis off
colorbar('SouthOutside')
h1 = gca;
h1.Position = [left+3*(axis_width+h_space)  top-2*(axis_height-v_space)-gap   axis_width    axis_height];

% 
subplot(4,6,17)
tmp = flipud(fliplr(eps_xz_fea)); 
tmp(~IN_kow) = NaN;
pcolor(Y_before_rot,Z,tmp)
shading interp
axis equal
caxis([-absval_xz absval_xz])
xlim(xlims)
ylim(ylims)
axis off
colorbar('SouthOutside')
h1 = gca;
h1.Position = [left+4*(axis_width+h_space)  top-2*(axis_height-v_space)-gap   axis_width    axis_height];

% 
subplot(4,6,18)
tmp = flipud(fliplr(eps_yz_fea)); 
tmp(~IN_kow) = NaN;
pcolor(Y_before_rot,Z,tmp)
shading interp
axis equal
caxis([-absval_yz absval_yz])
xlim(xlims)
ylim(ylims)
axis off
colorbar('SouthOutside')
h1 = gca;
h1.Position = [left+5*(axis_width+h_space)  top-2*(axis_height-v_space)-gap   axis_width    axis_height];


%%
[r c] = size(y2_save);
error = NaN(c,1);
for i = 1:c
   inds =  inpolygon(Y_before_rot,Z,y2_save(:,i),z2_save(:,i));
   error(i) = mean([epsxx_pred(inds);epsyy_pred(inds);epszz_pred(inds);epsxy_pred(inds);epsxz_pred(inds);epsyz_pred(inds)]-...
    [eps_xx_kow(inds);eps_yy_kow(inds);eps_zz_kow(inds);eps_xy_kow(inds);eps_xz_kow(inds);eps_yz_kow(inds)]);
    
end

fprintf('\n')
disp('Computed difference between KOWARI measurements and Reconstruction')
disp(['mean difference = ' num2str(mean(error))])
disp(['mean difference magnitude = ' num2str(mean(abs(error)))])
rel_error = mean(abs(error))/max(abs([eps_xx_kow(IN_kow);eps_yy_kow(IN_kow);eps_zz_kow(IN);eps_xy_kow(IN_kow);eps_xz_kow(IN_kow);eps_yz_kow(IN_kow)]));
disp(['mean relative difference = ' num2str(rel_error)])

%% plot comparison over full plane with FEA
absval_xx = max(abs([epsxx_pred(IN)]));
absval_yy = max(abs([epsyy_pred(IN)]));
absval_zz = max(abs([epszz_pred(IN)]));
absval_xy = max(abs([epsxy_pred(IN)]));
absval_xz = max(abs([epsxz_pred(IN)]));
absval_yz = max(abs([epsyz_pred(IN)]));

FigHandle = figure(3);
clf

FigHandle.Position = [150   240   920   654];
colormap(cmap)

axis_width = 0.25;
axis_height = 0.25;
top = 0.65;
left = 0.1;
v_space = 0.035;
h_space = 0.025;
% gap = 0.05
fontsize = 25;
tmp = load('RdYlBu_up.mat');
cmap = tmp.RdYlBu;



% reconstruction results
subplot 331
tmp = flipud(fliplr(epsxx_pred));       % because the kowari and other arent lined up the same
pcolor(Y_before_rot,Z,tmp)
title('$\epsilon_{xx}$','Interpreter','latex','FontSize',fontsize*1.5)
shading interp
axis equal
caxis([-absval_xx absval_xx])
axis off
h1 = gca;
% h1.Position
h1.Position = [left  top    axis_width    axis_height];
t1 = text(-0.015,-0.00,'RADEN','Interpreter','latex','FontSize',fontsize);
set(t1,'Rotation',90,'HorizontalAlignment','center')


subplot 332
tmp = flipud(fliplr(epsxy_pred)); 
pcolor(Y_before_rot,Z,tmp)
shading interp
title('$\epsilon_{xy}$','Interpreter','latex','FontSize',fontsize*1.5)
axis equal
caxis([-absval_xy absval_xy])
axis off
h1 = gca;
h1.Position = [left+(axis_width+h_space)  top    axis_width    axis_height];

subplot(3,3,3)
tmp = flipud(fliplr(epsyz_pred)); 
pcolor(Y_before_rot,Z,tmp)
title('$\epsilon_{yz}$','Interpreter','latex','FontSize',fontsize*1.5)
shading interp
axis equal
caxis([-absval_yz absval_yz])
axis off
h1 = gca;
h1.Position = [left+2*(axis_width+h_space)  top    axis_width    axis_height];

% FEA
subplot(3,3,4)
tmp = flipud(fliplr(eps_xx_fea)); 
pcolor(Y_before_rot,Z,tmp)
shading interp
axis equal
caxis([-absval_xx absval_xx])
axis off
colorbar('SouthOutside')
h1 = gca;
h1.Position = [left  top-1*(axis_height+v_space)   axis_width    axis_height];
t1 = text(-0.015,-0.00,'FEA','Interpreter','latex','FontSize',fontsize);
set(t1,'Rotation',90,'HorizontalAlignment','center')


subplot(3,3,5)
tmp = flipud(fliplr(eps_xy_fea)); 
pcolor(Y_before_rot,Z,tmp)
shading interp
axis equal
caxis([-absval_xy absval_xy])
axis off
colorbar('SouthOutside')
h1 = gca;
h1.Position = [left+1*(axis_width+h_space)  top-1*(axis_height+v_space)   axis_width    axis_height];

subplot(3,3,6)
tmp = flipud(fliplr(eps_yz_fea)); 
pcolor(Y_before_rot,Z,tmp)
shading interp
axis equal
caxis([-absval_yz absval_yz])
axis off
colorbar('SouthOutside')
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
