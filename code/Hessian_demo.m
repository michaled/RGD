

clearvars; close all; clc;

addpath(genpath('code/')); 
addpath('data/');

filename = 'spot_rr';
filename1 = 'spot_rr_holes';
x0 = 2273; % source point / set
x01 = 7331; % source point / set


% load mesh
Mm = MeshClass(filename);
Mm1 = MeshClass(filename1);


%% Regularized - Hessian Energy
alpha_hat0 = 0.001; % scale invariant, regularizer weight

u_H = rdg_ADMM(Mm, x0, 'reg', 'H', 'alpha_hat', alpha_hat0);
u_H_holes = rdg_ADMM(Mm1, x01, 'reg', 'H', 'alpha_hat', alpha_hat0);
u_D_holes = rdg_ADMM(Mm1, x01, 'reg', 'D', 'alpha_hat', 0.1);


%% Figures
cam = load('spot_rr_cam.mat'); cam = cam.cam;
u_all = [u_H(:); u_H_holes(:); u_D_holes(:)];
umin = min(u_all);
umax = max(u_all);
nlines = 15;

Mm.visualizeDistances(u_H, x0, nlines, [umin, umax], cam);
Mm1.visualizeDistances(u_H_holes, x01, nlines, [umin, umax], cam);
Mm1.visualizeDistances(u_D_holes, x01, nlines, [umin, umax], cam);


