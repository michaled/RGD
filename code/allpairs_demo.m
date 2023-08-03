

clearvars; close all; clc;

addpath(genpath('code/')); 
addpath('data/');

filename = 'cat_rr';

% load mesh
Mm = MeshClass(filename);

alpha_hat = 0.03;
U = rdg_allpairs_admm(Mm, alpha_hat);


%% Figures
umin = min(U(:));
umax = max(U(:));
nlines = 15;

x0 = 494;
Mm.visualizeDistances(U(x0,:), x0, nlines, [umin, umax]);

x0 = 242;
Mm.visualizeDistances(U(x0,:), x0, nlines, [umin, umax]);
