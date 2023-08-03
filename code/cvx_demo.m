
clearvars; close all; clc;

addpath(genpath('code/')); 
addpath('data/');

% add cvx to path
cvx_folder = '/path/to/cvx'; % update to cvx_folder
addpath(genpath(cvx_folder));

filename = 'bent_pipe_closed_lr';
x0 = 4883; % source point / set


% load mesh
Mm = MeshClass(filename);
Mm.center_mesh;

alpha_hat = 0.03; % regularizer weight

alpha = alpha_hat*sqrt(sum(Mm.va));


%% Regularized - Dirichlet Energy
cvx_begin
    variable u_D(Mm.nv,1)
    grad = reshape(Mm.G*u_D,Mm.nf,3);
    
    maximize sum( Mm.va .* u_D ) - alpha*quad_form(u_D, Mm.Ww)
    subject to
        u_D(x0) == 0
        norms(grad,2,2) <= 1
cvx_end



%% Regularized - Linf
cvx_begin
    variable u_Linf(Mm.nv,1)
    grad = reshape(Mm.G*u_Linf,Mm.nf,3);
    
    maximize sum( Mm.va .* u_Linf ) - alpha*sum(Mm.ta .* pow_pos(norms(grad,inf,2),2))
    subject to
        u_Linf(x0) == 0
        norms(grad,2,2) <= 1
cvx_end


% rotate mesh
A = rotz(45);
VV = Mm.vertices*A;
Mm_rot = MeshClass(VV, Mm.faces);


cvx_begin
    variable u_Linf_rot(Mm.nv,1)
    grad = reshape(Mm_rot.G*u_Linf_rot,Mm_rot.nf,3);
    
    maximize sum( Mm_rot.va .* u_Linf_rot ) - alpha*sum(Mm_rot.ta .* pow_pos(norms(grad,inf,2),2))
    subject to
        u_Linf_rot(x0) == 0
        norms(grad,2,2) <= 1
cvx_end



%% Figures
u_all = [u_D(:); u_Linf(:); u_Linf_rot(:)];
umin = min(u_all);
umax = max(u_all);
nlines = 15;

Mm.visualizeDistances(u_D, x0, nlines, [umin, umax]); axis on;
Mm.visualizeDistances(u_Linf, x0, nlines, [umin, umax]); axis on;
Mm_rot.visualizeDistances(u_Linf_rot, x0, nlines, [umin, umax]); axis on;

