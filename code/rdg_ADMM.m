function u = rdg_ADMM(Mm, x0, varargin)

    % ADMM algorithm for computing regularized geodesic distances
    %
    % Inputs:
    %   Mm - MeshClass
    %   x0 - source set - vertex indices
    %   Optional:
    %       reg - which regularizer to use.
    %           'D' - dirichlet
    %           'vfa' - vector field alignment
    %           'H' - dirichlet
    %       alpha_hat - regularizer weight, scale invariant
    %       beta_hat - vector field alignment weight, scale-invariant, relevant for reg = 'vfa'
    %       vf - |F|x3 - vector field to align to, relevant for reg = 'vfa'
    %
    % Outputs:
    %   u - the computed regularized distance


p = inputParser;
p.KeepUnmatched=true;
addOptional(p,'reg','D');
addOptional(p,'alpha_hat',0.1);
addOptional(p,'beta_hat',0);
addOptional(p,'vf',0);
parse(p,varargin{:});


reg = p.Results.reg;
alpha_hat = p.Results.alpha_hat;
beta_hat = p.Results.beta_hat;
vf = p.Results.vf;


%% mesh data
vertices = Mm.vertices;
faces = Mm.faces;
nv = Mm.nv;
nf = Mm.nf;
va = Mm.va;
ta = Mm.ta;
G = Mm.G;
Ww = Mm.Ww;
tasq = repmat(sqrt(ta),3,1); 


% set parameters according to the desired regularizer
if strcmp(reg,'D')
    alpha = alpha_hat*sqrt(sum(va));
    varRho = 1; % determine whether to use a varying penalty parameter

elseif strcmp(reg,'H')
    alpha = alpha_hat*sqrt(sum(va)^3);

    % This code uses the Hessian computed by 
    % "A Smoothness Energy without Boundary Distortion for Curved Surfaces" by [Stein et al., 2020]
    % Ensure that the mex `curved_hessian` is available.
    % If not available, please follow the instructions in
    % https://github.com/odedstein/ASmoothnessEnergyWithoutBoundaryDistortionForCurvedSurfaces
    Ww_s = curved_hessian(vertices,faces);
    varRho = 0; % determine whether to use a varying penalty parameter
elseif strcmp(reg,'vfa')
    alpha = alpha_hat*sqrt(sum(va));
    beta = beta_hat*sqrt(sum(va));
    if max(Mm.normv(vf))<1e-10
        error('Vector field for alignment is empty');
    end
    Vmat = [[spdiags(vf(:,1).*vf(:,1),0,nf,nf), spdiags(vf(:,1).*vf(:,2),0,nf,nf), spdiags(vf(:,1).*vf(:,3),0,nf,nf)]; ...
        [spdiags(vf(:,2).*vf(:,1),0,nf,nf), spdiags(vf(:,2).*vf(:,2),0,nf,nf), spdiags(vf(:,2).*vf(:,3),0,nf,nf)]; ...
        [spdiags(vf(:,3).*vf(:,1),0,nf,nf), spdiags(vf(:,3).*vf(:,2),0,nf,nf), spdiags(vf(:,3).*vf(:,3),0,nf,nf)]];
    Ww_s = G'*spdiags(repmat((ta),3,1),0,3*nf ,3*nf)*(speye(3*nf,3*nf)+beta*Vmat)*G;
    varRho = 0; % determine whether to use a varying penalty parameter
else
    error('Unrecognized regularizer');
end



%% ADMM parameters
rho     = 2*sqrt(sum(va));
niter   = 10000;
QUIET   = 1;
ABSTOL  = 1e-5/2;
RELTOL  = 1e-2;
mu      = 10;  % >1
tauinc  = 2;   % >1
taudec  = 2;   % >1
alphak  = 1.7; % over-relaxation

if strcmp(reg,'H')
    ABSTOL  = ABSTOL/20;
    RELTOL  = RELTOL/20;
end

thresh1 = sqrt(3*nf)*ABSTOL*sqrt(sum(va));
thresh2 = sqrt(nv)*ABSTOL*(sum(va));


u_p = zeros(nv-length(x0),1);       % distance to x in M \ x0
y = zeros(3*nf,1);                  % dual variable
z = zeros(3*nf,1);                  % auxiliary variable
div_y = zeros(nv-length(x0),1);
div_z = zeros(nv-length(x0),1);

history.r_norm = zeros(niter,1);
history.s_norm = zeros(niter,1);
history.eps_pri = zeros(niter,1);
history.eps_dual = zeros(niter,1);


%% 
% Eliminating x0 (b.c):
nv_p = 1:nv;
nv_p(x0) = [];
va_p = va; va_p(x0) = [];
Ww_p = Ww; Ww_p(:,x0) = []; Ww_p(x0,:) = [];
G_p = G; G_p(:,x0) = []; G_pt = G_p';
div_p = G_pt.*repmat(ta,3,1)';

if (strcmp(reg,'vfa') || strcmp(reg,'H'))
    Ww_s_p = Ww_s;
    Ww_s_p(:,x0) = [];
    Ww_s_p(x0,:) = [];
end


%%
if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t\n', ...
        'iter', 'r norm', 'eps pri', 's norm', 'eps dual');
end


% Pre-factorization
if strcmp(reg,'D')
    [L,~,P] = chol(Ww_p,'lower');
else % 'H', 'vfa'
    if ~varRho
        [L,~,P] = chol(alpha*Ww_s_p+rho*Ww_p,'lower');
    end
end


for ii = 1:niter
   
    % step 1 - u-minimization
    b = va_p-div_y+rho*div_z;

    if strcmp(reg,'D')
        u_p = P * (L'  \ (L \ (P' * b ))) /(alpha+rho);
    else % 'H', 'vfa'
        if ~varRho
            u_p = P * (L'  \ (L \ (P' * b )));
        else % strcmp(reg,'H') && varRho
            u_p = (alpha*Ww_s_p+rho*Ww_p) \ b;
        end
    end
    Gx = G_p*u_p;

    % step 2 - z-minimization
    zold = z;
    div_zold = div_z;
    z = (1/rho)*y + Gx;
    z = reshape(z,nf,3)'; 
    norm_z = sqrt(sum(z.^2));
    norm_z(norm_z<1) = 1;
    z = z./norm_z;
    z = z'; z = z(:);
    div_z = div_p*z;
    
    % step 3 - dual variable update
    y = y + rho*(alphak*Gx+(1-alphak)*zold-z);
    div_y = div_p*y;

    % residuals update
    tasqGx = tasq.*Gx;
    tasqZ = tasq.*z;
    history.r_norm(ii)  = norm(tasqGx-tasqZ,'fro');
    history.s_norm(ii)  = rho*norm(div_z-div_zold,'fro');
    history.eps_pri(ii)  = thresh1 + RELTOL*max(norm(tasqGx,'fro'), norm(tasqZ,'fro'));
    history.eps_dual(ii) = thresh2 + RELTOL*norm(div_y,'fro');

    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n', ii, ...
            history.r_norm(ii), history.eps_pri(ii), ...
            history.s_norm(ii), history.eps_dual(ii));
    end

    % stopping criteria
    if ii>1 && (history.r_norm(ii) < history.eps_pri(ii) && ...
       history.s_norm(ii) < history.eps_dual(ii))
        break;
    end

    % varying penalty parameter
    if varRho
        if history.r_norm(ii) > mu*history.s_norm(ii)
            rho = tauinc*rho;
        elseif history.s_norm(ii) > mu*history.r_norm(ii)
            rho =  rho/taudec;
        end
    end
    
end


u = zeros(nv,1);
u(nv_p) = u_p;


end

