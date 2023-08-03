function [U] = rdg_allpairs_admm(Mm, alpha_hat)

% ADMM algorithm for computing regularized all-pairs geodesic distances
%   regularized with the Dirichlet regularizer
%
% Inputs:
%   Mm - MeshClass
%   alpha_hat - regularizer weight
%
% Outputs:
%   U - the computed all-pairs regularized distance



%% mesh data
nv = Mm.nv;
nf = Mm.nf;
va = Mm.va;
ta = Mm.ta;
G = Mm.G;
Ww = Mm.Ww;
vasq = sqrt(va);
vaMat = spdiags(va,0,nv,nv);
vainv = 1./va;
vavatMat = va*va';
tasq = repmat(sqrt(ta),3,1); 
taMat = spdiags(repmat(ta,3,1),0,3*nf ,3*nf);
div = G'*taMat;
    
alpha = alpha_hat*sqrt(sum(va));


%% ADMM parameters
rho     = 2*sqrt(sum(va));
rho2    = 10/sqrt(sum(va));
niter   = 20000;
QUIET   = 0;
ABSTOL  = 1e-6;
RELTOL  = 2e-4;
mu      = 10;  % >1
tauinc  = 2;   % >1
taudec  = 2;   % >1
alphak  = 1.7; % over-relaxation

thresh1 = sqrt(3*nf)*ABSTOL*(sum(va));
thresh2 = sqrt(nv)*ABSTOL*(sum(va))^2;
thresh3 = sqrt(nv)*ABSTOL*sqrt(sum(va))^3;
thresh4 = sqrt(nv)*ABSTOL*(sum(va));

varRho = 1; % determine whether to use a varying penalty parameter
rho1or2changed = 1;


% initialization:
X = zeros(nv,nv);       % all-pairs distances, represent the gradient along the columns
R = zeros(nv,nv);       % all-pairs distances, represent the gradient along the rows
U = zeros(nv,nv);       % dual consensus variable, all-pairs distance matrix
Z = zeros(3*nf,nv);     % auxiliary variable, GX = Z
Q = zeros(3*nf,nv);     % auxiliary variable, GR = Q
Y = zeros(3*nf,nv);     % dual variable
S = zeros(3*nf,nv);     % dual variable
H = zeros(nv,nv);       % dual consensus variable
K = zeros(nv,nv);       % dual consensus variable
div_Z = sparse(nv,nv);
div_Q = sparse(nv,nv);
div_Y = sparse(nv,nv);
div_S = sparse(nv,nv);


history.r_norm = zeros(niter,1);
history.s_norm = zeros(niter,1);
history.eps_pri = zeros(niter,1);
history.eps_dual = zeros(niter,1);
history.r_norm2 = zeros(niter,1);
history.s_norm2 = zeros(niter,1);
history.eps_pri2 = zeros(niter,1);
history.eps_dual2 = zeros(niter,1);
history.r_xr1 = zeros(niter,1);
history.r_xr2 = zeros(niter,1);
history.s_xr = zeros(niter,1);
history.eps_pri_xr = zeros(niter,1);
history.eps_dual_xr = zeros(niter,1);

if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n', ...
        'iter', 'r norm', 'eps pri', 's norm', 'eps dual','r2 norm', 'eps2 pri', 's2 norm', 'eps2 dual', 'xr1 r norm', 'xr1 eps pri', 'xr2 r norm', 'xr2 eps pri', 'xr s norm', 'xr eps dual');
end

% Pre-factorization
if ~varRho
    A2in_p = (alpha+rho)*Ww+rho2*vaMat;
    [L,~,P] = chol(A2in_p,'lower');
end


for ii = 1:niter

    % step 1 - X,R-minimization
    bx = (0.5*vavatMat.*vainv' - div_Y + rho*div_Z - va.*H + rho2*va.*U );
    br = (0.5*vavatMat.*vainv' - div_S + rho*div_Q - va.*K + rho2*va.*U' );

    if varRho && rho1or2changed
        A2in_p = (alpha+rho)*Ww+rho2*vaMat;
        [L,g,P] = chol(A2in_p,'lower');
    end

    X = P * (L'  \ (L \ (P' * bx )));
    R = P * (L'  \ (L \ (P' * br )));

    Gx = G*X;
    Gr = G*R;
    
    % step 2 - Z,Q,U-minimization:
    Zold = Z;
    div_Zold = div_Z;
    Z = (1/rho)*Y + Gx;
    Z = reshape(Z,nf,3,nv); 
    Z = Z./max(1,sqrt(sum(Z.^2,2)));
    Z = reshape(Z,3*nf,nv);
    div_Z = div*Z;
    
    Qold = Q;
    div_Qold = div_Q;
    Q = (1/rho)*S + Gr;
    Q = reshape(Q,nf,3,nv); 
    Q = Q./max(1,sqrt(sum(Q.^2,2)));
    Q = reshape(Q,3*nf,nv); 
    div_Q = div*Q;

    Uold = U;
    U1 = 0.5*((1/rho2)*(H+K') + X+R');
    U = U1-diag(diag(U1));
    U(U<0) = 0;

    % step 3 - dual variable update
    Y = Y + rho*(alphak*Gx+(1-alphak)*Zold-Z);
    S = S + rho*(alphak*Gr+(1-alphak)*Qold-Q);
    H = H + rho2*(alphak*X+(1-alphak)*Uold-U);
    K = K + rho2*(alphak*R+(1-alphak)*Uold'-U');
    div_Y = div*Y;
    div_S = div*S;

    % residuals update
    GxW = tasq.*Gx.*vasq';
    ZW =  tasq.*Z.*vasq';
	history.r_norm(ii)  = norm(GxW-ZW,'fro');
    history.eps_pri(ii)  = thresh1 + RELTOL*max(norm(GxW,'fro'), norm(ZW,'fro'));
    history.s_norm(ii)  = rho*norm((div_Z - div_Zold).*va','fro');
    history.eps_dual(ii) = thresh2 + RELTOL*norm(div_Y.*va','fro');
    
    GrW = tasq.*Gr.*vasq';
    QW =  tasq.*Q.*vasq';
    history.r_norm2(ii)  = norm(GrW-QW,'fro');
    history.eps_pri2(ii)  = thresh1 + RELTOL*max(norm(GrW,'fro'), norm(QW,'fro'));
    history.s_norm2(ii)  = rho*norm((div_Q - div_Qold).*va','fro');
    history.eps_dual2(ii) = thresh2 + RELTOL*norm(div_S.*va','fro');
    
    history.r_xr1(ii) = norm(vasq.*(X-U).*vasq','fro');
    history.r_xr2(ii) = norm(vasq.*(R-U').*vasq','fro');
    history.eps_pri_xr1(ii) = thresh3 + RELTOL*min(norm(vasq.*X.*vasq','fro'), norm(vasq.*U.*vasq','fro'));
    history.eps_pri_xr2(ii) = thresh3 + RELTOL*min(norm(vasq.*R.*vasq','fro'), norm(vasq.*U'.*vasq','fro'));
    history.s_xr(ii) = sqrt(2)*rho2*norm(vasq.*(U-Uold).*vasq','fro');
    history.eps_dual_xr(ii) = thresh4 + RELTOL*0.5*(norm(vasq.*H.*vasq','fro')+norm(vasq.*K.*vasq','fro'));


    if ~QUIET && ~mod(ii,10)
        fprintf('%3d |%10.4f  %10.4f |%10.4f  %10.4f |%10.4f  %10.4f |%10.4f  %10.4f |%10.4f  %10.4f |%10.4f  %10.4f |%10.4f  %10.4f\n', ii, ...
            history.r_norm(ii), history.eps_pri(ii), ...
            history.s_norm(ii), history.eps_dual(ii), ...
            history.r_norm2(ii), history.eps_pri2(ii), ...
            history.s_norm2(ii), history.eps_dual2(ii), ...
            history.r_xr1(ii), history.eps_pri_xr1(ii), ...
            history.r_xr2(ii), history.eps_pri_xr2(ii), ...
            history.s_xr(ii), history.eps_dual_xr(ii));
    end

     if (history.r_norm(ii) < history.eps_pri(ii) && ...
       history.s_norm(ii) < history.eps_dual(ii) && ...
       history.r_norm2(ii) < history.eps_pri2(ii) && ...
       history.s_norm2(ii) < history.eps_dual2(ii) && ...
       history.r_xr1(ii) < history.eps_pri_xr1(ii) && ...
       history.r_xr2(ii) < history.eps_pri_xr2(ii) && ...
       history.s_xr(ii) < history.eps_dual_xr(ii))
         break;
    end

    % varying penalty parameter
    if varRho
        if (history.r_norm(ii)/history.eps_pri(ii) > mu*history.s_norm(ii)/history.eps_dual(ii)) && ...
            (history.r_norm2(ii)/history.eps_pri2(ii) > mu*history.s_norm2(ii)/history.eps_dual2(ii))
            rho = tauinc*rho;
            rho1or2changed = 1;
        elseif (history.s_norm(ii)/history.eps_dual(ii) > mu*history.r_norm(ii)/history.eps_pri(ii)) && ...
                (history.s_norm2(ii)/history.eps_dual2(ii) > mu*history.r_norm2(ii)/history.eps_pri2(ii))
            rho =  rho/taudec;
            rho1or2changed = 1;
        end
        if (history.r_xr1(ii)/history.eps_pri_xr1(ii) > mu*history.s_xr(ii)/history.eps_dual_xr(ii)) && ...
                (history.r_xr2(ii)/history.eps_pri_xr2(ii) > mu*history.s_xr(ii)/history.eps_dual_xr(ii))
            rho2 = tauinc*rho2;
            rho1or2changed = 1;
        elseif (history.s_xr(ii)/history.eps_dual_xr(ii) > mu*history.r_xr1(ii)/history.eps_pri_xr1(ii)) && ...
                (history.s_xr(ii)/history.eps_dual_xr(ii) > mu*history.r_xr2(ii)/history.eps_pri_xr2(ii))
            rho2 =  rho2/taudec;
            rho1or2changed = 1;
        end

    end
    
end

end





