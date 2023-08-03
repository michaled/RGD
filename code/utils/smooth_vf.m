function w = smooth_vf(Mm, vf, n)
% n - power vector, 2 for line fields

% Based on code written by
% "PH-CPF: Planar Hexagonal Meshing using Coordinate Power Fields" by [Pluta et al., 2021]

nf = Mm.nf;
we = [];
    
locs = find(MeshClass.normv(vf) > 1e-5);            
nl = size(locs,1);

Aeq = sparse([1:2*nl], [locs, nf+locs], ones(2*nl,1),2*nl, 2*nf);

% -> local basis -> power n -> locs
beq = reshape(ff(Mm.EB*vf(:),n),[],2);
beq = beq(locs,:); 
beq = beq(:);

C = Mm.godf(n);
vf = zeros(nf*2,1);

% If there are no constraints, use the eigenvector
if size(Aeq,1) > 0
    x = lsqlin(C, vf, [], [], Aeq, beq);
else
    [x,~] = eigs(C, 1, 'SM');
end

% -> sqrt n -> extrinsic
w = reshape(Mm.EBI*ff(x,1/n),[],3);
% Normalize
w = MeshClass.normalize_vf(w);


end



% in: 2n x 1
% out: 2n x 1, f(x,n) = x.^n
% grad: 2n x 2n
function e = ff(in,n)
    s = size(in,1)/2;
    a = in(1:s); b = in(s+1:end);
    c = a+1i*b;
    cn = c.^n;
    e = [real(cn); imag(cn)];
end
