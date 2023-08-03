classdef MeshClass < handle
    
    properties (Access='public')
        name            % input mesh name 

        vertices
        faces
        edges


        nv              % #vertices
        nf              % #faces
        ne              % #edges
        nie             % #interior edges
        
        va              % vertex areas
        ta              % face areas
        ea              % edge areas


        Nf              % normal-per-face
        Nv              % normal-per-vertex

        E1, E2, E3      % triangle edges
        E               % unique triangle edges
        R               % in-plane rotation operator wrt normal

        F1, F2          % frame per face
        EB, EBI         % an edge based basis for vfs.

        e2t, t2e, v2e   % map edges to triangles, triangles to edges  
                        % vertices to edges
        e2t1, e2t2      % sparse matrices mapping edges to left and 
                        % right triangles on inner edges.


        ie              % 1 if interior edge
        inner_edges     % indices of interior edges
        bv              % boundary vertices
        bf              % faces which have boundary vertices
        
        Lap             % cotangent laplacian
        Ww              % cotangent laplacian weights
        Aa              % diagonal matrix containing va
        
        G               % gradient on vertex functions
        D               % divergence on face fields

    end
    
    properties (Access = public) % FIX
    end
    
    methods
        % parameters options (bases sizes are optional):
        % #1: vertices, faces, LB_basis_size, LB_type, name
        % #2: meshname, LB_basis_size, LB_type
        function obj = MeshClass(varargin)
            
            if size(varargin{1},1) >= 3 && size(varargin{1},2) == 3
                X = varargin{1};
                T = varargin{2};
				
				if nargin > 2
                    obj.name = varargin{3};
                end
            else
				try
					[X, T] = readOff([varargin{1} '.off']);
				catch
					error('Problem reading file');
				end
                T = double(T);
                if isempty(find(varargin{1}=='/', 1, 'last'))
                    obj.name = varargin{1};
                else
                    obj.name = varargin{1}(find(varargin{1}=='/', 1, 'last')+1:end);
                end
            end
            obj.vertices = X;
            obj.faces = T;
            obj.compute_all;
        end
        
        function compute_all(obj)
            obj.nv = size(obj.vertices,1);
            obj.nf = size(obj.faces,1);
            Nf = cross(obj.vertices(obj.faces(:,1),:)-obj.vertices(obj.faces(:,2),:), ...
                obj.vertices(obj.faces(:,1),:) - obj.vertices(obj.faces(:,3),:));
            obj.Nf = MeshClass.normalize_vf(Nf);
            obj.ta = sqrt(sum(Nf.^2,2))/2;
            obj.Nv = vertex_normals( obj );
            obj.va = obj.calculatefvConnectivity()' * obj.ta / 3;
			obj.E1 = obj.vertices(obj.faces(:,2),:) - obj.vertices(obj.faces(:,3),:);
            obj.E2 = obj.vertices(obj.faces(:,1),:) - obj.vertices(obj.faces(:,3),:);
            obj.E3 = obj.vertices(obj.faces(:,1),:) - obj.vertices(obj.faces(:,2),:);
            obj.R = rot( obj );
            EE = sort([obj.faces(:,1) obj.faces(:,2) ; obj.faces(:,2) obj.faces(:,3) ; obj.faces(:,3) obj.faces(:,1) ],2);
            obj.E = unique(EE,'rows');
            [obj.edges,obj.e2t,obj.t2e, ...
                obj.e2t1, obj.e2t2, ...
                obj.v2e,obj.ie, ...
                obj.ne, obj.nie, obj.inner_edges, obj.bv, obj.bf] = ...
                nc_data( obj );
            obj.ea = edge_areas( obj );
            obj.edge_basis;
            obj.compute_LB;
            obj.GG;
            obj.DD;
        end


        function [ fv, I_F2V ] = interpulateFace2Vertices( mesh, fF)

            Af = mesh.ta;
            Av = mesh.va;

            I_F2V = sparse([mesh.faces(:,1); mesh.faces(:,2); mesh.faces(:,3)], ...
                [(1:mesh.nf)' ; (1:mesh.nf)' ; (1:mesh.nf)' ], ...
                [(1/3)*Af./Av(mesh.faces(:,1)) ; (1/3)*Af./Av(mesh.faces(:,2)); (1/3)*Af./Av(mesh.faces(:,3))]);
            
            if nargin>1
                fv = I_F2V*fF;
            else
                fv = [];
            end

        end
        

        function [ fF, I_V2F ] = interpulateVertices2Face( mesh, fV )

            Afinv = sparse( 1:mesh.nf, 1:mesh.nf, 1./mesh.ta);
            Av = sparse( 1:mesh.nv, 1:mesh.nv, mesh.va);
            [ ~, I_F2V ] = mesh.interpulateFace2Vertices( ones(mesh.nf,1) );
            I_V2F = Afinv*I_F2V'*Av;
            if nargin>1
                fF = I_V2F*fV;
            else
                fF = [];
            end
        end


        function [ NE1, NE2, EB, EBI ] = edge_basis( mesh )
            
            NE1 = MeshClass.normalize_vf( mesh.E1 );
            NE2 = reshape( mesh.R * NE1(:), [], 3 );

            I = repmat(1:mesh.nf,1,3); I = I(:); 
            J = (1:mesh.nf)'; J = [J; J+mesh.nf; J+2*mesh.nf];
            
            B1 = sparse(I,J,NE1,mesh.nf,3*mesh.nf);
            B2 = sparse(I,J,NE2,mesh.nf,3*mesh.nf);
            
            EB = [B1; B2];
            EBI = EB';

            mesh.F1 = NE1;
            mesh.F2 = NE2;
            mesh.EB = EB;
            mesh.EBI = EBI;
        end


        function [ R ] = rot( mesh )
            sf = mesh.nf; n = mesh.Nf;
            
            II = repmat((1:sf)',1,2); II = [II II+sf II+2*sf]'; II = II(:);
            JJ1 = (1:sf)'; JJ2 = JJ1+sf; JJ3 = JJ2+sf;
            JJ = [JJ2 JJ3 JJ1 JJ3 JJ1 JJ2]'; JJ = JJ(:);
            SS = [-n(:,3) n(:,2) n(:,3) -n(:,1) -n(:,2) n(:,1)]'; SS = SS(:);
            
            R = sparse(II,JJ,SS,3*sf,3*sf);
        end


        function fvConnectivity = calculatefvConnectivity(obj)
            fvConnectivity = sparse(repmat(1:obj.nf, 1, 3), ...
                obj.faces, ...
                ones(size(obj.faces)), ...
                obj.nf, obj.nv);
        end


        function [ baryCenters ] = baryCentersCalc( mesh )
            v1 = mesh.vertices(mesh.faces(:,1),:);
            v2 = mesh.vertices(mesh.faces(:,2),:);
            v3 = mesh.vertices(mesh.faces(:,3),:);
            baryCenters = (1/3)*(v1+v2+v3);
        end


        function [ ea ] = edge_areas( mesh )
            T = double( mesh.faces );
            I = [T(:,2);T(:,3);T(:,1)];
            J = [T(:,3);T(:,1);T(:,2)];
            S = repmat(mesh.ta/3,3,1);
            In = [I;J];
            Jn = [J;I];
            Sn = [S;S];
            W = sparse(In,Jn,Sn,mesh.nv,mesh.nv);
            ea = zeros(length(mesh.edges),1);
            for i=1:length(ea)
                ea(i) = W(mesh.edges(i,1),mesh.edges(i,2));
                ea(i) = ea(i) + W(mesh.edges(i,2),mesh.edges(i,1));
            end
        end


        function rvf = rotate_vf( mesh, vf )
            vf = reshape(vf,mesh.nf,3);
            rvf = cross( mesh.Nf, vf );
        end


        function compute_LB(mesh)
            mesh.Ww = cotLaplacian(mesh);
            laplacian = sparse(1:mesh.nv, 1:mesh.nv, 1./mesh.va) * mesh.Ww;
            mesh.Lap = laplacian;
            mesh.Aa = spdiags(mesh.va,0,mesh.nv,mesh.nv);
        end
        

        function grad_op = GG(mesh)            
            I = repmat(1:mesh.nf,3,1);
            II = [I(:); I(:)+mesh.nf; I(:)+2*mesh.nf];
            J = double( mesh.faces' );
            JJ = [J(:); J(:); J(:)];
            RE1 = mesh.rotate_vf( mesh.E1 );
            RE2 = mesh.rotate_vf( mesh.E2 );
            RE3 = mesh.rotate_vf( mesh.E3 );
            TA = mesh.ta;
            S = [-RE1(:) RE2(:) -RE3(:)]'; SS = S(:);
            
            G = sparse(II,JJ,SS,3*mesh.nf,mesh.nv);
            ITA = spdiags(.5*repmat(1./TA,3,1),0,3*mesh.nf,3*mesh.nf);
            
            grad_op = ITA*G;
            if any(isnan(grad_op(:)))
                warning('Grad: NANs exist');
            end
            grad_op(isnan(grad_op)) = 0;
            mesh.G = grad_op;
        end
              

        function [ D ] = DD( mesh )
            IVA = spdiags(1./mesh.va,0,mesh.nv,mesh.nv);
            TA = spdiags([mesh.ta; mesh.ta; mesh.ta],0,3*mesh.nf,3*mesh.nf);
            D = - IVA * mesh.G' * TA;
            mesh.D = D;
        end

        function [ Nv ] = vertex_normals( mesh )        
            
            I = double( repmat(mesh.faces(:),3,1) );
            J = repmat(1:3,3*mesh.nf,1); J = J(:);
            
            TA = spdiags([mesh.ta; mesh.ta; mesh.ta],0,3*mesh.nf,3*mesh.nf);
            S = repmat(TA(1:mesh.nf,1:mesh.nf)*mesh.Nf,3,1); S = S(:);
            
            Nv = full( sparse(I,J,S,mesh.nv,3) );
            Nv = MeshClass.normalize_vf( Nv );
        end



        function [ edges, e2t, t2e, ...
                e2t1, e2t2, ...
                v2e, ie, ...
                ne, nie, inner_edges, bv, bf] = nc_data( mesh )
            T = double( mesh.faces );
            
            
            I = [T(:,2);T(:,3);T(:,1)];
            J = [T(:,3);T(:,1);T(:,2)];
            S = [1:mesh.nf,1:mesh.nf,1:mesh.nf];
            E = sparse(I,J,S,mesh.nv,mesh.nv);
            
            Elisto = [I,J];
            sElist = sort(Elisto,2);
            s = (MeshClass.normv(Elisto - sElist) > 1e-12);
            t = S'.*(-1).^s;
            [edges,une] = unique(sElist, 'rows');
            
            ne = size(edges,1);
            e2t = zeros(ne,4);
            t2e = zeros(mesh.nf,3);
            ie = zeros(ne,1);
            for m=1:length(edges)
                i = edges(m,1); j = edges(m,2);
                t1 = t(une(m));
                t2 = -(E(i,j) + E(j,i) - abs(t1))*sign(t1);
                e2t(m,1:2) = [t1, t2];
                f = T(abs(t1),:); loc = find(f == (sum(f) - i - j));
                t2e(abs(t1),loc) = m*sign(t1);
                e2t(m,3) = loc;
                if t2 ~= 0
                    f = T(abs(t2),:); loc = find(f == (sum(f) - i - j));
                    t2e(abs(t2),loc) = m*sign(t2);
                    e2t(m,4) = loc;
                    ie(m) = 1;
                end
            end

            v2e = sparse(edges(:,1),edges(:,2),1:length(edges),mesh.nv,mesh.nv);            
            
            ne = size(edges,1);
            nie = sum(ie);
            inner_edges = find(ie);
            bv = zeros(mesh.nv,1);
            bv(edges(ie == 0,:)) = 1;
            bf = zeros(mesh.nf,1);
            bf(sum(ismember(mesh.faces, find(bv==1)),2)>0) = 1;
            
            t1 = abs(e2t(inner_edges,1)); t2 = abs(e2t(inner_edges,2));
            
            I = 1:2*nie;
            S = ones(2*nie,1);
            e2t1 = sparse(I, [t1; t1+mesh.nf], S, 2*nie, 2*mesh.nf);
            e2t2 = sparse(I, [t2; t2+mesh.nf], S, 2*nie, 2*mesh.nf);
        end

        function normalize_mesh( mesh, bbdO )
            if nargin<2
                bbdO = 1;
            end
            xx = mesh.vertices - mean(mesh.vertices);
            bbd = norm(max(mesh.vertices) - min(mesh.vertices));
            xx = mesh.vertices  / bbd * bbdO ;
            mesh.vertices = xx;

            mesh.compute_all;
        end


        function aa = center_mesh( mesh , aa)
            if nargin>1
                xx = mesh.vertices - aa;
            else
                xx = mesh.vertices - mean(mesh.vertices);
                aa = mean(mesh.vertices);
            end
            mesh.vertices = xx;

            mesh.compute_all;
        end

        function scale_mesh( mesh, scale_fac )
            
            xx = mesh.vertices * scale_fac ;
            mesh.vertices = xx;

            mesh.compute_all;
        end


        function p = visualizeMesh( mesh, f_vertices, f_faces, edgeColorFlag, figFlag )
            if nargin < 2
                f_vertices = []; f_faces = []; edgeColorFlag = 0; figFlag = 1;
            elseif nargin <3
                f_faces = []; edgeColorFlag = 0; figFlag = 1;
            elseif nargin <4
                edgeColorFlag = 0; figFlag = 1;
            elseif nargin <5
                figFlag = 1;
            end
            
            f_vertices = f_vertices(:); f_faces = f_faces(:);
            
            if figFlag
                figure;
            end
            p = patch('Faces',mesh.faces,'Vertices',mesh.vertices);
            if ~isempty(f_vertices) && ~isempty(f_faces)
                disp('Cannot display f_vertices and f_faces');
            elseif ~isempty(f_vertices) % && isempty(f_faces)
                p.FaceVertexCData = f_vertices; p.FaceColor = 'interp'; colorbar;
            elseif ~isempty(f_faces) % && isempty(f_vertices)
                p.FaceVertexCData = f_faces; p.FaceColor = 'flat'; colorbar;
            else  % isempty(f_vertices) && isempty(f_faces)
                p.FaceColor = 'w'; p.FaceAlpha = 0.5; title(mesh.name);
            end
            if edgeColorFlag ==1
                p.EdgeColor = 'none';
            end
            p.FaceAlpha = 1;
            axis equal; axis off;
            cameratoolbar; cameratoolbar('SetCoordSys','none');
        end
        
        function p = vectorFieldVisualization( mesh , ...
            vectorField , vectorFieldPos, f_vertices, f_faces, edgeColorFlag )
            
            if size(vectorField,2)==1
                vectorField = reshape(vectorField,[],3);
            end
            if nargin<3
                vectorFieldPos = []; f_vertices = []; f_faces = []; edgeColorFlag = 0;
            elseif nargin<4
                f_vertices = []; f_faces = []; edgeColorFlag = 0;
            end

            if isempty(vectorFieldPos)
               vectorFieldPos = mesh.baryCentersCalc; 
            end
            p = mesh.visualizeMesh( f_vertices,f_faces,edgeColorFlag );
            hold on;
            quiver3(vectorFieldPos(:,1),vectorFieldPos(:,2),vectorFieldPos(:,3),...
                vectorField(:,1),vectorField(:,2),vectorField(:,3));
        end
        
        function p = vectorFieldVisualization2( mesh , ...
            vectorField1 , vectorField2 , vectorFieldPos, f_vertices, f_faces, edgeColorFlag )

            p = mesh.visualizeMesh( f_vertices,f_faces,edgeColorFlag );
            hold on;
            quiver3(vectorFieldPos(:,1),vectorFieldPos(:,2),vectorFieldPos(:,3),...
                vectorField1(:,1),vectorField1(:,2),vectorField1(:,3),'b');
            quiver3(vectorFieldPos(:,1),vectorFieldPos(:,2),vectorFieldPos(:,3),...
                vectorField2(:,1),vectorField2(:,2),vectorField2(:,3),'r');
        end

        function p = visualizeDistances( mesh, u, x0, nisolines, urange, cam)
            if nargin < 6
                cam = [];
                if nargin < 5
                    urange = [min(u), max(u)];
                    if nargin < 4
                        nisolines = 0;
                    end
                end
            end

            p = mesh.visualizeMesh(u,[],1,1);
            caxis([urange(1), urange(2)]);
            % colorbar off;
            hold on;
            scatter3(mesh.vertices(x0,1),mesh.vertices(x0,2),mesh.vertices(x0,3),'r','filled');
            if nisolines > 0
                colormap(lines(2*(nisolines-1)));
                colormap jet;
                % [LS,LD,I] = isolines(Mm.vertices,Mm.faces,f,linspace(umin,umax,nisolines)); % use gptoolbox to display isoline
                % plot3([LS(:,1) LD(:,1)]',[LS(:,2) LD(:,2)]', [LS(:,3) LD(:,3)]','Color','k','LineWidth',4);
            else
                colormap jet;
            end

            if ~isempty(cam); MeshClass.set_camera(gca, cam); end

            l = camlight; lighting phong; material dull;
        end


        % op = oph' * oph so can be used with lsqnonlin
        function [ op, oph ] = godf( mesh, n )
            
            M = mesh;
            
            inner = M.inner_edges;
            
            t1 = abs(M.e2t(inner,1)); t2 = abs(M.e2t(inner,2));
            oe = -ones(M.nie,1); ze = zeros(M.nie,1);
            
            % Connection-specific computations
            EV = M.vertices(M.edges(inner,2),:) - M.vertices(M.edges(inner,1),:);
            EV = MeshClass.normalize_vf( EV );
            
            IN1 = atan2(dot(EV,M.F2(t1,:),2),dot(EV,M.F1(t1,:),2));
            IN2 = atan2(dot(EV,M.F2(t2,:),2),dot(EV,M.F1(t2,:),2));
            PT = n*(IN2-IN1);
            
            II = repmat((1:M.nie)',2,1); II = repmat([II; II+M.nie],2,1);
            JJ = [t1; t1+M.nf; t1; t1+M.nf; ...
                  t2; t2+M.nf; t2; t2+M.nf];
            SS = [cos(PT); -sin(PT); sin(PT); cos(PT); oe; ze; ze; oe]; 
            CovD = sparse(II,JJ,SS,2*M.nie,2*M.nf);
            
            Ws = spdiags(repmat(sqrt(M.ea(inner)),2,1), 0, 2*M.nie, 2*M.nie );
            oph = Ws * CovD;
            op = oph' * oph;

        end


    end

    methods (Static)
        function [ nv ] = normv( vf )
            nv = sqrt(sum(vf.^2,2));
        end

        function [ nnv ] = normalize_vf( vf )
            nnv = vf ./ repmat(MeshClass.normv(vf),1,3);
            nnv(MeshClass.normv(vf)<1e-15,:) = 0;
        end

        function cam = get_camera(ca)
            if nargin<1
                ca = gca;
            end

            cam.pba = get(ca, 'PlotBoxAspectRatio');
            cam.dar = get(ca, 'DataAspectRatio');
            cam.cva = get(ca, 'CameraViewAngle');
            cam.cuv = get(ca, 'CameraUpVector');
            cam.ct = get(ca, 'CameraTarget');
            cam.cp = get(ca, 'CameraPosition');
        end

        function set_camera(ca,cam)
            set(ca, 'PlotBoxAspectRatio',cam.pba);
            set(ca, 'DataAspectRatio',cam.dar);
            set(ca, 'CameraViewAngle',cam.cva);
            set(ca, 'CameraUpVector',cam.cuv);
            set(ca, 'CameraTarget',cam.ct);
            set(ca, 'CameraPosition',cam.cp);
        end
    end
end
        