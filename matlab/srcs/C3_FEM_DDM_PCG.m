%%
function [Phi,timeVal,error_list] = C3_FEM_DDM_PCG(mesh,N_,plot_logistics)
%% functions, macros and lambda expressions
Nel = @(x,y,ael,bel,cel,Delta_e) 1/(2*Delta_e)*(ael+bel*x+cel*y); 

globalMesh = mesh; 
% [p,e,t] = meshToPet(mesh);
Num_Nodes = size(globalMesh.Nodes,2);
Num_Elements = size(globalMesh.Elements,2);

eps = ones(1,Num_Elements); % relative permittivity of the elements


%% find source
Node_Source = [findNodes(globalMesh,'nearest',[-4;4]),findNodes(globalMesh,'nearest',[-4;-4]),findNodes(globalMesh,'nearest',[4;4]),findNodes(globalMesh,'nearest',[4;-4]) ...
    , findNodes(globalMesh,'nearest',[-1;1]),findNodes(globalMesh,'nearest',[-1;-1]),findNodes(globalMesh,'nearest',[1;1]),findNodes(globalMesh,'nearest',[1;-1])];
%% find nodes on the Boundary and the Interface of subdomains 
% Find the nodes associated with subdomain 1
Nodes_SubDomain1 = findNodes(globalMesh,'region','Face',[4,5,6,7]);
% Find the nodes associated with subdomain 2
Nodes_SubDomain2 = findNodes(globalMesh,'region','Face',[3,5,7,8]);
% Find the nodes associated with subdomain 3
Nodes_SubDomain3 = findNodes(globalMesh,'region','Face',[2,5,8,9]);
% Find the nodes associated with subdomain 4
Nodes_SubDomain4 = findNodes(globalMesh,'region','Face',[1,4,5,9]);
% Find the nodes associated with edges 1, 2 and 10 (the boundary of triangle object).
Nodes_onObj_Tri = findNodes(globalMesh,'region','Edge',[1,2,17]);
% Find the nodes associated with edges 3, 4, 11 and 12 (the boundary of square object).
Nodes_onObj_Square = findNodes(globalMesh,'region','Edge',[3, 4, 18, 19]);
% Find the nodes associated with edges 20, 21, 22 and 23 (the boundary of circle object).
Nodes_onObj_Circle = findNodes(globalMesh,'region','Edge',[32, 33, 34, 35]);
% Find the nodes associated with edges 24, 25, 26, 27 and 28 (the boundary of ellipse object).
Nodes_onObj_Ellipse = findNodes(globalMesh,'region','Edge',[36, 37, 38, 39, 40]);

%% find elements on different subdomains 
% Find the elements associated with subdomain 1 
Elements_SubDomain1 = findElements(globalMesh,'region','Face',[4,5,6,7]);
% Find the elements associated with subdomain 2 
Elements_SubDomain2 = findElements(globalMesh,'region','Face',[3,5,7,8]);
% Find the elements associated with subdomain 3 
Elements_SubDomain3 = findElements(globalMesh,'region','Face',[2,5,8,9]);
% Find the elements associated with subdomain 4 
Elements_SubDomain4 = findElements(globalMesh,'region','Face',[1,4,5,9]);

%% build subdomains local index
Num_Elements1 = size(Elements_SubDomain1,2);

Num_Elements2 = size(Elements_SubDomain2,2);

Num_Elements3 = size(Elements_SubDomain3,2);

Num_Elements4 = size(Elements_SubDomain4,2);

Num_Elements_List = {Num_Elements1,Num_Elements2,Num_Elements3,Num_Elements4};
Elements_List = {Elements_SubDomain1,Elements_SubDomain2,Elements_SubDomain3,Elements_SubDomain4};
%% Plot mesh geometry
mesh = globalMesh;
if plot_logistics
    figure;hold on;
    % pdeplot(model); 
    xlabel('x');
    xlabel('y');
    %title('Finite element discretization using first-order triangular elements')
    pdemesh(mesh,'ElementLabels','off')
    hold on
    pdemesh(mesh.Nodes,mesh.Elements(:,Elements_SubDomain1),'EdgeColor','blue')
    pdemesh(mesh.Nodes,mesh.Elements(:,Elements_SubDomain2),'EdgeColor','magenta')
    pdemesh(mesh.Nodes,mesh.Elements(:,Elements_SubDomain4),'EdgeColor','green')
    pdemesh(mesh.Nodes,mesh.Elements(:,Elements_SubDomain3),'EdgeColor','red')
set(gca,'fontsize',24);
    xlim([-6,6])
    ylim([-6,6])
    hold off;
end

%% other initializations
Nodes_onBoundary = findNodes(mesh,'region','Edge',[23,25,24,7,8,9,14,15,16,20,21,22]);
Num_Unknown_Nodes = Num_Nodes;%Num_Nodes - Num_Nodes_onBoundary;
N = N_;
X=linspace(-5,5,N);
Y=linspace(-5,5,N);
Interp_Points_List_X = {[1,N*0.55],[N*0.45,N],[N*0.45,N],[1,N*0.55]};
Interp_Points_List_Y = {[1,N*0.55],[1,N*0.55],[N*0.45,N],[N*0.45,N]};
Phi = zeros(N,N);

tic;
K = zeros(Num_Unknown_Nodes,Num_Unknown_Nodes);
Ke = zeros(3,3);
b = zeros(Num_Unknown_Nodes,1);


%% build FEM LHS matrix [K] and RHS vector [b]
for e=1:Num_Elements
    who_am_I = mesh.Elements(1,e); % the global node idx of node 1 in element e
    where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 1 in element e
    xe(1) = where_am_I(1);
    ye(1) = where_am_I(2);
    ne(1) = who_am_I;
    who_am_I = mesh.Elements(2,e); % the global node idx of node 2 in element e
    where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 2 in element e
    xe(2) = where_am_I(1);
    ye(2) = where_am_I(2);
    ne(2) = who_am_I;
    who_am_I = mesh.Elements(3,e); % the global node idx of node 3 in element e
    where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 3 in element e
    xe(3) = where_am_I(1);
    ye(3) = where_am_I(2);
    ne(3) = who_am_I;
    
    be(1) = ye(2) - ye(3);% y 2-3
    be(2) = ye(3) - ye(1);% y 3-1
    be(3) = ye(1) - ye(2);% y 1-2
    ce(1) = xe(3) - xe(2);% x 3-2
    ce(2) = xe(1) - xe(3);% x 1-3
    ce(3) = xe(2) - xe(1);% x 2-1
    Delta_e = (be(1)*ce(2)-be(2)*ce(1))/2; % area of element e 
    
    % Generate the elemental matrix [Ke]
    for i=1:3
        for j=1:3
            Ke(i,j) = eps(e)/(4*Delta_e) * (be(i)*be(j) + ce(i)*ce(j)); 
        end
    end
    % Add [Ke] to [K]
    for i=1:3
        for j=1:3
            K(ne(i),ne(j)) = K(ne(i),ne(j))+Ke(i,j); 
        end
    end    
end

Nodes_Dirichlet = [Nodes_onBoundary,Nodes_onObj_Circle,Nodes_onObj_Ellipse,Nodes_onObj_Square,Nodes_onObj_Tri];
% Impose the Dirichlet boundary condition
for i=1:length(Nodes_Dirichlet)
    b(Nodes_Dirichlet(i))=0;
    K(Nodes_Dirichlet(i),Nodes_Dirichlet(i))=1;
    for j=1:Num_Unknown_Nodes
        if j==Nodes_Dirichlet(i), continue; end
        b(j)=b(j)-K(j,Nodes_Dirichlet(i))*0;
        K(Nodes_Dirichlet(i),j)=0;
        K(j,Nodes_Dirichlet(i))=0;
    end
end
% Impose Source
Elements_Node_Source = findElements(mesh,'attached',Node_Source);
for i=1:length(Elements_Node_Source)
    who_am_I = mesh.Elements(1,Elements_Node_Source(i)); % the global node idx of node 1 in Elements_Node_Source(i)
    where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 1 in element Elements_Node_Source(i)
    xe(1) = where_am_I(1);
    ye(1) = where_am_I(2);
    ne(1) = who_am_I;
    who_am_I = mesh.Elements(2,Elements_Node_Source(i)); % the global node idx of node 2 in Elements_Node_Source(i)
    where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 2 in Elements_Node_Source(i)
    xe(2) = where_am_I(1);
    ye(2) = where_am_I(2);
    ne(2) = who_am_I;
    who_am_I = mesh.Elements(3,Elements_Node_Source(i)); % the global node idx of node 3 in Elements_Node_Source(i)
    where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 3 in Elements_Node_Source(i)
    xe3 = where_am_I(1);
    ye(3) = where_am_I(2);
    ne(3) = who_am_I;
    
    be(1) = ye(2) - ye(3);% y 2-3
    be(2) = ye(3) - ye(1);% y 3-1
    be(3) = ye(1) - ye(2);% y 1-2
    ce(1) = xe3 - xe(2);% x 3-2
    ce(2) = xe(1) - xe3;% x 1-3
    ce(3) = xe(2) - xe(1);% x 2-1
    Delta_e = (be(1)*ce(2)-be(2)*ce(1))/2;
    
    for j=1:length(Node_Source)
        if ismember(Node_Source(j),ne)
            b(Node_Source(j))=b(Node_Source(j))+1*1/3*Delta_e;
        end
    end
end
%% precondition (Incomplete LU factorization)
L = {0,0;0,0};
U = {0,0;0,0};
P = {0,0;0,0};
for col = 1:2
    for row = 1:2 
        subdomain_idx = eval(['Nodes_SubDomain',num2str((col-1)*2+row)]);        
        A = K(subdomain_idx, subdomain_idx); % the sub FEM matrix
        n = size(A,2);
        L{col,row} = eye(n);
        P{col,row} = zeros(n,n);
        U{col,row} = zeros(n,n);

        q = 1:n;
        
        for i=1:n-1
            % row interchanges (pivoting)
            [~, pivot]=max(abs(A(i:n,i)));
            pivot = pivot + i - 1;
            if pivot ~= i
                temp = A(i,:);
                A(i,:) = A(pivot,:);
                A(pivot,:) = temp;      
                temp = q(i);
                q(i)= q(pivot);
                q(pivot) = temp;
            end
            
            % Gaussian Elimination
            for j=i+1:n
                m = -A(j,i) / A(i,i);
                A(j,i) = -m;
                A(j,i+1:n) = A(j,i+1:n) + m * A(i,i+1:n);
            end
        end
        % force sparsity 
        A(abs(A)<1e-3)=0;
        % Partition (stores)
        for i = 1:n
            for j = 1:n
                if i <= j
                    U{col,row}(i,j) = A(i,j);
                else
                    L{col,row}(i,j) = A(i,j);
                end
            end

            P{col,row}(i, q(i)) = 1;
        end
        % force ILU(0)
        U{col,row}(A==0)=0;
        L{col,row}(A==0)=0;
    end
end

%% solve (Preconditioned CG), start with an initial guess of zero vector
phi = zeros(Num_Unknown_Nodes,1);
max_iter_num = 100; % we do not set threshold for it

norm_b2 = norm(b);
if  ( norm_b2 == 0.0 ), norm_b2 = 1.0; end

r = b - K*phi;
error = norm( r ) / norm_b2;

error_list = zeros(1,max_iter_num+1);
for iter = 1:max_iter_num
   
    z = zeros(Num_Nodes,1);
    for col = 1:2
        for row = 1:2 
            subdomain_idx = eval(['Nodes_SubDomain',num2str((col-1)*2+row)]);
            z(subdomain_idx)=z(subdomain_idx)+U{col,row}^-1*(L{col,row}^-1*P{col,row})*r(subdomain_idx);
        end
    end

    rho = (r'*z);

    if ( iter > 1 )
    beta = rho / rho_1;
    p = z + beta*p;
    else
    p = z;
    end

    q = K*p;
    alpha = rho / (p'*q );

    phi = phi + alpha * p;

    r = r - alpha*q;

    error_list(iter) = error;
    error = norm( r ) / norm_b2;  
    if mod(iter,10)==0
        disp(['Error',' at iteration ', num2str(iter),  ' = ', num2str(error)])
    end
    rho_1 = rho;
end
timeVal = toc;
%% show result, also in a DDM fashion.

for subDomain=1:4
    Num_Elements_Subdomain = Num_Elements_List{subDomain};
    for i=1:Num_Elements_Subdomain
        if mod(i,100)==0
            disp(['Interpolating ',num2str(i),'th element'])
        end

        who_am_I = mesh.Elements(1,Elements_List{subDomain}(i)); % the global node idx of node 1 in Elements(i)
        where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 1 in Elements(i)
        xe1 = where_am_I(1);
        ye1 = where_am_I(2);
        ne1 = who_am_I;

        who_am_I = mesh.Elements(2,Elements_List{subDomain}(i)); % the global node idx of node 1 in Elements(i)
        where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 1 in Elements(i)
        xe2 = where_am_I(1);
        ye2 = where_am_I(2);
        ne2 = who_am_I;

        who_am_I = mesh.Elements(3,Elements_List{subDomain}(i)); % the global node idx of node 1 in Elements(i)
        where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 1 in Elements(i)
        xe3 = where_am_I(1);
        ye3 = where_am_I(2);
        ne3 = who_am_I;

        ae1 = xe2*ye3-xe3*ye2;
        ae2 = xe3*ye1-xe1*ye3;
        ae3 = xe1*ye2-xe2*ye1;
        be1 = ye2 - ye3;% y 2-3
        be2 = ye3 - ye1;% y 3-1
        be3 = ye1 - ye2;% y 1-2
        ce1 = xe3 - xe2;% x 3-2
        ce2 = xe1 - xe3;% x 1-3
        ce3 = xe2 - xe1;% x 2-1
        Delta_e = (be1*ce2-be2*ce1)/2;
        
        start_q = Interp_Points_List_X{subDomain}(1);
        end_q = Interp_Points_List_X{subDomain}(2);
        start_p = Interp_Points_List_Y{subDomain}(1);
        end_p = Interp_Points_List_Y{subDomain}(2);
        for q=start_q:end_q
            for p=start_p:end_p
                [in,on]=inpolygon(X(q),Y(p),[xe1,xe2,xe3],[ye1,ye2,ye3]);

                if (in || on) %indicating if inside the element or on the edge of the element
                    Phi(q,p)=Nel(X(q),Y(p),ae1,be1,ce1,Delta_e)*phi(ne1)+Nel(X(q),Y(p),ae2,be2,ce2,Delta_e)*phi(ne2)+Nel(X(q),Y(p),ae3,be3,ce3,Delta_e)*phi(ne3);
                end
            end
        end
    end
end

end
