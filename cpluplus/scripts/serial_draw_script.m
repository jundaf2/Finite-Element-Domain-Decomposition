clear all; close all;
N = 100; % interpolation size
Hmax=0.07; % element size

model = createpde;
% rect domain x 4 
R1 = [3,4,-5,0.5,0.5,-5,-5,-5,0.5,0.5]'; % Rectangle(square) subdomain SW
R2 = [3,4,-0.5,5,5,-0.5,-5,-5,0.5,0.5]'; % Rectangle(square) subdomain SE 
R3 = [3,4,-0.5,5,5,-0.5,-0.5,-0.5,5,5]'; % Rectangle(square) subdomain NE
R4 = [3,4,-5,0.5,0.5,-5,-0.5,-0.5,5,5]'; % Rectangle(square) subdomain NW
% obj x 4 
C1 = [1,-2.5,-2.5,1]'; % Circle object in R1
C1 = [C1;zeros(length(R1) - length(C1),1)]; % Append extra zeros to the circles so they have the same number of rows as the rectangle
E1 = [4,2.5,-2.5,1.5,1,0]'; % Ellipse object in R2
E1 = [E1;zeros(length(R1) - length(E1),1)];
r1 = [3,4,1.5,3.5,3.5,1.5,1.5,1.5,3.5,3.5]'; % Rectangle object (square) in R3
T1 = [2,3,-3.5,-1.5,-2.5,2.5-sqrt(3)/2,2.5-sqrt(3)/2,2.5+sqrt(3)/2]'; % Triangle object in R4
T1 = [T1;zeros(length(R1) - length(T1),1)];

gm = [R1,R2,R3,R4,C1,E1,T1,r1];
sf = 'R1-C1+R2-E1+R3-r1+R4-T1';
ns = char('R1','R2','R3','R4','C1','E1','T1','r1');
ns = ns';
g = decsg(gm,sf,ns);
geometryFromEdges(model,g);
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',Hmax);
figure;
hold on;
pdegplot(model,'FaceLabels','on','EdgeLabels','on')
xlim([-6,6])
ylim([-6,6])
xlabel('x');
xlabel('y');
hold off;
Num_Nodes = size(mesh.Nodes,2);
Num_Elements = size(mesh.Elements,2);

%%
% Find the elements associated with subdomain 1 
Elements_SubDomain1 = findElements(mesh,'region','Face',[4,5,6,7]);
% Find the elements associated with subdomain 2 
Elements_SubDomain2 = findElements(mesh,'region','Face',[3,5,7,8]);
% Find the elements associated with subdomain 3 
Elements_SubDomain3 = findElements(mesh,'region','Face',[2,5,8,9]);
% Find the elements associated with subdomain 4 
Elements_SubDomain4 = findElements(mesh,'region','Face',[1,4,5,9]);
Num_Elements1 = size(Elements_SubDomain1,2);
Num_Elements2 = size(Elements_SubDomain2,2);
Num_Elements3 = size(Elements_SubDomain3,2);
Num_Elements4 = size(Elements_SubDomain4,2);
Elements_List = {Elements_SubDomain1,Elements_SubDomain2,Elements_SubDomain3,Elements_SubDomain4};
Num_Elements_List = {Num_Elements1,Num_Elements2,Num_Elements3,Num_Elements4};

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


%%
file_read = ['naive_fe_ddm_',num2str(Hmax),'out.txt'];
first_line = dlmread(file_read, ' ');
phi = zeros(Num_Nodes,1);
for i=2:length(first_line)
    phi(i-1)=first_line(i);
end
%%
Nel = @(x,y,ael,bel,cel,Delta_e) 1/(2*Delta_e)*(ael+bel*x+cel*y);
X=linspace(-5,5,N);
Y=linspace(-5,5,N);
Interp_Points_List_X = {[1,N*0.55],[N*0.45,N],[N*0.45,N],[1,N*0.55]};
Interp_Points_List_Y = {[1,N*0.55],[1,N*0.55],[N*0.45,N],[N*0.45,N]};
Phi = zeros(N,N);
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

%%
figure; hold on;
title('\Phi');
h1 = imagesc(X,Y(end:-1:1),rot90(Phi));
set(gca,'fontsize',24);
axis tight;
axis equal; hold off;