clear all; close all;
%%
H_List=0.095:-0.005:0.03;%0.19:-0.01:0.11; % 1.5:-0.1:0.1;
%%
for Hmax=H_List
f=fopen(['naive_fe_ddm_',num2str(Hmax),'.txt'],'w');
model = createpde;
% rect domain x 4 
R1 = [3,4,-5.5,0.5,0.5,-5.5,-5.5,-5.5,0.5,0.5]'; % Rectangle(square) subdomain SW
R2 = [3,4,-0.5,5.5,5.5,-0.5,-5.5,-5.5,0.5,0.5]'; % Rectangle(square) subdomain SE 
R3 = [3,4,-0.5,5.5,5.5,-0.5,-0.5,-0.5,5.5,5.5]'; % Rectangle(square) subdomain NE
R4 = [3,4,-5.5,0.5,0.5,-5.5,-0.5,-0.5,5.5,5.5]'; % Rectangle(square) subdomain NW
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
%% generate information
Num_Nodes = size(mesh.Nodes,2);
Num_Elements = size(mesh.Elements,2);
fprintf(f,'%d %d\n',[Num_Nodes,Num_Elements]);
Element2NodeList = mesh.Elements';
GlobalNodeCoord = mesh.Nodes';
for i=1:Num_Nodes
    fprintf(f,'%f %f ',GlobalNodeCoord(i,:));
end
fprintf(f,'\n');
for i=1:Num_Elements
    fprintf(f,'%d %d %d ',Element2NodeList(i,:)-1);
end
fprintf(f,'\n');
eps = ones(1,Num_Elements); % relative permittivity of the elements
%% find source
Node_Source = [findNodes(mesh,'nearest',[-4;4]),findNodes(mesh,'nearest',[-4;-4]),findNodes(mesh,'nearest',[4;4]),findNodes(mesh,'nearest',[4;-4]) ...
    , findNodes(mesh,'nearest',[-1;1]),findNodes(mesh,'nearest',[-1;-1]),findNodes(mesh,'nearest',[1;1]),findNodes(mesh,'nearest',[1;-1])];
for i=1:length(Node_Source)
    fprintf(f,'%d ',Node_Source(i)-1);
end
fprintf(f,'\n');

Elements_Node_Source = findElements(mesh,'attached',Node_Source);
for i=1:length(Elements_Node_Source)
    fprintf(f,'%d ',Elements_Node_Source(i)-1);
end
fprintf(f,'\n');
% Node_Source_SubDomain1 = ()
%% find nodes on the Boundary and the Interface of subdomains 
% Find the nodes associated with subdomain 1
Nodes_SubDomain1 = findNodes(mesh,'region','Face',[4,5,6,7]);
% Find the nodes associated with subdomain 2
Nodes_SubDomain2 = findNodes(mesh,'region','Face',[3,5,7,8]);
% Find the nodes associated with subdomain 3
Nodes_SubDomain3 = findNodes(mesh,'region','Face',[2,5,8,9]);
% Find the nodes associated with subdomain 4
Nodes_SubDomain4 = findNodes(mesh,'region','Face',[1,4,5,9]);

% when counting the artificial dirichlet boundary condition we should be
% careful about the repetitive counts
% Find the nodes associated with the Interface of subdomain 1.
Nodes_Interface_SubDomain41 = findNodes(mesh,'region','Edge',[5]);
Nodes_Interface_SubDomain31 = findNodes(mesh,'region','Edge',[6]);
Nodes_Interface_SubDomain21 = findNodes(mesh,'region','Edge',[30,13]);
Nodes_Interface_SubDomain31 = setdiff(Nodes_Interface_SubDomain31,Nodes_Interface_SubDomain41);
Nodes_Interface_SubDomain21 = setdiff(Nodes_Interface_SubDomain21,[Nodes_Interface_SubDomain41,Nodes_Interface_SubDomain31]);
% Find the nodes associated with the boundary of subdomain 2.
Nodes_Interface_SubDomain42 = findNodes(mesh,'region','Edge',[10]);
Nodes_Interface_SubDomain32 = findNodes(mesh,'region','Edge',[26,6]);
Nodes_Interface_SubDomain12 = findNodes(mesh,'region','Edge',[27]);
Nodes_Interface_SubDomain32 = setdiff(Nodes_Interface_SubDomain32,Nodes_Interface_SubDomain42);
Nodes_Interface_SubDomain42 = setdiff(Nodes_Interface_SubDomain42,[Nodes_Interface_SubDomain32,Nodes_Interface_SubDomain12]);
% Find the nodes associated with the boundary of subdomain 3.
Nodes_Interface_SubDomain43 = findNodes(mesh,'region','Edge',[28,10]);
Nodes_Interface_SubDomain23 = findNodes(mesh,'region','Edge',[29]);
Nodes_Interface_SubDomain13 = findNodes(mesh,'region','Edge',[12]);
Nodes_Interface_SubDomain13 = setdiff(Nodes_Interface_SubDomain13,Nodes_Interface_SubDomain23);
Nodes_Interface_SubDomain43 = setdiff(Nodes_Interface_SubDomain43,[Nodes_Interface_SubDomain23,Nodes_Interface_SubDomain13]);
% Find the nodes associated with the boundary of subdomain 4.
Nodes_Interface_SubDomain34 = findNodes(mesh,'region','Edge',[31]);
Nodes_Interface_SubDomain24 = findNodes(mesh,'region','Edge',[13]);
Nodes_Interface_SubDomain14 = findNodes(mesh,'region','Edge',[11,12]);
Nodes_Interface_SubDomain24 = setdiff(Nodes_Interface_SubDomain24,Nodes_Interface_SubDomain34);
Nodes_Interface_SubDomain14 = setdiff(Nodes_Interface_SubDomain14,[Nodes_Interface_SubDomain24,Nodes_Interface_SubDomain34]);
% what you need to do about the above indices is to find the local indices
% of each interfaces
% only the FEM coefficient need this local indices, because only they will
% be updated in each interation

% these are all global indices
% Find the nodes associated within subdomain 1 (except those on the interface)
Nodes_InteriorSubDomain1 = setxor(Nodes_SubDomain1,[Nodes_Interface_SubDomain41,Nodes_Interface_SubDomain31,Nodes_Interface_SubDomain21]);
% Find the nodes associated within subdomain 2 (except those on the interface)
Nodes_InteriorSubDomain2 = setxor(Nodes_SubDomain2,[Nodes_Interface_SubDomain42,Nodes_Interface_SubDomain32,Nodes_Interface_SubDomain12]);
% Find the nodes associated within subdomain 3 (except those on the interface)
Nodes_InteriorSubDomain3 = setxor(Nodes_SubDomain3,[Nodes_Interface_SubDomain43,Nodes_Interface_SubDomain23,Nodes_Interface_SubDomain13]);
% Find the nodes associated within subdomain 4 (except those on the interface)
Nodes_InteriorSubDomain4 = setxor(Nodes_SubDomain4,[Nodes_Interface_SubDomain34,Nodes_Interface_SubDomain24,Nodes_Interface_SubDomain14]);


% we only count the interior as the node list of the subdomain, those on the boundaries are only used for enforcing the RHS
% send: select from the subdomain local phi array
% recv: use appropriative array * 3
Nodes_InteriorSubDomain_List={Nodes_InteriorSubDomain1,Nodes_InteriorSubDomain2,Nodes_InteriorSubDomain3,Nodes_InteriorSubDomain4};
for ii=1:4
    for i=1:length(Nodes_InteriorSubDomain_List{ii})
        fprintf(f,'%d ',Nodes_InteriorSubDomain_List{ii}(i)-1);
    end
    fprintf(f,'\n');
end

Nodes_InterfaceSubDomain_List={Nodes_Interface_SubDomain41,Nodes_Interface_SubDomain31,Nodes_Interface_SubDomain21;
    Nodes_Interface_SubDomain42,Nodes_Interface_SubDomain32,Nodes_Interface_SubDomain12;
    Nodes_Interface_SubDomain43,Nodes_Interface_SubDomain23,Nodes_Interface_SubDomain13;
    Nodes_Interface_SubDomain34,Nodes_Interface_SubDomain24,Nodes_Interface_SubDomain14};
for iii=1:4
    for ii=1:3
        for i=1:length(Nodes_InterfaceSubDomain_List{iii,ii})
            fprintf(f,'%d ',Nodes_InterfaceSubDomain_List{iii,ii}(i)-1);
        end
        fprintf(f,'\n');
    end
end


Nodes_Interface_InWhichSubDomain=[4,3,2;4,3,1;4,2,1;3,2,1];
% the subdomain ABC to local array: store the local indices that resides others boundary
ThisSubABCNodes2ThatSubLocalNodes = {0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;}; % 4 subdomains x 3 subdomains
for this_dm=1:4
    for those_dm=1:3
        ABC_InWhichSubDomain = Nodes_Interface_InWhichSubDomain(this_dm,those_dm);
        this2thatLocal = zeros(size(Nodes_InterfaceSubDomain_List{this_dm,those_dm}));
        for i=1:length(Nodes_InterfaceSubDomain_List{this_dm,those_dm})%length(Nodes_InteriorSubDomain_List{ABC_InWhichSubDomain})
            this2thatLocal(i) = find(Nodes_InteriorSubDomain_List{ABC_InWhichSubDomain}==Nodes_InterfaceSubDomain_List{this_dm,those_dm}(i)); % the local node index is the index of the global index array 
        end
        ThisSubABCNodes2ThatSubLocalNodes{this_dm,Nodes_Interface_InWhichSubDomain(this_dm,those_dm)}=this2thatLocal;
    end
end

% Find the nodes associated with edges 1, 2 and 10 (the boundary of triangle object).
Nodes_onObj_Tri = findNodes(mesh,'region','Edge',[1,2,17]);
for i=1:length(Nodes_onObj_Tri)
    fprintf(f,'%d ',Nodes_onObj_Tri(i)-1);
end
fprintf(f,'\n');
% Find the nodes associated with edges 3, 4, 11 and 12 (the boundary of square object).
Nodes_onObj_Square = findNodes(mesh,'region','Edge',[3, 4, 18, 19]);
for i=1:length(Nodes_onObj_Square)
    fprintf(f,'%d ',Nodes_onObj_Square(i)-1);
end
fprintf(f,'\n');
% Find the nodes associated with edges 20, 21, 22 and 23 (the boundary of circle object).
Nodes_onObj_Circle = findNodes(mesh,'region','Edge',[32, 33, 34, 35]);
for i=1:length(Nodes_onObj_Circle)
    fprintf(f,'%d ',Nodes_onObj_Circle(i)-1);
end
fprintf(f,'\n');
% Find the nodes associated with edges 24, 25, 26, 27 and 28 (the boundary of ellipse object).
Nodes_onObj_Ellipse = findNodes(mesh,'region','Edge',[36, 37, 38, 39, 40]);
for i=1:length(Nodes_onObj_Ellipse)
    fprintf(f,'%d ',Nodes_onObj_Ellipse(i)-1);
end
fprintf(f,'\n');
%% find elements on different subdomains 
% Find the elements associated with subdomain 1 
Elements_SubDomain1 = findElements(mesh,'region','Face',[4,5,6,7]);
% Find the elements associated with subdomain 2 
Elements_SubDomain2 = findElements(mesh,'region','Face',[3,5,7,8]);
% Find the elements associated with subdomain 3 
Elements_SubDomain3 = findElements(mesh,'region','Face',[2,5,8,9]);
% Find the elements associated with subdomain 4 
Elements_SubDomain4 = findElements(mesh,'region','Face',[1,4,5,9]);

%% build subdomains local index
Num_Elements1 = size(Elements_SubDomain1,2);
Num_Elements2 = size(Elements_SubDomain2,2);
Num_Elements3 = size(Elements_SubDomain3,2);
Num_Elements4 = size(Elements_SubDomain4,2);

Num_Elements_List = {Num_Elements1,Num_Elements2,Num_Elements3,Num_Elements4};
for i=1:4
        fprintf(f,'%d ',Num_Elements_List{i});
end
fprintf(f,'\n');
Elements_List = {Elements_SubDomain1,Elements_SubDomain2,Elements_SubDomain3,Elements_SubDomain4};
for ii=1:4
    for i=1:length(Elements_List{ii})
        fprintf(f,'%d ',Elements_List{ii}(i)-1);
    end
    fprintf(f,'\n');
end

Nodes_onBoundary = findNodes(mesh,'region','Edge',[23,25,24,7,8,9,14,15,16,20,21,22]);
Nodes_Dirichlet = [Nodes_onBoundary,Nodes_onObj_Circle,Nodes_onObj_Ellipse,Nodes_onObj_Square,Nodes_onObj_Tri];

for i=1:length(Nodes_Dirichlet)
        fprintf(f,'%d ',Nodes_Dirichlet(i)-1);
end
fclose(f);
end
fclose('all');