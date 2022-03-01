%% define geometry
clear all; close all;
N = 100; % interpolation size
Hmax=0.1; % element size
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

%% generate information
Num_Nodes = size(mesh.Nodes,2);
Num_Elements = size(mesh.Elements,2);
GlobalElement2NodeList = mesh.Elements';
GlobalNodeCoord = mesh.Nodes';

% Find the nodes associated with subdomain 1
Nodes_SubDomain1 = findNodes(mesh,'region','Face',[4,5,6,7]);
% Find the nodes associated with subdomain 2
Nodes_SubDomain2 = findNodes(mesh,'region','Face',[3,5,7,8]);
% Find the nodes associated with subdomain 3
Nodes_SubDomain3 = findNodes(mesh,'region','Face',[2,5,8,9]);
% Find the nodes associated with subdomain 4
Nodes_SubDomain4 = findNodes(mesh,'region','Face',[1,4,5,9]);
Nodes_SubDomain_List = {Nodes_SubDomain1,Nodes_SubDomain2,Nodes_SubDomain3,Nodes_SubDomain4};

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

Nodes_InterfaceSubDomain_List={Nodes_Interface_SubDomain41,Nodes_Interface_SubDomain31,Nodes_Interface_SubDomain21;
    Nodes_Interface_SubDomain42,Nodes_Interface_SubDomain32,Nodes_Interface_SubDomain12;
    Nodes_Interface_SubDomain43,Nodes_Interface_SubDomain23,Nodes_Interface_SubDomain13;
    Nodes_Interface_SubDomain34,Nodes_Interface_SubDomain24,Nodes_Interface_SubDomain14};

% Find the nodes associated within subdomain 1 (except those on the interface)
Nodes_InteriorSubDomain1 = setxor(Nodes_SubDomain1,[Nodes_Interface_SubDomain41,Nodes_Interface_SubDomain31,Nodes_Interface_SubDomain21]);
% Find the nodes associated within subdomain 2 (except those on the interface)
Nodes_InteriorSubDomain2 = setxor(Nodes_SubDomain2,[Nodes_Interface_SubDomain42,Nodes_Interface_SubDomain32,Nodes_Interface_SubDomain12]);
% Find the nodes associated within subdomain 3 (except those on the interface)
Nodes_InteriorSubDomain3 = setxor(Nodes_SubDomain3,[Nodes_Interface_SubDomain43,Nodes_Interface_SubDomain23,Nodes_Interface_SubDomain13]);
% Find the nodes associated within subdomain 4 (except those on the interface)
Nodes_InteriorSubDomain4 = setxor(Nodes_SubDomain4,[Nodes_Interface_SubDomain34,Nodes_Interface_SubDomain24,Nodes_Interface_SubDomain14]);

Nodes_InteriorSubDomain_List={Nodes_InteriorSubDomain1,Nodes_InteriorSubDomain2,Nodes_InteriorSubDomain3,Nodes_InteriorSubDomain4};

LocalInteriorNodesIdx_List = {0,0,0,0}; % the local index of the interior nodes
for this_dm=1:4
        LocalInteriorNodesIdx = zeros(size(Nodes_InteriorSubDomain_List{this_dm}));
        for i=1:length(Nodes_InteriorSubDomain_List{this_dm})%length(Nodes_InteriorSubDomain_List{ABC_InWhichSubDomain})
            LocalInteriorNodesIdx(i) = find(Nodes_SubDomain_List{this_dm}==Nodes_InteriorSubDomain_List{this_dm}(i)); % the local node index is the index of the global index array 
        end
        LocalInteriorNodesIdx_List{this_dm}=LocalInteriorNodesIdx;
end

% for assembly
LocalInterefaceNodesIdx_List = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}}; % the local index of the interface nodes
for this_dm=1:4
    for those_dm=1:3
        LocalInterefaceNodesIdx = zeros(size(Nodes_InterfaceSubDomain_List{this_dm,those_dm}));
        for i=1:length(Nodes_InterfaceSubDomain_List{this_dm,those_dm})%length(Nodes_InteriorSubDomain_List{ABC_InWhichSubDomain})
            LocalInterefaceNodesIdx(i) = find(Nodes_SubDomain_List{this_dm}==Nodes_InterfaceSubDomain_List{this_dm,those_dm}(i)); % the local node index is the index of the global index array 
        end
        LocalInterefaceNodesIdx_List{this_dm}{those_dm}=LocalInterefaceNodesIdx;
    end
end


% for send
Nodes_Interface_InWhichSubDomain=[4,3,2;4,3,1;4,2,1;3,2,1];
% the subdomain ABC to local array: store the local indices that resides others boundary
ThisSubABCNodes2ThatSubLocalNodes = {0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;}; % 4 subdomains x 3 subdomains
for this_dm=1:4
    for those_dm=1:3
        ABC_InWhichSubDomain = Nodes_Interface_InWhichSubDomain(this_dm,those_dm);
        this2thatLocal = zeros(size(Nodes_InterfaceSubDomain_List{this_dm,those_dm}));
        for i=1:length(Nodes_InterfaceSubDomain_List{this_dm,those_dm})%length(Nodes_InteriorSubDomain_List{ABC_InWhichSubDomain})
            this2thatLocal(i) = find(Nodes_SubDomain_List{ABC_InWhichSubDomain}==Nodes_InterfaceSubDomain_List{this_dm,those_dm}(i)); % the local node index is the index of the global index array 
        end
        ThisSubABCNodes2ThatSubLocalNodes{this_dm,Nodes_Interface_InWhichSubDomain(this_dm,those_dm)}=this2thatLocal;
    end
end

% Find the nodes associated with edges 1, 2 and 10 (the boundary of triangle object).
Nodes_onObj_Tri = findNodes(mesh,'region','Edge',[1,2,17]);
% Find the nodes associated with edges 3, 4, 11 and 12 (the boundary of square object).
Nodes_onObj_Square = findNodes(mesh,'region','Edge',[3, 4, 18, 19]);
% Find the nodes associated with edges 20, 21, 22 and 23 (the boundary of circle object).
Nodes_onObj_Circle = findNodes(mesh,'region','Edge',[32, 33, 34, 35]);
% Find the nodes associated with edges 24, 25, 26, 27 and 28 (the boundary of ellipse object).
Nodes_onObj_Ellipse = findNodes(mesh,'region','Edge',[36, 37, 38, 39, 40]);

% Find the elements associated with subdomain 1 
Elements_SubDomain1 = findElements(mesh,'region','Face',[4,5,6,7]);
% Find the elements associated with subdomain 2 
Elements_SubDomain2 = findElements(mesh,'region','Face',[3,5,7,8]);
% Find the elements associated with subdomain 3 
Elements_SubDomain3 = findElements(mesh,'region','Face',[2,5,8,9]);
% Find the elements associated with subdomain 4 
Elements_SubDomain4 = findElements(mesh,'region','Face',[1,4,5,9]);

Elements_List = {Elements_SubDomain1,Elements_SubDomain2,Elements_SubDomain3,Elements_SubDomain4};

Node_Source_SubDomain4 = [findNodes(mesh,'nearest',[-4;4]), findNodes(mesh,'nearest',[-1;1])];
Node_Source_SubDomain1 = [findNodes(mesh,'nearest',[-4;-4]), findNodes(mesh,'nearest',[-1;-1])];
Node_Source_SubDomain3 = [findNodes(mesh,'nearest',[4;4]), findNodes(mesh,'nearest',[1;1])];
Node_Source_SubDomain2 = [findNodes(mesh,'nearest',[4;-4]), findNodes(mesh,'nearest',[1;-1])];
Node_Source_SubDomain_List = {Node_Source_SubDomain1,Node_Source_SubDomain2,Node_Source_SubDomain3,Node_Source_SubDomain4};
LocalSourceNodesIdx_List = {0,0,0,0}; % the local index of the source nodes
for this_dm=1:4
        LocalSourceNodesIdx = zeros(size(Node_Source_SubDomain_List{this_dm}));
        for i=1:length(Node_Source_SubDomain_List{this_dm})
            LocalSourceNodesIdx(i) = find(Nodes_SubDomain_List{this_dm}==Node_Source_SubDomain_List{this_dm}(i)); % the local node index is the index of the global index array 
        end
        LocalSourceNodesIdx_List{this_dm}=LocalSourceNodesIdx;
end

% element2node local
LocalElement2LocalNode_List = {0,0,0,0};
for this_dm=1:4
    locallist_with_globalIdx = GlobalElement2NodeList(Elements_List{this_dm},:);
    LocalElement2LocalNode = zeros(size(locallist_with_globalIdx));
    for i=1:size(GlobalElement2NodeList(Elements_List{this_dm},:),2)
        for ii=1:size(GlobalElement2NodeList(Elements_List{this_dm},:),1)
            LocalElement2LocalNode(ii,i) = find(Nodes_SubDomain_List{this_dm}==locallist_with_globalIdx(ii,i)); % the local node index is the index of the global index array 
        end
    end
    LocalElement2LocalNode_List{this_dm}=LocalElement2LocalNode;
end

% local element that contains the node with source
LocalElementsSourceNodesIdx_List = {0,0,0,0}; % the local index of the element
for this_dm=1:4
    GlobalElements_Node_Source = findElements(mesh,'attached',Node_Source_SubDomain_List{this_dm});
    LocalSourceElementsIdx = zeros(size(GlobalElements_Node_Source));
    for i=1:length(GlobalElements_Node_Source)
        LocalSourceElementsIdx(i) = find(Elements_List{this_dm}==GlobalElements_Node_Source(i)); % the local node index is the index of the global index array 
    end
    LocalElementsSourceNodesIdx_List{this_dm}=LocalSourceElementsIdx;
end

Nodes_onBoundary1 = findNodes(mesh,'region','Edge',[23,24,20,21]);
Nodes_onBoundary2 = findNodes(mesh,'region','Edge',[25,24,7,8]);
Nodes_onBoundary3 = findNodes(mesh,'region','Edge',[8,9,15,16]);
Nodes_onBoundary4 = findNodes(mesh,'region','Edge',[14,15,21,22]);
Nodes_Dirichlet1 = [Nodes_onBoundary1,Nodes_onObj_Circle];
Nodes_Dirichlet2 = [Nodes_onBoundary2,Nodes_onObj_Ellipse];
Nodes_Dirichlet3 = [Nodes_onBoundary3,Nodes_onObj_Square];
Nodes_Dirichlet4 = [Nodes_onBoundary4,Nodes_onObj_Tri];
Nodes_DirichletSubDomain_List = {Nodes_Dirichlet1,Nodes_Dirichlet2,Nodes_Dirichlet3,Nodes_Dirichlet4};
LocalDirichletNodesIdx_List = {0,0,0,0}; % the local index of the interior nodes
for this_dm=1:4
        LocalDirichletNodesIdx = zeros(size(Nodes_DirichletSubDomain_List{this_dm}));
        for i=1:length(Nodes_DirichletSubDomain_List{this_dm})
            LocalDirichletNodesIdx(i) = find(Nodes_SubDomain_List{this_dm}==Nodes_DirichletSubDomain_List{this_dm}(i)); % the local node index is the index of the global index array 
        end
        LocalDirichletNodesIdx_List{this_dm}=LocalDirichletNodesIdx;
end

%% file
for rank=1:4
    Num_Nodes=size(Nodes_SubDomain_List{rank},2);
    Num_Elements=size(Elements_List{rank},2);
    NodeCoord=GlobalNodeCoord(Nodes_SubDomain_List{rank},:);
    Element2NodeList=LocalElement2LocalNode_List{rank};
    Node_Source=LocalSourceNodesIdx_List{rank};
    Elements_Node_Source=LocalElementsSourceNodesIdx_List{rank};
    Nodes_InteriorSubDomain=LocalInteriorNodesIdx_List{rank};
    Nodes_InterfaceSubDomain=LocalInterefaceNodesIdx_List{rank};
    Nodes_Dirichlet = LocalDirichletNodesIdx_List{rank};
    
    f=fopen(['naive_fe_ddm_rank',num2str(rank),'_',num2str(Hmax),'.txt'],'w');
    fprintf(f,'%d \n',size(mesh.Nodes,2)); % total nodes
    fprintf(f,'%d %d\n',[Num_Nodes,Num_Elements]);
    
    for i=1:Num_Nodes
        fprintf(f,'%f %f ',NodeCoord(i,:));
    end
    fprintf(f,'\n');
    for i=1:Num_Elements
        fprintf(f,'%d %d %d ',Element2NodeList(i,:)-1); % should be local
    end
    fprintf(f,'\n');
    
    for i=1:length(Nodes_InteriorSubDomain)
        fprintf(f,'%d ',Nodes_InteriorSubDomain(i)-1); % should be local
    end
    fprintf(f,'\n');
    
    for i=1:length(Nodes_InteriorSubDomain_List{rank})
        fprintf(f,'%d ',Nodes_InteriorSubDomain_List{rank}(i)-1);% should be global
    end
    fprintf(f,'\n');
    % x3, in which other 3 ranks
    for ii=1:3
        for i=1:length(Nodes_InterfaceSubDomain{ii})
            fprintf(f,'%d ',Nodes_InterfaceSubDomain{ii}(i)-1); % should be local
        end
        fprintf(f,'\n');
    end

    for i=1:length(Node_Source)
        fprintf(f,'%d ',Node_Source(i)-1); % should be local
    end
    fprintf(f,'\n');

    for i=1:length(Elements_Node_Source)
        fprintf(f,'%d ',Elements_Node_Source(i)-1); % should be local
    end
    fprintf(f,'\n');

    for i=1:length(Nodes_Dirichlet)
            fprintf(f,'%d ',Nodes_Dirichlet(i)-1);
    end
    fprintf(f,'\n');
    
    for ii=1:4 % to the other three
        if ii==rank
            continue;
        end
        for i=1:length(ThisSubABCNodes2ThatSubLocalNodes{ii,rank})
            fprintf(f,'%d ',ThisSubABCNodes2ThatSubLocalNodes{ii,rank}(i)-1); % which nodes to send, should be local
        end
        fprintf(f,'\n');
    end
    fclose(f);
end
fclose('all');