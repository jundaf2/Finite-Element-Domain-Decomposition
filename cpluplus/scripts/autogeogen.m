clear all; close all;
%% define geometry
R_List=[2];%[2,3,4,5,6];
H_List=0.095:-0.005:0.03;%H_List=0.19:-0.01:0.11; % 1.5:-0.05:0.1
for R=R_List % R x R grid
    for Hmax=H_List
% R=6; 
N = 100; % interpolation size
% Hmax=1; % element size
model = createpde;

Domain_length = 5;
overlap = 0.5;
x_starts = Domain_length*(0:R-1)-overlap;
x_ends = Domain_length*(1:R)+overlap;
y_starts = Domain_length*(0:R-1)-overlap;
y_ends = Domain_length*(1:R)+overlap;

obj_center_x = Domain_length*(0:R-1)+Domain_length/2;
obj_center_y = Domain_length*(0:R-1)+Domain_length/2;
obj_radius = 1;

gm=[];
sf=[];
ns=[];
for Row=1:R
    for Col=1:R
        domainID = (Row-1)*R+Col;
        if(domainID>=10)
            dmstr=num2str(domainID);
        else
            dmstr=['0',num2str(domainID)];
        end
        eval(['Rect',dmstr,'=[3,4,x_starts(Col),x_ends(Col),x_ends(Col),x_starts(Col),y_starts(Row),y_starts(Row),y_ends(Row),y_ends(Row)]''']);
        eval(['Circle',dmstr,'=[1,obj_center_x(Col),obj_center_y(Row),obj_radius]''']);
        eval(['Circle',dmstr,'=[eval([''Circle'',dmstr]);zeros(length(eval([''Rect'',dmstr])) - length(eval([''Circle'',dmstr])),1)]']);
        eval(['gm=','[gm,eval([''Rect'',dmstr]),eval([''Circle'',dmstr])]']);
        if(domainID==1)
            sf = 'Rect01-Circle01';
        else
            eval(['sf=','[sf,','''+Rect',dmstr,'-Circle',dmstr,''']']);
        end
        ns = [ns;char(['Rect',dmstr],['Circle',dmstr])];
    end
end

%% plot
ns = ns';
g = decsg(gm,sf,ns);
gme = geometryFromEdges(model,g);
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',Hmax);
% figure(1)
% pdegplot(gme,'FaceLabels','on','EdgeLabels','on')

%% generate information
Num_Nodes = size(mesh.Nodes,2);
Num_Elements = size(mesh.Elements,2);
GlobalElement2NodeList = mesh.Elements';
GlobalNodeCoord = mesh.Nodes';

x_pos_mid=(0:R*2)*Domain_length/2;
x_pos_back=-overlap+(0:R-1)*Domain_length;
x_pos_for=+overlap+(1:R)*Domain_length;
[GridCol,GridRow] = meshgrid(x_pos_mid,x_pos_mid);
% Find the nodes and elements associated with subdomain (Row,Col)
Nodes_SubDomain=cell(R,R);
Elements_SubDomain=cell(R,R);
for Row=1:R
    for Col=1:R
        xypos = [];
        for colcol=1:3
            for rowrow=1:3
                xypos=[xypos;[GridCol((Row-1)*2+rowrow,(Col-1)*2+colcol),GridRow((Row-1)*2+rowrow,(Col-1)*2+colcol)]];
            end
        end
        faceIDs = nearestFace(gme,xypos);
        Nodes_SubDomain{Row,Col}=unique(findNodes(mesh,'region','Face',faceIDs));
        Elements_SubDomain{Row,Col}=unique(findElements(mesh,'region','Face',faceIDs));
    end
end
figure;
pdemesh(mesh,'ElementLabels','off')
saveas(gcf,['R',num2str(R),'H',num2str(length(Nodes_SubDomain{1,1})),'.png']); 
% Objs
Nodes_Object=cell(R,R);
for Row=1:R
    for Col=1:R
        xypos = [];
        for colcol=-1:2:1  % -1 0 1 does not matter
            for rowrow=-1:2:1 
                xypos=[xypos;[x_pos_mid(Col*2)+colcol*obj_radius/2,x_pos_mid(2*Row)+rowrow*obj_radius/2]];
            end
        end
        edgeIDs = nearestEdge(gme,xypos);
        Nodes_Object{Row,Col}=findNodes(mesh,'region','Edge',edgeIDs);
    end
end
% Outer Boundary
Nodes_Boundary=cell(R,R);
for Row=1:R
    for Col=1:R
        xypos = [];
        if Row==R
            xypos=[xypos;[x_pos_mid(Col*2-1),x_pos_for(Row)]];
            xypos=[xypos; [x_pos_mid(Col*2),x_pos_for(Row)]];
            xypos=[xypos;[x_pos_mid(Col*2+1),x_pos_for(Row)]];
        end
        if Col==R
            xypos=[xypos;[x_pos_for(Col),x_pos_mid(2*Row+1)]];
            xypos=[xypos; [x_pos_for(Col),x_pos_mid(2*Row)]];
            xypos=[xypos;[x_pos_for(Col),x_pos_mid(2*Row-1)]];
        end
        if Row==1
            xypos=[xypos;[x_pos_mid(2*Col+1),x_pos_back(Row)]];
            xypos=[xypos; [x_pos_mid(Col*2),x_pos_back(Row)]];
            xypos=[xypos;[x_pos_mid(Col*2-1),x_pos_back(Row)]];
        end
        if Col==1
            xypos=[xypos;[x_pos_back(Col),x_pos_mid(2*Row-1)]];
            xypos=[xypos; [x_pos_back(Col),x_pos_mid(2*Row)]];
            xypos=[xypos;[x_pos_back(Col),x_pos_mid(2*Row+1)]];
        end
        edgeIDs = [];
        if ~isempty(xypos)
            edgeIDs = nearestEdge(gme,xypos);
        end
        
        Nodes_Boundary{Row,Col}=findNodes(mesh,'region','Edge',edgeIDs);
    end
end
% Dirichlet
Nodes_Dirichlet=cell(R,R);
for Row=1:R
    for Col=1:R
        Nodes_Dirichlet{Row,Col}=[Nodes_Object{Row,Col},Nodes_Boundary{Row,Col}];
    end
end
% Sources
Nodes_Source=cell(R,R); 
for Row=1:R
    for Col=1:R
        xypos = [];
        for colcol=-1:2:1
            for rowrow=-1:2:1 
                xypos=[xypos,[x_pos_mid(Col*2)+colcol*3/2;x_pos_mid(2*Row)+rowrow*3/2]];
            end
        end
        Nodes_Source{Row,Col}=findNodes(mesh,'nearest',xypos);
    end
end

% 1  2  3  4  5  6  7  8 
% NW N  NE E  SE S  SW W
Nodes_Interface = cell(R,R,R^2-1);
Nodes_Interior = Nodes_SubDomain;
for Row=1:R
    for Col=1:R
        xypos = [];
        xypos=[xypos;[x_pos_back(Col),x_pos_mid(2*Row+1)]];
        xypos=[xypos;[x_pos_mid(Col*2-1),x_pos_for(Row)]];
        edgeIDs = nearestEdge(gme,xypos);
        Nodes_Interface{Row,Col,1}=findNodes(mesh,'region','Edge',edgeIDs);
        
        xypos = [];
        xypos=[xypos; [x_pos_mid(Col*2),x_pos_for(Row)]];
        edgeIDs = nearestEdge(gme,xypos);
        Nodes_Interface{Row,Col,2}=findNodes(mesh,'region','Edge',edgeIDs);
        
        xypos = [];
        xypos=[xypos;[x_pos_mid(Col*2+1),x_pos_for(Row)]];
        xypos=[xypos;[x_pos_for(Col),x_pos_mid(2*Row+1)]];
        edgeIDs = nearestEdge(gme,xypos);
        Nodes_Interface{Row,Col,3}=findNodes(mesh,'region','Edge',edgeIDs);
        
        xypos = [];
        xypos=[xypos; [x_pos_for(Col),x_pos_mid(2*Row)]];
        edgeIDs = nearestEdge(gme,xypos);
        Nodes_Interface{Row,Col,4}=findNodes(mesh,'region','Edge',edgeIDs);

        xypos = [];
        xypos=[xypos;[x_pos_for(Col),x_pos_mid(2*Row-1)]];
        xypos=[xypos;[x_pos_mid(2*Col+1),x_pos_back(Row)]];
        edgeIDs = nearestEdge(gme,xypos);
        Nodes_Interface{Row,Col,5}=findNodes(mesh,'region','Edge',edgeIDs);

        xypos = [];
        xypos=[xypos; [x_pos_mid(Col*2),x_pos_back(Row)]];
        edgeIDs = nearestEdge(gme,xypos);
        Nodes_Interface{Row,Col,6}=findNodes(mesh,'region','Edge',edgeIDs); 

        xypos = [];
        xypos=[xypos;[x_pos_mid(Col*2-1), x_pos_back(Row)]];
        xypos=[xypos;[x_pos_back(Col),x_pos_mid(2*Row-1)]];
        edgeIDs = nearestEdge(gme,xypos);
        Nodes_Interface{Row,Col,7}=findNodes(mesh,'region','Edge',edgeIDs);

        xypos = [];
        xypos=[xypos; [x_pos_back(Col),x_pos_mid(2*Row)]];
        edgeIDs = nearestEdge(gme,xypos);
        Nodes_Interface{Row,Col,8}=findNodes(mesh,'region','Edge',edgeIDs); 
        
        i=1;
        for j=[2,8]
            Nodes_Interface{Row,Col,i}=setdiff(Nodes_Interface{Row,Col,i},Nodes_Interface{Row,Col,j});
        end
        
        i=3;
        for j=[2,4]
            Nodes_Interface{Row,Col,i}=setdiff(Nodes_Interface{Row,Col,i},Nodes_Interface{Row,Col,j});
        end
        i=5;
        for j=[4,6]
            Nodes_Interface{Row,Col,i}=setdiff(Nodes_Interface{Row,Col,i},Nodes_Interface{Row,Col,j});
        end
        i=7;
        for j=[6,8]
            Nodes_Interface{Row,Col,i}=setdiff(Nodes_Interface{Row,Col,i},Nodes_Interface{Row,Col,j});
        end
        
        for i=1:8
            Nodes_Interface{Row,Col,i}=setdiff(Nodes_Interface{Row,Col,i},Nodes_Dirichlet{Row,Col});
        end
        for i=1:8
            Nodes_Interior{Row,Col}=setdiff(Nodes_Interior{Row,Col},Nodes_Interface{Row,Col,i});
        end
    end
end

%% to local
OneDRank = @(Row,Col) (Row-1)*R+Col;
LocalInteriorNodesIdx_List = cell(R^2,1); % the local index of the interior nodes
for Row=1:R
    for Col=1:R
        LocalInteriorNodesIdx = zeros(size(Nodes_Interior{Row,Col}));
        for i=1:length(Nodes_Interior{Row,Col})%length(Nodes_InteriorSubDomain_List{ABC_InWhichRank})
            LocalInteriorNodesIdx(i) = find(Nodes_SubDomain{Row,Col}== Nodes_Interior{Row,Col}(i)); % the local node index is the index of the global index array 
        end
        LocalInteriorNodesIdx_List{OneDRank(Row,Col)}=LocalInteriorNodesIdx;
    end
end

% for assembly
% search for those local nodes that have nodes on AB put them at the front
% and put other ranks that do not have at the back
% LocalInterefaceNodesIdx_List = cell(R^2,R^2); % the local index of the interface nodes
% for Row=1:R
%     for Col=1:R
%         count=0;
%         for those_dm=1:R^2
%             % those_dm: 1-8 nearby clock-wise start from up-left  but does not necessarily have nodes (isempty)
%             %           9-15 non overlap with Nodes_Interface isempty
%             if ~isempty(Nodes_Interface{Row,Col,those_dm})
%                 count=count+1;
%                 LocalInterefaceNodesIdx = zeros(size(Nodes_Interface{Row,Col,those_dm}));
%                 for i=1:length(Nodes_Interface{Row,Col,those_dm})%length(Nodes_InteriorSubDomain_List{ABC_InWhichRank})
%                     %disp([num2str(Row),' ',num2str(Col),' ',num2str(those_dm),' ',num2str(i)]);
%                     LocalInterefaceNodesIdx(i) = find(Nodes_SubDomain{Row,Col}==Nodes_Interface{Row,Col,those_dm}(i)); % the local node index is the index of the global index array 
%                 end
%                 LocalInterefaceNodesIdx_List{OneDRank(Row,Col),count}=LocalInterefaceNodesIdx; 
%             end    
%         end
%         for i=(count+1):R^2
%             LocalInterefaceNodesIdx_List{OneDRank(Row,Col),i} = 0; % do not have interface to that rank
%         end
%             
%     end
% end

% for send
Nodes_Interface_InWhichRank=zeros(R,R,R^2); % To which rank
Nodes_Interface_New = cell(R,R,R^2);
for Row=1:R
    for Col=1:R
        count=1;
        mark=zeros(8,1);
        temp_mark=0;
        % keep the sequence of 1-8
        if Row~=R && Col~=1
            %count=count+1;
            %Nodes_Interface_InWhichRank(Row,Col,1)=OneDRank(Row+1,Col-1);
            if ~isempty(Nodes_Interface{Row,Col,1}) && ~mark(1)
                mark(1)=1;
                Nodes_Interface_InWhichRank(Row,Col,count)=OneDRank(Row+1,Col-1);
                Nodes_Interface_New{Row,Col,count}=[Nodes_Interface_New{Row,Col,count},Nodes_Interface{Row,Col,1}];
                temp_mark=1;
            end
            if temp_mark==1
                count=count+1;
                temp_mark=0;
            end
        end
        if Row~=R
            %count=count+1;
            %Nodes_Interface_InWhichRank(Row,Col,2)=OneDRank(Row+1,Col);
            for i=[1,2]
                if ~isempty(Nodes_Interface{Row,Col,i}) && ~mark(i)
                    mark(i)=1;
                    Nodes_Interface_InWhichRank(Row,Col,count)=OneDRank(Row+1,Col);
                    Nodes_Interface_New{Row,Col,count}=[Nodes_Interface_New{Row,Col,count},Nodes_Interface{Row,Col,i}];
                    temp_mark=1;
                end
            end
            if temp_mark==1
                count=count+1;
                temp_mark=0;
            end
        end
        
        if Row~=R && Col~=R
            if ~isempty(Nodes_Interface{Row,Col,3}) && ~mark(3)
                mark(3)=1;
                Nodes_Interface_InWhichRank(Row,Col,count)=OneDRank(Row+1,Col+1);
                Nodes_Interface_New{Row,Col,count}=[Nodes_Interface_New{Row,Col,count},Nodes_Interface{Row,Col,3}];
                temp_mark=1;
            end
            if temp_mark==1
                count=count+1;
                temp_mark=0;
            end
        end
        
        if Col~=R   
            %count=count+1;
            %Nodes_Interface_InWhichRank(Row,Col,4)=OneDRank(Row,Col+1);
            for i=[3,4]
                if ~isempty(Nodes_Interface{Row,Col,i}) && ~mark(i)
                    mark(i)=1;
                    Nodes_Interface_InWhichRank(Row,Col,count)=OneDRank(Row,Col+1);
                    Nodes_Interface_New{Row,Col,count}=[Nodes_Interface_New{Row,Col,count},Nodes_Interface{Row,Col,i}];
                    temp_mark=1;
                end
            end
            if temp_mark==1
                count=count+1;
                temp_mark=0;
            end
        end
        
        if Row~=1 && Col~=R
            %count=count+1;
            %Nodes_Interface_InWhichRank(Row,Col,5)=OneDRank(Row-1,Col+1);
            if ~isempty(Nodes_Interface{Row,Col,5}) && ~mark(5)
                mark(5)=1;
                temp_mark=1;
                Nodes_Interface_InWhichRank(Row,Col,count)=OneDRank(Row-1,Col+1);
                Nodes_Interface_New{Row,Col,count}=[Nodes_Interface_New{Row,Col,count},Nodes_Interface{Row,Col,5}];
            end
            if temp_mark==1
                count=count+1;
                temp_mark=0;
            end
        end
        if Row~=1   
            %count=count+1;
            %Nodes_Interface_InWhichRank(Row,Col,6)=OneDRank(Row-1,Col);
            for i=[5,6]
                if ~isempty(Nodes_Interface{Row,Col,i}) && ~mark(i)
                    mark(i)=1;
                    Nodes_Interface_InWhichRank(Row,Col,count)=OneDRank(Row-1,Col);
                    Nodes_Interface_New{Row,Col,count}=[Nodes_Interface_New{Row,Col,count},Nodes_Interface{Row,Col,i}];
                    temp_mark=1;
                end
            end
            if temp_mark==1
                count=count+1;
                temp_mark=0;
            end
        end
        if Row~=1 && Col~=1
            %count=count+1;
            %Nodes_Interface_InWhichRank(Row,Col,7)=OneDRank(Row-1,Col-1);
            if ~isempty(Nodes_Interface{Row,Col,7}) && ~mark(7)
                mark(7)=1;
                Nodes_Interface_InWhichRank(Row,Col,count)=OneDRank(Row-1,Col-1);
                Nodes_Interface_New{Row,Col,count}=[Nodes_Interface_New{Row,Col,count},Nodes_Interface{Row,Col,7}];
                temp_mark=1;
            end
            if temp_mark==1
                count=count+1;
                temp_mark=0;
            end
        end        
        if Col~=1
            %count=count+1;
            %Nodes_Interface_InWhichRank(Row,Col,8)=OneDRank(Row,Col-1);
            for i=[7,8]
                if ~isempty(Nodes_Interface{Row,Col,i}) && ~mark(i)
                    mark(i)=1;
                    Nodes_Interface_InWhichRank(Row,Col,count)=OneDRank(Row,Col-1);
                    Nodes_Interface_New{Row,Col,count}=[Nodes_Interface_New{Row,Col,count},Nodes_Interface{Row,Col,i}];
                    temp_mark=1;
                end
            end
            if temp_mark==1
                count=count+1;
                temp_mark=0;
            end
        end
        
        
        % other are all zero ,so there sequence does not matter
        all_others=1:R^2;
        others = ismember(all_others,[squeeze(Nodes_Interface_InWhichRank(Row,Col,1:(count-1)))']);
        Nodes_Interface_InWhichRank(Row,Col,count:end)=all_others(~others);
    end
end

% for recv
LocalInterefaceNodesIdx_List = cell(R^2,R^2); % the local index of the interface nodes
for Row=1:R
    for Col=1:R
        count=1;
        for those_dm=1:R^2
            % those_dm: 1-8 nearby clock-wise start from up-left  but does not necessarily have nodes (isempty)
            %           9-15 non overlap with Nodes_Interface isempty
            if ~isempty(Nodes_Interface_New{Row,Col,count})
                
                LocalInterefaceNodesIdx = zeros(size(Nodes_Interface_New{Row,Col,count}));
                for i=1:length(Nodes_Interface_New{Row,Col,count})%length(Nodes_InteriorSubDomain_List{ABC_InWhichRank})
                    %disp([num2str(Row),' ',num2str(Col),' ',num2str(those_dm),' ',num2str(i)]);
                    LocalInterefaceNodesIdx(i) = find(Nodes_SubDomain{Row,Col}==Nodes_Interface_New{Row,Col,count}(i)); % the local node index is the index of the global index array 
                end
                LocalInterefaceNodesIdx_List{OneDRank(Row,Col),count}=LocalInterefaceNodesIdx; 
                count=count+1;
            end    
        end
        for i=count:R^2
            LocalInterefaceNodesIdx_List{OneDRank(Row,Col),i} = 0; % do not have interface to that rank
        end
            
    end
end


% the above two: 1:R^2-1, one array tell the local nodes one array tell the rank
% it resides, we expect the following follows the same patern

% the subdomain ABC to local array: store the local indices that resides others boundary
ThisSubABCNodes2ThatSubLocalNodes = cell(R^2,R^2); % R*R subdomains x 8 subdomains
Nodes_SubDomain_T=Nodes_SubDomain';
for Row=1:R
    for Col=1:R
        count = 1;
        for those_dm=1:R^2 % all others, the first ones are those have overlap, the laters are non-overlap but need send empty messages
            ABC_InWhichRank = Nodes_Interface_InWhichRank(Row,Col,count); % overlap in which rank
            %Nodes_Interface_New{Row,Col,those_dm}
            if ~isempty(Nodes_Interface_New{Row,Col,count})
                this2thatLocal = zeros(size(Nodes_Interface_New{Row,Col,count}));
                %temp1=length(Nodes_Interface{Row,Col,those_dm})
                for i=1:length(Nodes_Interface_New{Row,Col,count})%length(Nodes_InteriorSubDomain_List{ABC_InWhichRank})
                    %temp2=Nodes_Interface{Row,Col,those_dm}(i)
                    temp2=find(Nodes_SubDomain_T{ABC_InWhichRank}==Nodes_Interface_New{Row,Col,count}(i)); %--> error is when it is node no 0
                    this2thatLocal(i) = find(Nodes_SubDomain_T{ABC_InWhichRank}==Nodes_Interface_New{Row,Col,count}(i)); % the local node index is the index of the global index array 
                end
                ThisSubABCNodes2ThatSubLocalNodes{OneDRank(Row,Col),ABC_InWhichRank}=this2thatLocal;
                count = count+1;
            end
        end
        for those_dm=1:R^2
            if isempty(ThisSubABCNodes2ThatSubLocalNodes{OneDRank(Row,Col),those_dm})
                ThisSubABCNodes2ThatSubLocalNodes{OneDRank(Row,Col),those_dm}=0;
            end
        end
    end
end

% the local index of the source nodes
LocalSourceNodesIdx_List = cell(R^2,1);
for Row=1:R
    for Col=1:R
        LocalSourceNodesIdx = zeros(size(Nodes_Source{Row,Col}));
        for i=1:length(Nodes_Source{Row,Col})
            LocalSourceNodesIdx(i) = find(Nodes_SubDomain{Row,Col}==Nodes_Source{Row,Col}(i)); % the local node index is the index of the global index array 
        end
        LocalSourceNodesIdx_List{OneDRank(Row,Col)}=LocalSourceNodesIdx;
    end
end

% local element that contains the node with source
LocalElementsSourceNodesIdx_List = cell(R^2,1); % the local index of the element
for Row=1:R
    for Col=1:R
        GlobalElements_Node_Source = findElements(mesh,'attached',Nodes_Source{Row,Col});
        LocalSourceElementsIdx = zeros(size(GlobalElements_Node_Source));
        for i=1:length(GlobalElements_Node_Source)
            LocalSourceElementsIdx(i) = find(Elements_SubDomain{Row,Col}==GlobalElements_Node_Source(i)); % the local node index is the index of the global index array 
        end
        LocalElementsSourceNodesIdx_List{OneDRank(Row,Col)}=LocalSourceElementsIdx;
    end
end

% element2node local
LocalElement2LocalNode_List = cell(R^2,1);
for Row=1:R
    for Col=1:R
        locallist_with_globalIdx = GlobalElement2NodeList(Elements_SubDomain{Row,Col},:);
        LocalElement2LocalNode = zeros(size(locallist_with_globalIdx));
        for i=1:size(locallist_with_globalIdx,2)
            for ii=1:size(locallist_with_globalIdx,1)
                LocalElement2LocalNode(ii,i) = find(Nodes_SubDomain{Row,Col}==locallist_with_globalIdx(ii,i)); % the local node index is the index of the global index array 
            end
        end
        LocalElement2LocalNode_List{OneDRank(Row,Col)}=LocalElement2LocalNode;
    end
end


LocalDirichletNodesIdx_List = cell(R^2,1); % the local index of the interior dirichlet nodes
for Row=1:R
    for Col=1:R
        LocalDirichletNodesIdx = zeros(size(Nodes_Dirichlet{Row,Col}));
        for i=1:length(Nodes_Dirichlet{Row,Col})
            LocalDirichletNodesIdx(i) = find(Nodes_SubDomain{Row,Col}==Nodes_Dirichlet{Row,Col}(i)); % the local node index is the index of the global index array 
        end
        LocalDirichletNodesIdx_List{OneDRank(Row,Col)}=LocalDirichletNodesIdx;
    end
end


%% file
% those have -1 means there is no node for that case
for Row=1:R
    for Col=1:R
        rank = OneDRank(Row,Col);
        Num_Nodes=size(Nodes_SubDomain{Row,Col},2);
        Num_Elements=size(Elements_SubDomain{Row,Col},2);
        NodeCoord=GlobalNodeCoord(Nodes_SubDomain{Row,Col},:);
        Element2NodeList=LocalElement2LocalNode_List{rank};
        Node_Source=LocalSourceNodesIdx_List{rank};
        Elements_Node_Source=LocalElementsSourceNodesIdx_List{rank};
        Nodes_InteriorSubDomain=LocalInteriorNodesIdx_List{rank};
        
        Nodes_Dirichlet = LocalDirichletNodesIdx_List{rank};

        f=fopen(['naive_fe_ddm',num2str(R^2),'_rank',num2str(rank),'_',num2str(Hmax),'.txt'],'w');
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

        for i=1:length(Nodes_Interior{Row,Col})
            fprintf(f,'%d ',Nodes_Interior{Row,Col}(i)-1);% should be global
        end
        fprintf(f,'\n');
        % x15, in which other 15 ranks (at least 7 are empty)
        for ii=1:R^2
            Nodes_InterfaceSubDomain=LocalInterefaceNodesIdx_List{rank,ii};
            for i=1:length(Nodes_InterfaceSubDomain)
                fprintf(f,'%d ',Nodes_InterfaceSubDomain(i)-1); % should be local
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

        for ii=1:R^2 % to the other 15
            for i=1:length(ThisSubABCNodes2ThatSubLocalNodes{ii,rank})
                fprintf(f,'%d ',ThisSubABCNodes2ThatSubLocalNodes{ii,rank}(i)-1); % which nodes to send, should be local
            end
            fprintf(f,'\n');
        end
        
        for ii=1:R^2 % to the other 15 
            fprintf(f,'%d ',ii-1); % which nodes to send, should be local, Nodes_Interface_ToWhichSubDomain
        end
        
        fprintf(f,'\n');
        for i=1:R^2
            fprintf(f,'%d ',Nodes_Interface_InWhichRank(Row,Col,i)-1); % Nodes_Interface_FromWhichSubDomain
        end
        fprintf(f,'\n');
        
        fclose(f);
    end
end
    end
end
fclose('all');