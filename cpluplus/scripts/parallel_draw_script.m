%% define geometry
% you only need to modify the R (grid size) and Hmax (mesh quality)
clear all; close all;
R=2; % R x R grid
Hmax=0.5; % element size

H_List = [1.5,1.2,0.9,0.6,0.3,0.1];
R_List = [2,3,4,5,6];
for Hmax=H_List
for R=R_List % R x R grid
    
N = 50*R; % interpolation size
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
figure(1)
pdegplot(gme,'FaceLabels','on','EdgeLabels','on')
figure(2)
pdemesh(mesh,'ElementLabels','off')


%%
Num_Nodes = size(mesh.Nodes,2);
Num_Elements = size(mesh.Elements,2);
GlobalElement2NodeList = mesh.Elements';
GlobalNodeCoord = mesh.Nodes';

x_pos_mid=(0:R*2)*Domain_length/2;
x_pos_back=-overlap+(0:R-1)*Domain_length;
x_pos_for=+overlap+(1:R)*Domain_length;
[GridCol,GridRow] = meshgrid(x_pos_mid,x_pos_mid);
% Find the nodes and elements associated with Row,Col (Row,Col)
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

%%
file_read = ['naive_fe_ddm',num2str(R^2),'_','out_',num2str(Hmax),'.txt'];
first_line = dlmread(file_read, ' ');
phi = zeros(Num_Nodes,1);
for i=1:length(first_line)
    phi(i)=first_line(i);
end
%%
Nel = @(x,y,ael,bel,cel,Delta_e) 1/(2*Delta_e)*(ael+bel*x+cel*y);
X=linspace(-overlap,Domain_length*R+overlap,N);
Y=linspace(-overlap,Domain_length*R+overlap,N);

Phi = zeros(N,N);
Num_Elements = size(mesh.Elements,2);
for i=1:Num_Elements
    if mod(i,100)==0
        disp(['Interpolating ',num2str(i),'th element'])
    end
    
    who_am_I = mesh.Elements(1,i); % the global node idx of node 1 in Elements(i)
    where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 1 in Elements(i)
    xe1 = where_am_I(1);
    ye1 = where_am_I(2);
    ne1 = who_am_I;
    
    who_am_I = mesh.Elements(2,i); % the global node idx of node 1 in Elements(i)
    where_am_I = mesh.Nodes(:,who_am_I); % the global coordinate of node 1 in Elements(i)
    xe2 = where_am_I(1);
    ye2 = where_am_I(2);
    ne2 = who_am_I;
    
    who_am_I = mesh.Elements(3,i); % the global node idx of node 1 in Elements(i)
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
    for q=1:N
        for p=1:N
            [in,on]=inpolygon(X(q),Y(p),[xe1,xe2,xe3],[ye1,ye2,ye3]);

            if (in || on) %indicating if inside the element or on the edge of the element
                Phi(q,p)=Nel(X(q),Y(p),ae1,be1,ce1,Delta_e)*phi(ne1)+Nel(X(q),Y(p),ae2,be2,ce2,Delta_e)*phi(ne2)+Nel(X(q),Y(p),ae3,be3,ce3,Delta_e)*phi(ne3);
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
axis equal; 
saveas(gcf,['Phi_R',num2str(R),'H',num2str(length(Nodes_SubDomain{1,1})),'.png']); 
hold off;
end
end