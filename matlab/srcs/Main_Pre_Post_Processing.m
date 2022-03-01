%% settings
clear all; close all;
N = 100; % interpolation size
Hmax=0.5; % element size
use_parallel_for = false; % if use parfor to parallel the ASM procedure, might be slow due to overhead
plot_logistics = true; % if to plot FEM matrices and meshes
%% generate geometry & create mesh using meshpde toolbox 
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
stlwrite('naive_fe_ddm.stl', mesh.Elements', [mesh.Nodes',zeros(Num_Nodes,1)]);
disp(['Mesh generated and save to',' naive_fe_ddm.stl'])
%% benchmarking
X=linspace(-5,5,N);
Y=linspace(-5,5,N);
[Phi1,timeVal_CGM,error_list_CGM] = C1_FEM_CGM(mesh,N,plot_logistics);
[Phi2,timeVal_ASM_DDM,error_list_ASM_DDM] = C2_FEM_DDM_ASM(mesh,N,use_parallel_for,plot_logistics);
[Phi3,timeVal_PCG,error_list_PCG] = C3_FEM_DDM_PCG(mesh,N,plot_logistics); 
clc;
%% interpolated results
figure; hold on;
title('\Phi_1');
h1 = imagesc(X,Y(end:-1:1),rot90(Phi1));
pdemesh(mesh.Nodes,mesh.Elements);
set(gca,'fontsize',24);
axis tight;
axis equal; hold off;

figure; hold on;
title('\Phi_2');
h2 = imagesc(X,Y(end:-1:1),rot90(Phi2));
pdemesh(mesh.Nodes,mesh.Elements);
set(gca,'fontsize',24);
axis tight;
axis equal; hold off;

figure; hold on;
title('\Phi_3');
h3 = imagesc(X,Y(end:-1:1),rot90(Phi3));
pdemesh(mesh.Nodes,mesh.Elements);
set(gca,'fontsize',24);
axis tight;
axis equal; hold off;


diff1=abs(Phi1-Phi2);
diff2=abs(Phi1-Phi3);
diff3=abs(Phi2-Phi3);
figure;imagesc(X,Y(end:-1:1),rot90(diff1)); % FEM_CGM vs FEM_ASM_DDM
title('Abs(\Phi_1-\Phi_2)');
set(gca,'fontsize',24);
axis equal;axis tight;
figure;imagesc(X,Y(end:-1:1),rot90(diff2)); % FEM_CGM vs FEM_ASM_PCG
title('Abs(\Phi_1-\Phi_3)');
set(gca,'fontsize',24);
axis equal;axis tight;
figure;imagesc(X,Y(end:-1:1),rot90(diff3)); % FEM_ASM_DDM vs FEM_ASM_PCG
title('Abs(\Phi_2-\Phi_3)');
set(gca,'fontsize',24);
axis equal;axis tight;

%% iterative method convergence
figure;hold on;
plot(error_list_CGM,'-g','Linewidth',4)
plot(error_list_ASM_DDM,'-b','Linewidth',4)
plot(error_list_PCG,'-m','Linewidth',4)

legend('FEM-CG','FE-DDM-ASM','FE-DDM-PCG')
xlabel('iteration')
ylabel('error (residual)')
set(gca, 'YScale', 'log')
set(gca,'fontsize',24);grid on;
hold off;

%% timing
timeVal_CGM
timeVal_ASM_DDM
timeVal_PCG