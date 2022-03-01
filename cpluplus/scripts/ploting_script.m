%% serial
% close all;
% unkn = [4580,5116,6094,7191,8294,10179,12370,15587,20277];
% time = [8.36184 ,10.2035,16.1075 ,24.8173 ,35.4795 ,62.3584,100.883 ,181.851 ,330.56  ];
% 
% figure;hold on;
% plot(unkn,time,'-b','Linewidth',4)
% 
% xlabel('Number of Unknowns')
% ylabel('Simulation Time (seconds)')
% title('Time Consumption')
% %set(gca, 'YScale', 'log')
% set(gca,'fontsize',24);grid on;
% hold off;
%% parallel result
num_proc = [4,9,16,25,36];
Serial_Result_Time = zeros(1,24);
filename = 'benchmark_serial.txt';
fileID = fopen(filename);
for unknIdx=1:length(Serial_Result_Time)
    if unknIdx==1
        textscan(fileID , 'Testing Serial', 1);
    end
    textscan(fileID , 'Doing Hmax = %d ', 1);
    C=textscan(fileID , ' ... ASM serial version costs time %f s ... ', 1);
    Serial_Result_Time(unknIdx) = C{1};
end
fclose(fileID);

Parallel_Result_Unknown = zeros(length(num_proc),length(Serial_Result_Time));
Parallel_Result_Time = zeros(length(num_proc),length(Serial_Result_Time));
filename = 'benchmark_mpi_openmp.txt';
fileID = fopen(filename);
for unknIdx=1:length(Serial_Result_Time)
    if unknIdx==1
        textscan(fileID , 'BEGIN _VARIES_', 1);
    end
    textscan(fileID , 'Doing Hmax = %d ', 1);
    for procIdx=1:length(num_proc)
        textscan(fileID , 'Doing MPI CPUS = %d ', 1);
        C=textscan(fileID , 'ASM end ... start gathering ... Total_Num_Nodes %f ...end gathering ... ASM MPI end ...  ... MPI ASM FEM costs time %f s ... ', 1);
        Parallel_Result_Unknown(procIdx,unknIdx) = C{1};
        Parallel_Result_Time(procIdx,unknIdx) = C{2};
    end
end
fclose(fileID);
%% Iso efficiency
% X=repmat(num_proc,[length(Serial_Result_Time),1]);
% ProblemSize=Parallel_Result_Unknown';
% Y=ProblemSize';
% Z=Parallel_Result_Time'/
% 
% n = 100; % number of contour levels
% efficiency_map = contourf(X,Y,Z);
% colorbar;
% xlabel('{\it\Omega}_r [rev/min]')
% ylabel('{\itT}_e [N*m]')
% legend
%% serial vs parallel
figure(1);hold on;
plot(Parallel_Result_Unknown(1,:),Parallel_Result_Time(1,:),'-r','Linewidth',4)
plot(Parallel_Result_Unknown(1,:),Serial_Result_Time(1,:),'-b','Linewidth',4)
legend('Parallel Result','Serial Result')
xlabel('Unknowns')
ylabel('Time (s)')
% set(gca, 'YScale', 'log')
set(gca,'fontsize',24);grid on;
hold off;

%% parallel 

figure(2);hold on;
plot(Parallel_Result_Unknown(1,:),Parallel_Result_Time(1,:),'-r','Linewidth',4)
plot(Parallel_Result_Unknown(2,:),Parallel_Result_Time(2,:),'-g','Linewidth',4)
plot(Parallel_Result_Unknown(3,:),Parallel_Result_Time(3,:),'-b','Linewidth',4)
plot(Parallel_Result_Unknown(4,:),Parallel_Result_Time(4,:),'-k','Linewidth',4)
plot(Parallel_Result_Unknown(5,:),Parallel_Result_Time(5,:),'-m','Linewidth',4)
legend('4 Ranks','9 Ranks','16 Ranks','25 Ranks','36 Ranks')
 legend('Location','best')
xlabel('Unknowns')
ylabel('Time (s)')
% set(gca, 'YScale', 'log')
set(gca,'fontsize',24);grid on;
hold off;

%%

N_nodeSize = 31;
Num_Unknowns=zeros(1,N_nodeSize);
FourDM_Parallel_Result_Time=zeros(1,N_nodeSize);
FourDM_Serial_Result_Time=zeros(1,N_nodeSize);

filename = 'benchmark_serial_more.txt';
fileID = fopen(filename);
for unknIdx=1:N_nodeSize
    if unknIdx==1
        textscan(fileID , 'Testing Serial', 1);
    end
    textscan(fileID , 'Doing Hmax = %d ', 1);
    C=textscan(fileID , ' ... ASM serial version costs time %f s ... ', 1);
    FourDM_Serial_Result_Time(unknIdx) = C{1};
end
fclose(fileID);

filename = 'benchmark_mpi_openmp_4dm.txt';
fileID = fopen(filename);
for unknIdx=1:N_nodeSize
    if unknIdx==1
        textscan(fileID , 'BEGIN _VARIES_', 1);
    end
    textscan(fileID , 'Doing MPI CPUS = %d ', 1);
    textscan(fileID , 'Doing Hmax = %d ', 1);
    C=textscan(fileID , 'ASM end ... start gathering ... Total_Num_Nodes %f ...end gathering ... ASM MPI end ...  ... MPI ASM FEM costs time %f s ... ', 1);
    Num_Unknowns(unknIdx) = C{1};
    FourDM_Parallel_Result_Time(unknIdx) = C{2};
end
fclose(fileID);
figure(3);hold on;
plot(Num_Unknowns,FourDM_Serial_Result_Time,'-b','Linewidth',4)
plot(Num_Unknowns,FourDM_Parallel_Result_Time,'-r','Linewidth',4)
legend('Serial','Parallel')
legend('Location','best')
xlabel('Unknowns')
ylabel('Time (s)')
% set(gca, 'YScale', 'log')
set(gca,'fontsize',24);grid on;
hold off;
%%
figure(4);hold on;
yyaxis left;
plot(Num_Unknowns,FourDM_Serial_Result_Time./FourDM_Parallel_Result_Time,'-k','Linewidth',4)
xlabel('Unknowns')
ylabel('Speedup')
yyaxis right;
plot(Num_Unknowns,FourDM_Serial_Result_Time./FourDM_Parallel_Result_Time/4,'-k','Linewidth',4)
ylabel('Parallel Efficiency')
set(gca,'fontsize',24);grid on;
hold off;
