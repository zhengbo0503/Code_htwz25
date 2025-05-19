clc; clear; close all;

% Load data from Julia
data_order = readmatrix("../../Julia/example/timing_order.csv");
data_kappa = readmatrix("../../Julia/example/timing_cnd.csv");

% Separate data_order
N = data_order(:,1);
order_tj = data_order(:,2);
order_tm2 = data_order(:,3);
order_tm3 = data_order(:,4);
order_tm2Prec = data_order(:,5);
order_tm2Apply = data_order(:,6);
order_tm2Jacobi = data_order(:,7);
order_tm2Else = data_order(:,8);
order_tm3Prec = data_order(:,9);
order_tm3Apply = data_order(:,10);
order_tm3Jacobi = data_order(:,11);
order_tm3Else = data_order(:,12);

% Separate data_kappa 
kappa = data_kappa(:,1);
kappa_tj = data_kappa(:,2);
kappa_tm2 = data_kappa(:,3);
kappa_tm3 = data_kappa(:,4);
kappa_tm2Prec = data_kappa(:,5);
kappa_tm2Apply = data_kappa(:,6);
kappa_tm2Jacobi = data_kappa(:,7);
kappa_tm2Else = data_kappa(:,8);
kappa_tm3Prec = data_kappa(:,9);
kappa_tm3Apply = data_kappa(:,10);
kappa_tm3Jacobi = data_kappa(:,11);
kappa_tm3Else = data_kappa(:,12);

%% Figure 1 for grand timings 
close all; figure(1); 
subplot(1,2,1);
plot(N, order_tj,"-o");
hold on;
plot(N, order_tm2, "-x");
plot(N, order_tm3, "->");
xlim([100,1e3]);
xlabel("$n$", Interpreter="latex")
xticks([100,300,500,700,1000]);
axis square
grid on

subplot(1,2,2);
semilogx(kappa, kappa_tj, "-o");
hold on;
semilogx(kappa, kappa_tm2, "-x");
semilogx(kappa, kappa_tm3, "->");
xlim([1e3,1e14]);
xlabel("$\kappa_2(A)$", Interpreter="latex")
xticks([1e3,1e6,1e9,1e12,1e14])
axis square
grid on

legend("Jacobi", "MP2Jacobi", "MP3Jacobi", Location="north", NumColumns=3);

export_fig("./output/timing_union.pdf");

%% Figure 2 for Decomposing the timing for Orders

% Specical colors such that one can recognize the difference in grayscale
segmentColors=[0,0,255;255,0,0;0,255,0;255,255,0]./255;

close all; figure(2)
subplot(1,2,1); 
plt1=bar([order_tm2Prec, order_tm2Apply, order_tm2Jacobi, order_tm3Else] , "stacked");
for k=1:numel(plt1)
    plt1(k).FaceColor = segmentColors(k,:); 
end
xlim([1 21]);
xticks([1,4,8,12,16,20]);
xticklabels(N([1,4,8,12,16,20]))
xlabel("$n$", Interpreter="latex")
title("MP2Jacobi", Interpreter="latex")
axis square 
grid on

subplot(1,2,2); 
plt2=bar([order_tm3Prec, order_tm3Apply, order_tm3Jacobi, order_tm3Else] , "stacked");
for k=1:numel(plt2)
    plt2(k).FaceColor = segmentColors(k,:); 
end
xlim([1 21]);
xticks([1,4,8,12,16,20]);
xticklabels(N([1,4,8,12,16,20]))
xlabel("$n$", Interpreter="latex")
title("MP3Jacobi", Interpreter="latex")
axis square 
grid on

% Add grand legend
legend("Form preconditioner", "Apply preconditioner", "Jacobi", "Everything else", ...
    Location="northwest", NumColumns=4, Direction='normal')

export_fig("./output/timing_decompose.pdf");