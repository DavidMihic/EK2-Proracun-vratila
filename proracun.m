% Elementi Konstrukcija II MiR
% Program 1 - Proracun Vratila
% David Mihic, 2-mir-2

clear;
close all;
clc;

% Helper functions
function plotGraph(xData, yData, titleText, color, l, l3, l6)
    figure('Name', titleText);
    grid;
    hold on;
    yline(0, 'w');

    xline(l3, 'LineStyle', '--');
    xline(l6, 'LineStyle', '--');
    xline(l, 'LineStyle', '--');

    for i=1:length(xData)
        curX = xData{i};
        curY = yData{i};
        plot(curX, curY(curX), color, 'LineWidth', 2)
        fill([curX, fliplr(curX)], [curY(curX), zeros(1, length(curY(curX)))], color, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
    end

    hold off;
end

% Constants
T = 420000; % Nmm
J2 = 0.0500; % kgm^2 = Nms^2
J3 = 1.050; % kgm^2 = Nms^2
l = 350; % mm
G_z2 = 210; % N
G_z3 = 80; % N
b2 = 110; % mm
b3 = 110; % mm
r2 = 165; % mm
r3 = 57.8; % mm

n = 5.33; % rps
L_h = 8000; % hrs
ALPHA = 20; % deg
ALPHA_N = 20; % deg
BETA = 18; % deg

l3 = 100; % mm
l6 = 250; % mm

% Forces and moments (N, Nmm)
F_t2 = T / r2;
F_r2 = F_t2 * tand(ALPHA);
T2 = F_t2 * r2;

F_t3 = T / r3;
F_r3 = F_t3 * tand(ALPHA_N) / cosd(BETA);
F_a3 = F_t3 * tand(BETA);
T3 = F_t3 * r3;

F_Ah = (F_r3*(l-l6) + F_a3*r3 - F_r2*(l-l3)) / l;
F_Av = ((F_t2 + G_z2)*(l-l3) + (F_t3 + G_z3)*(l-l6))/l;
F_A = sqrt(F_Ah^2 + F_Av^2);

F_Bh = -F_Ah - F_r2 + F_r3;
F_Bv = F_t2 + G_z2 + F_t3 + G_z3 - F_Av;
F_B = sqrt(F_Bh^2 + F_Bv^2);

% Graphs
x = 0:1:350;
x1 = x(1:l3+1);
x2 = x(l3+1:l6+1);
x3 = x(l6+1:end);

% Horizontal plane

% Qy
Qy1 = @(x) F_Ah + x*0;
Qy2 = @(x) F_Ah + F_r2 + x*0;
Qy3 = @(x) -F_Bh + x*0;
plotGraph({x1, x2, x3}, {Qy1, Qy2, Qy3}, 'Qy', 'r', l, l3, l6);

% Mz
Mz1 = @(x) -F_Ah * x;
Mz2 = @(x) -(F_Ah + F_r2) * (x-l3) + Mz1(l3);
Mz3 = @(x) F_Bh * (x-l6) + F_a3*r3 + Mz2(l6);
plotGraph({x1, x2, x3}, {Mz1, Mz2, Mz3}, 'Mz', 'g', l, l3, l6);

% Tx
Tx1 = @(x) x*0;
Tx2 = @(x) -T2 + x*0;
Tx3 = @(x) x*0;
plotGraph({x1, x2, x3}, {Tx1, Tx2, Tx3}, 'Tx', 'b', l, l3, l6);


% Vertical plane

% Qz
Qz1 = @(x) F_Av + x*0;
Qz2 = @(x) F_Av - F_t2 - G_z2 + x*0;
Qz3 = @(x) -F_Bv + x*0;
plotGraph({x1, x2, x3}, {Qz1, Qz2, Qz3}, 'Qz', 'c', l, l3, l6);

% My
My1 = @(x) F_Av * x;
My2 = @(x) (F_Av - F_t2 - G_z2) * (x - l3) + My1(l3);
My3 = @(x) -F_Bv * (x-l6) + My2(l6);
plotGraph({x1, x2, x3}, {My1, My2, My3}, 'My', 'm', l, l3, l6);


