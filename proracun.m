classdef Proracun
    properties
        % Sample values
        % T = 610000; % Nmm
        % J2 = 0.0500; % kgm^2 = Nms^2
        % J3 = 0.8750; % kgm^2 = Nms^2
        % l = 360; % mm
        % G_z2 = 240; % N
        % G_z3 = 110; % N
        % b2 = 120; % mm
        % b3 = 120; % mm
        % r2 = 180; % mm
        % r3 = 63.1; % mm

        % sigmaF_max = 50; % N/mm^2
        % sigma_fDN = 240; % N/mm^2
        % tau_fDI = 190; % N/mm^2
        % alpha0;

        % n = 5.33; % rps
        % L_h = 8000; % hrs
        % ALPHA = 20; % deg
        % ALPHA_N = 20; % deg
        % BETA = 18; % deg

        % l3 = 108; % mm
        % l6 = 252; % mm
        
        % Custom values
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

        sigmaF_max = 50; % N/mm^2
        sigma_fDN = 240; % N/mm^2
        tau_fDI = 190; % N/mm^2
        alpha0;

        n = 5.33; % rps
        L_h = 8000; % hrs
        ALPHA = 20; % deg
        ALPHA_N = 20; % deg
        BETA = 18; % deg

        l3 = 100; % mm
        l6 = 250; % mm

        resolution = 0.5;

        F_t2; F_r2; T2;

        F_t3; F_r3; F_a3; T3;

        F_Ah; F_Av; F_A;

        F_Bh; F_Bv; F_B;

        x; x1; x2; x3;

        Qy1; Qy2; Qy3;
        Mz1; Mz2; Mz3;
        Tx1; Tx2; Tx3;
        Qz1; Qz2; Qz3;
        My1; My2; My3;
    end

    methods
        function obj = Proracun()
            obj.F_t2 = obj.T / obj.r2;
            obj.F_r2 = obj.F_t2 * tand(obj.ALPHA);
            obj.T2 = obj.F_t2 * obj.r2;

            obj.F_t3 = obj.T / obj.r3;
            obj.F_r3 = obj.F_t3 * tand(obj.ALPHA_N) / cosd(obj.BETA);
            obj.F_a3 = obj.F_t3 * tand(obj.BETA);
            obj.T3 = obj.F_t3 * obj.r3;

            obj.F_Ah = (obj.F_r3*(obj.l-obj.l6) + obj.F_a3*obj.r3 - obj.F_r2*(obj.l-obj.l3)) / obj.l;
            obj.F_Av = ((obj.F_t2 + obj.G_z2)*(obj.l-obj.l3) + (obj.F_t3 + obj.G_z3)*(obj.l-obj.l6))/obj.l;
            obj.F_A = sqrt(obj.F_Ah^2 + obj.F_Av^2);

            obj.F_Bh = -obj.F_Ah - obj.F_r2 + obj.F_r3;
            obj.F_Bv = obj.F_t2 + obj.G_z2 + obj.F_t3 + obj.G_z3 - obj.F_Av;
            obj.F_B = sqrt(obj.F_Bh^2 + obj.F_Bv^2);
            
            obj.x = 0:obj.resolution:obj.l;
            l3_index = find(obj.x == obj.l3);
            l6_index = find(obj.x == obj.l6);
            obj.x1 = obj.x(1:l3_index+1);
            obj.x2 = obj.x(l3_index+1:l6_index+1);
            obj.x3 = obj.x(l6_index+1:end);

            obj.Qy1 = @(x) obj.F_Ah + x*0;
            obj.Qy2 = @(x) obj.F_Ah + obj.F_r2 + x*0;
            obj.Qy3 = @(x) -obj.F_Bh + x*0;

            obj.Mz1 = @(x) -obj.F_Ah * x;
            obj.Mz2 = @(x) -(obj.F_Ah + obj.F_r2) * (x-obj.l3) + obj.Mz1(obj.l3);
            obj.Mz3 = @(x) obj.F_Bh * (x-obj.l6) + obj.F_a3*obj.r3 + obj.Mz2(obj.l6);
            
            obj.Tx1 = @(x) x*0;
            obj.Tx2 = @(x) -obj.T2 + x*0;
            obj.Tx3 = @(x) x*0;

            obj.Qz1 = @(x) obj.F_Av + x*0;
            obj.Qz2 = @(x) obj.F_Av - obj.F_t2 - obj.G_z2 + x*0;
            obj.Qz3 = @(x) -obj.F_Bv + x*0;

            obj.My1 = @(x) obj.F_Av * x;
            obj.My2 = @(x) (obj.F_Av - obj.F_t2 - obj.G_z2) * (x-obj.l3) + obj.My1(obj.l3);
            obj.My3 = @(x) -obj.F_Bv * (x-obj.l6) + obj.My2(obj.l6);

            obj.alpha0 = obj.sigma_fDN / (sqrt(3) * obj.tau_fDI);
        end

        function plotIdealAxle(obj)
            figure('Name', 'Idealni oblik vratila')
            grid;
            hold on;

            xline(obj.l3, '--');
            xline(obj.l6, '--');
            xline(obj.l, '--');

            plot(obj.x, obj.getDiameter(obj.x));
            xlabel('Udaljenost (mm)');
            ylabel('Promjer vratila (mm)');

            hold off
        end

        function d = getDiameter(obj, x)
            k = (10 / obj.sigmaF_max)^(1/3);
            d = zeros(size(x));

            for i=1:length(x)
                if x(i) <= obj.l3
                    M = sqrt(obj.Mz1(x(i))^2 + obj.My1(x(i))^2);
                elseif x(i) <= obj.l6
                    Mf = sqrt(obj.Mz2(x(i))^2 + obj.My2(x(i))^2);
                    M = sqrt(Mf^2 + 0.75 * (obj.alpha0 * obj.Tx2(x(i)))^2);
                else
                    M = sqrt(obj.Mz3(x(i))^2 + obj.My3(x(i))^2);
                end

                d(i) = k * M^(1/3);
            end
        end

        function plotForcesMoments(obj)
            figure('Name', 'Sile i momenti')

            graphs = {
                {obj.My1, obj.My2, obj.My3}, 'm', 'My', 'Moment (Nmm)'
                {obj.Mz1, obj.Mz2, obj.Mz3}, 'g', 'Mz', 'Moment (Nmm)';
                {obj.Tx1, obj.Tx2, obj.Tx3}, 'b', 'Tx', 'Moment (Nmm)';
                {obj.Qy1, obj.Qy2, obj.Qy3}, 'r', 'Qy', 'Sila (N)';
                {obj.Qz1, obj.Qz2, obj.Qz3}, 'c', 'Qz', 'Sila (N)';
            };
            for i=1:size(graphs, 1)
                subplot(2, 3, i)
                obj.plotGraph(graphs{i, 1}, graphs{i, 2});
                title(graphs{i, 3});
                xlabel('Udaljenost (mm)');
                ylabel(graphs{i, 4});
                grid;
            end
        end

        function plotGraph(obj, y, color)
            hold on;
        
            xline(obj.l3, '--');
            xline(obj.l6, '--');
            xline(obj.l, '--');

            xData = {obj.x1, obj.x2, obj.x3};

            for i=1:length(xData)
                curX = xData{i};
                curY = y{i};
                plot(curX, curY(curX), color, 'LineWidth', 2)
                fill([curX, fliplr(curX)], [curY(curX), zeros(1, length(curY(curX)))], color, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
            end
        
            hold off;
        end
    end
end