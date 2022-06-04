close all
clear, clc
% alpha = 10; T = 0.01;
% alpha = 5; T = 0.01;
% alpha = 3; T = 0.001;
% alpha = 20; T = 0.001;

alpha = 0.1; T = 2;
alpha = 0.1; T = 10;

%reachsetdyn(alpha, 1, 2, 20, 'C:\Users\User\Desktop\edu\edu_3\оу\task2\mymovie.avi');
%alpha = 1; T = 1;
%alpha = 0.1; T = 1;



%return


[bord, traj] = reachset (alpha, T, 20);

cur_figure = figure(1);
cur_figure.Position = [50 50 900 900];
ax = axes(cur_figure, 'NextPlot', 'add', 'XGrid', 'on', 'YGrid', 'on');

for i = 1:size(traj, 1)
    y1 = reshape(traj(i, 1, :), size(traj, 3), 1);
    y2 = reshape(traj(i, 2, :), size(traj, 3), 1);
    ind1 = find(y1 == 0);
    ind2 = find(y2 == 0);
    ind = find(ind1 == ind2);
    %plot(ax, traj(i, 1, :), traj(i, 2, :));
    col = [i/size(traj, 1), i/size(traj, 1), (size(traj, 1) - i)/size(traj, 1) ];
    plot(ax, y1(1: ind1(ind(2))-1), y2(1: ind1(ind(2))-1), 'Color', col)
end
return

if T>=1
    my_i = find (bord(:,1) > 0)
    %plot(ax, bord(my_i, 1), bord(my_i, 2), 'g', 'Linewidth', 10);
    bord(my_i,:) = bord(my_i - 1,:);
end 
plot(ax, bord(:, 1), bord(:, 2), 'g', 'Linewidth', 2);
%fill(bord(:, 1), bord(:, 2), 'g')
bord
i
p = fill(bord(:, 1), bord(:, 2), [0.7 1 0.8]);
p.LineWidth = 1.5;
p.EdgeColor = [0.2 0.8 0.55];
p.FaceAlpha = 0.5;
ylabel('x_1(t) = x''')
xlabel('x_2(t) = x')
str = strcat('Множество достижимости при \alpha = ', num2str(alpha));
str = strcat(str, ', T = ');
str = strcat(str, num2str(T));
title(str)

return






% ---------------------------- ФУНКЦИИ ----------------------------

function [border, trajectory] = reachset(alpha, T, N)
eps = 0.001;
t0 = 0;
%N = 10; % перебор переключений
trajectory = zeros(N*N, 2, 200);
cur_traj = [1, 1];

myconv = zeros(61, 2);
iconv = 1;
global cur_u
%cur_figure = figure(1);
%cur_figure.Position = [50 50 800 500];
%ax = axes(cur_figure, 'NextPlot', 'add', 'XGrid', 'on', 'YGrid', 'on');
thresh = 0;      % Threshold value
tol = 1e-2;      % Tolerance on threshold
opts = odeset('Events',@(t,y)event_switch(t,y,thresh,tol)); % Create events function

for twocycles = 1:2
    cur_u = alpha;
    [t, y, te, yend, ie] = ode45(@odeS,[t0, T],[0,0], opts);
    if size(te, 1) >1
        te = te(1);
        yend = yend(1, :);
    end
    
    if (size(te, 1) == 0)
        %plot(ax, y(:, 1), y(:, 2),  'Color', [0; 0; 1]);
        %trajectory(cur_traj(1), :, cur_traj(2):cur_traj(2)+size(y, 1)) = y;
        trajectory(cur_traj(1), 1, cur_traj(2):cur_traj(2)+size(y, 1)-1) = y(:, 1);
        trajectory(cur_traj(1), 2, cur_traj(2):cur_traj(2)+size(y, 1)-1) = y(:, 2);
        %if cur_traj(2) < 200
         %   trajectory(cur_traj(1), 2, cur_traj(2):end) = ones(1, 1, 200 - cur_traj(2)+1).*trajectory(cur_traj(1), 2, cur_traj(2):end);
          %  trajectory(cur_traj(1), 1, cur_traj(2):end) = ones(1, 1, 200 - cur_traj(2)+1).*trajectory(cur_traj(1), 1, cur_traj(2):end);   
        %end
        cur_traj(2)=1;
        cur_traj(1) = cur_traj(1)+1;
        
        myconv(iconv, :) = y(end, :);
        iconv = iconv + 1;
    else
        % Перебор по времени преключения
        for tswitch = t0+te/100:te/N:T
            iend = find(t == tswitch);
            [t, y] = ode45(@odeS,[t0, tswitch],[0,0]);
            col_traj = [tswitch./te; 0; (te-tswitch)/te];
            
            %plot(ax, y(:, 1), y(:, 2), '-',  'Color', col_traj);
            %cur_traj(2)
            %cur_traj(2)+size(y, 1)
            %size(y)
            %size(trajectory(cur_traj(1), :, cur_traj(2):cur_traj(2)+size(y, 1)-1))
            trajectory(cur_traj(1), 1, cur_traj(2):cur_traj(2)+size(y, 1)-1) = y(:, 1);
            trajectory(cur_traj(1), 2, cur_traj(2):cur_traj(2)+size(y, 1)-1) = y(:, 2);
            cur_traj(2) = cur_traj(2)+size(y, 1);
            
            tend1 = tswitch;
            while tend1 < T
                [tnew, psinew, tend, psiend, iend] = ode45(@odefun_psi,[tend1, T],[sign(cur_u),0], opts);
                if (size(tend,1) == 0)
                    tend = T;
                    iend = size(tnew, 1);
                end
                if (size(tend,1) > 1)
                    tend = tend(1);
                    iend = iend(1);
                    yend = yend(1, :);
                end
                
                cur_u = cur_u*(-1);
                [tnew, y] = ode45(@odeS,[tend1, tend],y(end, :), opts);
                %plot(ax, y(:, 1), y(:, 2), 'Color', col_traj);
                trajectory(cur_traj(1), 1, cur_traj(2):cur_traj(2)+size(y, 1)-1) = y(:, 1);
                trajectory(cur_traj(1), 2, cur_traj(2):cur_traj(2)+size(y, 1)-1) = y(:, 2);
                cur_traj(2) = cur_traj(2)+size(y, 1);
            
                tend1 = tend;
            end
            
            %if cur_traj(2) < 200
             %   trajectory(cur_traj(1), 2, cur_traj(2):end) = ones(1, 1, 200 - cur_traj(2)+1).*trajectory(cur_traj(1), 2, cur_traj(2):end);
              %  trajectory(cur_traj(1), 1, cur_traj(2):end) = ones(1, 1, 200 - cur_traj(2)+1).*trajectory(cur_traj(1), 1, cur_traj(2):end);   
            %end
            cur_traj(1) = cur_traj(1) + 1;
            cur_traj(2) = 1;
            
            myconv(iconv, :) = y(end, :);
            iconv = iconv + 1;
        end
    end
    alpha = -alpha;
    iconv = 31;
end

%myconv
%find(myconv ~= 0, 1)
%if size(find(myconv ~= 0, 10, ) < 2)
 %   find(myconv ~= 0, 1)
  %  border = myconv(find(myconv ~= 0, 1), :);
%else 

    [k, ~] = convhull(myconv);
    border = myconv(k, 1:2);
%end
%p = fill(myconv(k, 1), myconv(k, 2), [0.7 1 0.8]);
%p.LineWidth = 1.5;
%p.EdgeColor = [0.2 0.8 0.55];
%p.FaceAlpha = 0.1;
%ylabel('x_1(t) = x''')
%xlabel('x_2(t) = x')
%str = strcat('Множество достижимости при \alpha = ', num2str(alpha));
%str = strcat(str, ', T = ');
%str = strcat(str, num2str(T));
%title(str)

end


function [value,isterminal,direction] = reachsetdyn(alpha,t1, t2, N, filename)
v = VideoWriter(filename);
open(v);
v1 = matfile('C:\Users\User\Desktop\edu\edu_3\прак\my_animation','Writable',true);
cur_figure = figure(1);
cur_figure.Position = [50 50 900 900];
ax = axes(cur_figure, 'NextPlot', 'add', 'XGrid', 'on', 'YGrid', 'on');


for T = t1:(t2-t1)/N:t2
clf
cur_figure = figure(1);
cur_figure.Position = [50 50 900 900];
ax = axes(cur_figure, 'NextPlot', 'add', 'XGrid', 'on', 'YGrid', 'on');
[bord, traj] = reachset (alpha, T, 20);


for i = 1:size(traj, 1)
    y1 = reshape(traj(i, 1, :), size(traj, 3), 1);
    y2 = reshape(traj(i, 2, :), size(traj, 3), 1);
    ind1 = find(y1 == 0);
    ind2 = find(y2 == 0);
    ind = find(ind1 == ind2);
    %plot(ax, traj(i, 1, :), traj(i, 2, :));
    col = [i/size(traj, 1), i/size(traj, 1), (size(traj, 1) - i)/size(traj, 1) ];
    plot(ax, y1(1: ind1(ind(2))-1), y2(1: ind1(ind(2))-1), 'Color', col)
end
plot(ax, bord(:, 1), bord(:, 2), 'g', 'Linewidth', 2);

y1lim = max(abs(reshape(traj(1, 1, :), size(traj, 3), 1)));
y2lim = max(abs(reshape(traj(1, 1, :), size(traj, 3), 1)));

xlim([-y1lim y1lim].*2)
ylim([-y2lim y2lim].*2)

%fill(bord(:, 1), bord(:, 2), 'g')
p = fill(bord(:, 1), bord(:, 2), [0.7 1 0.8]);
p.LineWidth = 1.5;
p.EdgeColor = [0.2 0.8 0.55];
p.FaceAlpha = 0.5;
ylabel('x_1(t) = x''')
xlabel('x_2(t) = x')
str = strcat('Множество достижимости при \alpha = ', num2str(alpha));
str = strcat(str, ', T = ');
str = strcat(str, num2str(T));
title(str)
pause(0.02)

tmpmov = getframe();
size(tmpmov.cdata)
tmpmov.cdata = tmpmov.cdata (1:734, 1:690, :)
writeVideo(v,tmpmov);
end

close(v);
end




% psi_2 = 0
% x_2 = 0
function [value,isterminal,direction] = event_switch(~,y,thresh,tol)
value = [y(2); y(2)]-thresh-tol;
isterminal = [0; 0]; % Change termination condition
direction = [1; 1];
end

function dydt = odeS1(~,y)
global cur_u;
dydt = zeros(2,1);
dydt(1) = y(2);
dydt(2) = cur_u - y(2) - y(1).*(2 - sin(y(1).^2)) + 2.*y(1).^2.*cos(y(1));
end


function dydt = odefun_psi1(~,y)
dydt = zeros(2,1);
dydt(1) = 2.*y(2) + sin(y(1).^2) + 2.*y(1).^2.*cos(y(1).^2) + 2.*y(1).^2.*sin(y(1)) - 4.*y(1).*cos(y(1));
dydt(2) = y(2) - y(1);
end


function dydt = odeS(~,x)
global cur_u;
dydt = [x(2); cur_u - x(1) + sin(x(2)) - x(1).*sin(x(1).^3)];
end


function dydt = odefun_psi(~,x)
global cur_u;
dydt = [x(2); cur_u - x(1) + sin(x(2)) - x(1).*sin(x(1).^3)];
end




