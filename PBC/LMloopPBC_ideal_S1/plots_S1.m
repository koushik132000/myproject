T = readtable("mderend.csv");
T_array = table2array(T); % removes the header 

X = (1:1:41)'; 
Y  = T_array(1:41,1); % extracts the column
err  = T_array(1:41,5)';

Y2  = T_array(1:41,2);
err2 = T_array(1:41,6);
figure('Position',[100 100 900 400]);% make figure wider\

% Error bars only
e1 = errorbar(X, Y, err, 'LineStyle', 'none');
e1.Color = [1 0 0];
e1.LineWidth = 1;

hold on;

% Bold line
p1 = plot(X, Y, 'Color', [1 0 0], 'LineWidth', 2);

% Axes limits
ylim([0 20])
xlim([0 45])

% Axis labels
xlabel("MC runs");
ylabel("Rrms");

%Custom x-axis ticks at 0,2,4,...,40
xticks(0:2:40);
yticks(0:2:20);

% Darken the plot box
ax = gca;                % get current axes
ax.LineWidth = 1.5;      % make box border thicker/darker
ax.XColor = 'k';         % set x-axis color (black)
ax.YColor = 'k';         % set y-axis color (black)
ax.Position = [0.1 0.15 0.85 0.75];   % stretch axes


figure('Position',[100 100 900 400]);% make figure wider\
Y  = T_array(1:41,2); % extracts the column
err  = T_array(1:41,6)';
% Error bars only
e1 = errorbar(X, Y, err, 'LineStyle', 'none');
e1.Color = [1 0 0];
e1.LineWidth = 1;

hold on;

% Bold line
p2 = plot(X, Y, 'Color', [1 0 0], 'LineWidth', 2);

% Axes limits
ylim([0 20])
xlim([0 45])

% Axis labels
xlabel("MC runs");
ylabel("Rrms");

%Custom x-axis ticks at 0,2,4,...,40
xticks(0:2:40);
yticks(0:2:20);

% Darken the plot box
ax = gca;                % get current axes
ax.LineWidth = 1.5;      % make box border thicker/darker
ax.XColor = 'k';         % set x-axis color (black)
ax.YColor = 'k';         % set y-axis color (black)
ax.Position = [0.1 0.15 0.85 0.75];   % stretch axes
