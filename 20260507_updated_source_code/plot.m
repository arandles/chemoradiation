cd ../sub
data = readmatrix('avg_fractional_change.csv');
T  = data(:,1);   % timesteps
fracChange = data(:,2);   % cell count
dt_per_day = 24 * 60 * 2;   % 2880 timesteps/day
days = T / dt_per_day;
figure;
plot(days, fracChange, 'LineWidth', 2);
xlabel('Days since therapy start');
ylabel('Fractional volume (relative to therapy start)');
title('Tumor Volume Change (Cell Count Proxy)');
grid on;
xlim([0 max(days)]);
ylim([-1 25]);

cd ../opt
data = readmatrix('avg_fractional_change.csv');
T  = data(:,1);   % timesteps
fracChange = data(:,2);   % cell count
dt_per_day = 24 * 60 * 2;   % 2880 timesteps/day
days = T / dt_per_day;
hold all;
plot(days, fracChange, 'r');
ylim([-1 5]);

legend('Suboptimal','Optimal');
