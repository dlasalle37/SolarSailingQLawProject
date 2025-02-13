clear; close all

withJ2 = readmatrix("data\plotting_data\CASE_E_WITH_J2.txt");
noJ2 = readmatrix("data\plotting_data\CASE_E_NO_J2.txt");
withJ2_oe = readmatrix("data\plotting_data\CASE_E_WITH_J2_oe.txt");
noJ2_oe = readmatrix("data\plotting_data\CASE_E_NO_J2_oe.txt");
t = 0:60:(length(noJ2)-1)*60;
N = length(withJ2);
diff = withJ2(:, 1:3) - noJ2(1:N,1:3);

oet = [26700.0, 0.75, 0.01*pi/180, 0.0, 90];
Woe = [1, 1, 1, 0, 0];

%normalize each oe with oet
for row = 1:length(noJ2_oe)
   noJ2_oe(row, 1:5) = noJ2_oe(row, 1:5) ./ oet;
end
for row = 1:length(withJ2_oe)
   withJ2_oe(row, 1:5) = withJ2_oe(row, 1:5) ./ oet;
end

figure
subplot(2,1,1)
plot(t, noJ2_oe(:,1:2) ,'linewidth', 1.5)
legend("a", "e", "i", "\omega", "location", "eastoutside")
title("Normalized Targeted elements without Perturbations")
ylabel("oe/oe_t")
xlabel("Time(s)")
xlim([t(1), t(end)]);

subplot(2,1,2)
plot(t(1:N), withJ2_oe(:,1:2), 'linewidth', 1.5)
legend("a", "e", "i", "\omega", "location", "eastoutside")
title("Normalized Targeted elements with Perturbations")
ylabel("oe/oe_t")
xlabel("Time(s)")
xlim = [t(1), t(end)];
