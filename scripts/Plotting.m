clear; close all

ecl = readmatrix("data\eclipsed.txt"); % coasting arcs
cart = readmatrix("data\cart.txt"); % non-coasting arcs


figure
hold on
plot3(cart(:,1), cart(:,2), cart(:,3))
plot3(ecl(:,1), ecl(:,2), ecl(:,3), '.r', 'markersize', 3)
axis equal
