z_thres = 0.01; 

c_o2 = zeros(size(sol_o2, 1), 1); 
c_co2 = zeros(size(sol_co2, 1), 1);

r_o2 = zeros(size(sol_o2, 1), 1); 
r_co2 = zeros(size(sol_co2, 1), 1); 

for i = 1:length(sol_o2)
    r_o2(i)  =  sqrt(sol_o2(i, 2)^2  + sol_o2(i, 3)^2  ); 
    r_co2(i) =  sqrt(sol_co2(i, 2)^2 + sol_co2(i, 3)^2 ); 
    
   
    c_o2(i) = sol_o2(i, 4); 
    c_co2(i) = sol_co2(i, 4); 
end


t = linspace(0, max(r_o2), 200); 

%ft_o2 = (t.^2-radius^2)+8.5;  % example 1
ft_co2 = -(t.^2-radius^2) + 0.0163;  % example 1 and example 3

ft_o2 = t.^4 - radius^4 + 8.5; % example 3

figure()
subplot(1, 2, 1)
scatter(r_o2, c_o2)
hold on
plot(t, ft_o2)
title('oxygen')
xlabel('radius')
ylabel('concentration')
subplot(1, 2, 2)
scatter(r_co2, c_co2)
hold on
plot(t, ft_co2)
title('carbon dioxide')
xlabel('radius')
ylabel('concentration')



