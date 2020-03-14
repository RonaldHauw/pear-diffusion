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


figure()
subplot(1, 2, 1)
scatter(r_o2, c_o2)
title('oxygen')
xlabel('radius')
ylabel('concentration')
subplot(1, 2, 2)
scatter(r_co2, c_co2)
title('carbon dioxide')
xlabel('radius')
ylabel('concentration')



