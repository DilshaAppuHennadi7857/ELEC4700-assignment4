%% DC Sweep

index = 1;

for V_in = -10:10
    
    F = [0 0 0 0 0 V_in 0 0];
    X = G\F';
    
    Vo(index) = X(8);
    V3(index) = X(5)*R3; %X(5) contains I3
    
    index = index + 1;
    
end

V_in = linspace(-10,10,21);

figure(2)
plot(V_in,Vo)
hold on
plot(V_in,V3)
xlabel('Input voltage (V)')
ylabel('Voltage at nodes (V)')
title('Node voltages through DC Sweep')
legend('Vo','V3')
hold off