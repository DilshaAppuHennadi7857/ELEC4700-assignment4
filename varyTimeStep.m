% A Gaussian pulse excitation
Cols = hsv(2);

mu = dt*60;
sigma = dt*30;
variance = sigma^2;

prevX = zeros(8,1);

V = [];
index = 1;

f0 = 1/0.03;

for t = 0:numSteps
    if (t < 60)
        V_in = 0;
    else
        V_in = exp(-(t*dt - mu)^2/(2*variance));
    end
%     V_in = sin(2*pi*(1/0.03)*t*dt);

    I_n = 0.001*randn;
    
    F = [0 0 -I_n 0 0 V_in 0 0].';
    
    X = (C/dt + G)\(F + (C*prevX/dt));
    
    prevX = X;
    
    V(1,index) = V_in;
    V(2,index) = X(8);
    V(3,index) = t*dt;
    index = index + 1;
    
end

figure(figNum)
plot(V(3,:),V(1,:),'color',Cols(1,:))
hold on
plot(V(3,:),V(2,:),'color',Cols(2,:))
xlabel('Time')
ylabel('Voltage')
title('Input and Output Voltages with dt = ', dt)
legend('Input','Output')