%% AC Sweep

% Let V_in = 1V
V_in = 1;

Vo = [];
V3 = [];

index = 1;

% frequency sweep
for f = 0:100
    
    w = 2*pi*f;
    X = (G + 1i*w*C)\F';
    
    Vo(index) = real(X(8));
    V3(index) = real(X(5))*R3; %X(5) contains I3
    
    index = index + 1;
    
end

f = linspace(0,100,101);

figure(3)
plot(f,Vo)
hold on
plot(f,V3)
xlabel('Frequency (Hz)')
ylabel('Voltage(V)')
title('Voltage at Nodes through Frequency Sweep')
legend('Vo','V3')
hold off

figure(4)
plot(f,20*log10(Vo/V_in))
xlabel('Frequency (Hz)')
ylabel('Gain (dB)')
title('Frequency Response')