%% Variable Capacitance

w = pi;

for index = 1:100

    Cap = 0.05*randn + 0.25;
    
    cPerturb(index) = Cap;
    
    C(1,:) = [0 Cap -Cap 0 0 0 0 0];
    C(2,:) = [0 -Cap Cap 0 0 0 0 0];
    C(8,:) = [0 0 0 -L_induct 0 0 0 0];
    
    X = (G + 1i*w*C)\F';
    
    gain(index) = 20*log10(abs(X(8))/V_in);
    
    index = index + 1;

end

figure(5)
hist(cPerturb)
title('Distribution of Capacitance')
xlabel('C')

figure(6)
hist(gain)
title('Spread of Gain')
xlabel('20log_{10}|V_o/V_{in}|')