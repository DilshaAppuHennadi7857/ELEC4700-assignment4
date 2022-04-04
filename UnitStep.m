% A unit step that steps up from 0 to 1 at 0.03 seconds
Cols = hsv(2);

prevX = zeros(8,1);

V = [];
index = 1;

for t = 0:numSteps
    
    % \[V(t) = \frac{\left[F + \frac{C}{\Delta t} V(t - \Delta t) \right]}{(\frac{C}{\Delta t} + G)}\]
    % where V_in = 0 when t < 0.03s, then jumps to V_in = 1 at 0.03s
    if (t*dt < 0.03)
        V_in = 0;
    else
        V_in = 1;
    end
    
    F = [0 0 0 0 0 V_in 0 0].';
    
    X = (C/dt + G)\(F + (C*prevX/dt));
    
    prevX = X;
    
    V(1,index) = V_in;
    V(2,index) = X(8);
    index = index + 1;
    
    figure(7)
    plot(t*dt,V_in,'*','color',Cols(1,:))
    hold on
    plot(t*dt,X(8),'O','color',Cols(2,:))
    xlabel('Time (s)')
    ylabel('Voltage (V)')
    title('Input and Output Voltages')
    legend('Input','Output')
    
    pause(0.01)
    
end

hold off

fs = 1/dt;
n = 2^nextpow2(numSteps);
f = (-n/2:n/2-1)*(fs/n);

VinX = fft(V(1,:),n);
VinY = fftshift(VinX);
VinPow = abs(VinY).^2/n;

VoX = fft(V(2,:),n);
VoY = fftshift(VoX);
VoPow = abs(VoY).^2/n;

% fft()

figure(8)
subplot(2,1,1); plot(f,VinX);
title('Single-Sided Amplitude Spectrum of Vin(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

subplot(2,1,2); plot(f,VoX);
title('Single-Sided Amplitude Spectrum of VO(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% fftshift()

figure(9)
subplot(2,1,1); plot(f,VinPow);
title('Amplitude Spectrum of Vin(t)')
xlabel('f (Hz)')
ylabel('|Pin(f)|')

subplot(2,1,2); plot(f,VoPow);
title('Amplitude Spectrum of VO(t)')
xlabel('f (Hz)')
ylabel('|Po(f)|')