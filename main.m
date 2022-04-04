%% ELEC 4700 Assignment 4: Circuit Modeling
%
% Dilsha Appu-Hennadi, 101107857
% Apr. 3, 2022
% Latest rev. Apr. 3, 2022

clear all
clearvars
clearvars -GLOBAL
close all
format shorte

global regL regW  nx ny L W

regL = 200e-9; % length of area
regW = 100e-9; % width of area

nx = 50; % number of space in x direction
ny = 50; % number of spaces in y direction

L = linspace(0, regL, nx); % nx spaces from 0 to region length
W = linspace(0, regW, ny); % ny spaces from 0 to region width

q = 1.60217662e-19; %C % elementary charge
m_e = 9.10938356e-31; %kg % rest mass of an electron
m_n = 0.26*m_e; %kg % mass of an electron
kb = 1.38064852e-23; %m^2 kg s^-2 K^-1 % Boltzmann's Constant

T = 300; %K % temperature

numElec = 800;%50;
numDispElec = 5;
t_mn = 0.2e-12; %s % mean time between collisions

numTimeStep = 800; %50;
dt = 1e-16;

%% Part 1
%
% The circuit shown in Figure 1 is the same circuit simulated in PA 7.
% However, this time there is no value for R3 given. Instead, a value for
% R3 will be determined through taking a voltage sweep of the bottle-neck
% simulation from Assignment 3, plotting the current-voltage
% characteristics, and taking the slope of a linear fit.
%

index = 1;

for V_0 = 0.1:0.1:10 % run voltage sweep
    
    icCharacteristics
    Voltage(index) = V_0;
    deviceCurrent(index) = Current;
    
    index = index + 1;
    
end

% conduct linear fit
p = polyfit(Voltage, deviceCurrent,1);
v = linspace(0.1,10,100);
c = polyval(p,v);

figure(1)
plot(Voltage, deviceCurrent, '*')
hold on
plot(v,c)
hold off
title('I-V Characteristics of R3')
xlabel('Voltage (V)')
ylabel('Current (A)')
legend('Data points', 'Linear fit')

%%
%
% The figure above depicts the currents through R3 as voltage is swept from
% 0.1V to 10V using 0.1V increments. By plotting this data, we see that
% there is a clear linear relationship between the voltage and current
% through R3. By fitting the data to a 1st-order polynomial, we obtain the
% coefficients for a linear plot describing this relationship. A value for
% R3 can be obtained by taking the inverse of this slope such that:
%
% \[R3~=~ slope^{-1} ~=~ \frac{voltage}{current} ~\approx~ 16 \Omega\]
%
% Below are the values for all components in the circuit:
%

V_in = 1;

R1 = 1
Cap = 0.25
R2 = 2
L_induct = 0.2
R3 = 1/p(1)
alpha = 100
R4 = 0.1
Ro = 1000

G1 = 1/R1
G2 = 1/R2
G3 = 1/R3
G4 = 1/R4
Go = 1/Ro

%% DC Sweep
%
% To begin analysis of the circuit, Kirchhoff's Current Law may be used to
% determine a set of differential equations that represent the network:
% \[1\:~I_{in} + G1V1 + C dV1/dt - G1V2 - C dV2/dt = 0\]
% \[2\:~-G1V1 - C dV1/dt + (G1+G2)V2 + C dV2/dt + I_L = 0\]
% \[3\:~- I_L + I3 = 0\]
% \[4\:~-I + G4V4 - G4VO = 0\]
% \[5\:~-G4V4 + (G4+GO)VO = 0\]
% \[6\:~V1 = V_{in}\]
% \[7\:~- \alpha I3 + V4 = 0\]
% \[8\:~V2 - L dI_L/dt - R3I3 = 0\]
%
% These equations can be writen in the frequency domain as well:
% \[1\:~I_{in} + G1V1 + Cj \omega V1 - G1V2 - Cj \omega V2 = 0\]
% \[2\:~-G1V1 - Cj \omega V1 + (G1+G2)V2 + Cj \omega V2 + I_L = 0\]
% \[3\:~- I_L + I3 = 0\]
% \[4\:~-I + G4V4 - G4VO = 0\]
% \[5\:~-G4V4 + (G4+GO)VO = 0\]
% \[6\:~V1 = V_{in}\]
% \[7\:~- \alpha I3 + V4 = 0\]
% \[8\:~V2 - Lj \omega I_L - R3I3 = 0\]
%
% Our unknowns are: X = [I_in, V1, V2, I_L, I3, I, V4, Vo];

X = [];

G = zeros(8);
G(1,:) = [1 G1 -G1 0 0 0 0 0];
G(2,:) = [0 -G1 (G1+G2) 1 0 0 0 0];
G(3,:) = [0 0 0 -1 1 0 0 0];
G(4,:) = [0 0 0 0 0 -1 G4 -G4];
G(5,:) = [0 0 0 0 0 0 -G4 (G4+Go)];
G(6,:) = [0 1 0 0 0 0 0 0];
G(7,:) = [0 0 0 0 -alpha 0 1 0];
G(8,:) = [0 0 1 0 -R3 0 0 0];

C = zeros(8);
C(1,:) = [0 Cap -Cap 0 0 0 0 0];
C(2,:) = [0 -Cap Cap 0 0 0 0 0];
C(8,:) = [0 0 0 -L_induct 0 0 0 0];

F = [0 0 0 0 0 V_in 0 0];

%%
%
% Thus our C and G matricies are:
%

C

G

%%
%
% With matrices C G and vector, F, determined, we can run a DC sweep with
% voltage, \(V_{in}\) sweeping from -10V to 10V. This can be done using the
% frequency domain equations with \( \omega = 0\).
%

DCsweep

%%
%
% In Figure 2 above, we can see how the voltages V3 and VO change through
% the voltage sweep. Both voltages vary linearly as the input voltage
% increases. Since VO is dependdent on \(\alpha I3\), it makes sense that
% if V3 is linear than VO must be linear. We also see that the slope of VO
% is steeper than V3, indicating that the signal is being amplified by the
% device.
%

%% AC Sweep
%
% Similarly, we can also conduct a frequency sweep. To do this, we hold
% \(V_{in}\) at 1V then loop through different frequency values. VO and V3
% are plotted with the change in frequency. The gain, in dB, can also be
% plotted.
%

ACsweep

%%
%
% We can also add some random perturbations to the
% capacitor to simulate a variable capacitor. For this, the standard
% deviation is assumed to be 0.05 and the frequency is assumed as \( \omega
% = \pi \).
%

Cperturb

%% Transient Circuit Simulation
%
% Looking back at the circuit depicted in Figure 1 of the assignment
% instructions, given the components within the circuit and their
% configuration, it seems that the circuit being modelled is a filter in
% series with an amplifier. Where R1 and C seem to model a capacitor with a
% lossy dielectric. This capacitor along with R2 and L would create a
% bandpass filter. R3, the current dependent voltage source and R4 would
% model the transistor amplifier and RO is the load.
%
% As mentioned earlier, as this circuit has both a capacitor and inductor,
% we can expect to see a bandpass type behaviour as both the capacitor and
% inductor will have an associateed cut-off frequency.

%%
%
% In the time domain, the circuit can be represented by the equation:
% \[C \frac{dV}{dt} + GV = F\]
% As a first order ODE this equation is a bit difficult to solve as it
% currently is. Instead, we can write it in its finite difference form and
% determine a numerical solution.
%
% \[C \frac{V(t) - V(t - \Delta t)}{\Delta t} + G V(t) = F\]
% \[\frac{C}{\Delta t} V(t) - \frac{C}{\Delta t} V(t - \Delta t) + G V(t) = F\]
% \[(\frac{C}{\Delta t} + G) V(t) = F + \frac{C}{\Delta t} V(t - \Delta t)\]
% \[V(t) = \left[F + \frac{C}{\Delta t} V(t - \Delta t) \right](\frac{C}{\Delta t} + G)^{-1}\]
%

%%
%
% For the purposes of the FDTD, we will use 1000 times steps to simulate 1
% second - i.e. each time step will be a milisecond. Three different inputs
% will be used:
% \being{itemize}
%     \item A unit step that steps up from 0 to 1 at 0.03 seconds
%     \item A sinusoid, \(\sin(2 \pi f t)\), with f = 1/0.03 Hz and various
%     other frequencies
%     \item A gaussian pulse with a magnitude of 1 and std. dev. of 0.03
%     seconds and delay of 0.06 seconds.
% \end{itemize}
% It should also be noted, that as previous values are being used in the
% calculation, all initial values will be assumed to be 0.
%

numSteps = 1000;
dt = 1e-3; %50e-3;

UnitStep

%%
%
% Figure 7, above, depicts the input and output signals of the circuit over
% time. The input signal, which is a unit step at 0.06s, starts at 0V and
% at t = 0.06s, the signal jumps to 1V. This causes a transient response in
% the graph as the output tries to catch up to the input resulting in an
% overshoot which decreases and settles off at a final value after some
% time. As we saw earlier in the DC sweep graphs the output signal has a
% much larger amplitude than the input signal as the circuit amplifies the
% input.
%
% Figure 8 and 9 show the single-sided amplitude spectrum and the two-sided
% spectrum respectively. As the unit step is mainly a jump at some time t,
% we only see a single delta function within the amplitude spectrum graphs
% for both signals.
%

Sinusoid

%%
%
% Figures 10 shows the input and output wave forms of the circuit when a
% sinusoidal excitations is used. This input sinusiod has an amplitude of 1
% and a frequency of about 33.33Hz. The output response from the circuit
% overshoots at the beginning as the circuit rushes to catch up to the
% input but quickly levels out at a steady state. Once again, we see the
% amplitude of the output is larger than the input, indicating an
% amplification of the input signal.
%
% Figure 11 is the single-sided amplitude spectrum, while Figure 12 is the
% two-siedd spectrum. While it's a bit unclear in Figure 11, Figure 12
% indicates the frequency of both the input and output as about 33Hz, which
% is to be expected.
%

SinusoidV2

%%
%
% Figures 13 to 15 depict a similar situation to the previous excitation.
% This time, however, the sinusoid has a frequency of 10Hz. The lower
% frequency allows us to see the signals more clearly and we also see that
% there is a delay between the input and the output as the peaks/troughs of
% the two singnals do not occur at the same time.
%
% Additionally, looking at the results of the fft in Figures 14 and 15, the
% frequency of the input and output is found to be about 10Hz. This is what
% we expect it to be.
%

GaussianPulse

%%
%
% Figure 16 depicts the input Gaussian excitation and the corresponding
% output from the circuit. As seen earlier, there is a delay in the output.
% Then, after the Gaussian pulse ends, the output overshoots and dips below
% 0V before settle out to its steady state of 0V.
%
% From the FFT of the signals, the frequency seems to be about 2Hz
%

%%
%
% So far the FFTs of the signals have been reasonably accurate to the
% fequency expected. However, as the step time, dt, is increased, we see
% the FFTs become less accurate.
%
% We can think of the step time as the period of the sampling frequency -
% i.e.:
% \[f_s = \frac{1}{dt}\]
% So we see, as the step time increases the sampling frequency decreases.
% According to the Nyquist-Shannon Sampling Theorem, in order to recreate a
% signal from a samples of the original signal, the sampling frequency must
% be more than double the frequency of the original signal:
% \[f < 2 f_s\]
%

%% Part 2
%
% For this part, the previous circuit is now slightly modified. It now
% includes a current source, \( I_n \) to model thermal noise and a
% capacitor, \( C_n \), to bandwidth limit the noise. With these additions,
% we must make some slight modifications to our set of equations:
% \[1\:~I_{in} + G1V1 + C dV1/dt - G1V2 - C dV2/dt = 0\]
% \[2\:~-G1V1 - C dV1/dt + (G1+G2)V2 + C dV2/dt + I_L = 0\]
% \[3\:~- I_L + G3V3 + I_n + Cn dV3/dt = 0\]
% \[4\:~-I + G4V4 - G4VO = 0\]
% \[5\:~-G4V4 + (G4+GO)VO = 0\]
% \[6\:~V1 = V_{in}\]
% \[7\:~- \alpha G3V3 + V4 = 0\]
% \[8\:~V2 - L dI_L/dt - V3 = 0\]
%
% The frequency domain equations can then be written as:
% \[1\:~I_{in} + G1V1 + Cj \omega V1 - G1V2 - Cj \omega V2 = 0\]
% \[2\:~-G1V1 - Cj \omega V1 + (G1+G2)V2 + Cj \omega V2 + I_L = 0\]
% \[3\:~- I_L + G3V3 + I_n + Cn j \omega V3 = 0\]
% \[4\:~-I + G4V4 - G4VO = 0\]
% \[5\:~-G4V4 + (G4+GO)VO = 0\]
% \[6\:~V1 = V_{in}\]
% \[7\:~- \alpha G3V3 + V4 = 0\]
% \[8\:~V2 - Lj \omega I_L - V3 = 0\]
%
% Our unknowns are: X = [I_in, V1, V2, I_L, V3, I, V4, Vo];

V_in = 1;
C_n = 0.00001;
I_n = rand;

X = [];

G = zeros(8);
G(1,:) = [1 G1 -G1 0 0 0 0 0];
G(2,:) = [0 -G1 (G1+G2) 1 0 0 0 0];
G(3,:) = [0 0 0 -1 G3 0 0 0];
G(4,:) = [0 0 0 0 0 -1 G4 -G4];
G(5,:) = [0 0 0 0 0 0 -G4 (G4+Go)];
G(6,:) = [0 1 0 0 0 0 0 0];
G(7,:) = [0 0 0 0 -alpha*G3 0 1 0];
G(8,:) = [0 0 1 0 -1 0 0 0];

C = zeros(8);
C(1,:) = [0 Cap -Cap 0 0 0 0 0];
C(2,:) = [0 -Cap Cap 0 0 0 0 0];
C(3,:) = [0 0 0 0 C_n 0 0 0];
C(8,:) = [0 0 0 -L_induct 0 0 0 0];

F = [0 0 -I_n 0 0 V_in 0 0];

%%
%
% We can once again use a FDTD numerical approach to simulate the circuit
% using a Gaussian excitation. The thermal noise can be simulated using a
% random numbers from a Gaussian distribution witha a mean of 0.001A.
%

Circ2Sim

%%
%
% The plot in Figure 19 is similar to Figure 16, however, the output signal
% is noisier due to the inclusion of In. The FFT plot in Figure 20, on the
% other hand, still appears to be the sam as the FFT plot in Figure 18 (for
% the Gaussian excitations without noise). This is likely because the
% computation for the FFT is able to filter out the noise in the signal to
% isolate the actual ferquency of the signal.
%

%%
%
% Next, we may vary Cn and observe the bandwidth to see how it affects the
% frequency response of the circuit.
%

figNum = 21;

for C_n = [1e-15 1e-6 1]
    
    C(3,:) = [0 0 0 0 C_n 0 0 0];
    
    VaryCn
    
    figNum = figNum + 1;
    
end

%%
%
% Figures 21 through 23, above, depict the input and output amplitude
% spectrum of the circuit as Cn changes. In Figure 21, \(C_n = 1 \times
% 10^{-15}F\), in Figure 22, \(C_n = 1 \times 10^{-6}F\) and in Figure 23,
% \(C_n = 1F\). Though these capacitance values are very different, there
% is not much difference in the amplitude spectrum with the only exception
% being the 0Hz response when \(C_n = 1F\).
%

%%
%
% Next, we can also study the effect of dt on the simulation. Note, we do
% have to make some adjustments to the code to resonably scale the input
% signal to changes to dt.
%

C_n = 1e-5;
C(3,:) = [0 0 0 0 C_n 0 0 0];

for dt = [1e-4 1e-3 1e-2]
    
    varyTimeStep
    
    figNum = figNum + 1;
    
end

%%
%
% As mentioned earlier, dt effects the sampling frequency used to measure
% the input and output signals. Because a smaller dt means the signal is
% sampled more frequently, more noise is introduced to the output waveform,
% whereas at higher a dt there is more time between each sample creating a
% smoother looking curve as less noise is picked up.

%% Introducing a Non-Linear Transconductance
%
% 
%
