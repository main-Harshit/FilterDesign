close all;
clearvars;
clc;

% Specs
m = 72;
q_m = floor(0.1*m);
r_m = m - 10*q_m;
BL = 25+1.7*q_m + 6.1*r_m;
BH = BL + 20;
trans_bw = 4*10^3;

% Band Edge Specifications
fs1 = BL*10^3-trans_bw;
fp1 = BL*10^3;
fp2 = BH*10^3;
fs2 = BH*10^3+trans_bw;
f_samp = 330e3;
ws1=fs1*2*pi/f_samp;
wp1=fp1*2*pi/f_samp;
wp2=fp2*2*pi/f_samp;
ws2=fs2*2*pi/f_samp;
delta = 0.15;

% Paramaters of the Kaiser Window
A = -20*log10(delta);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end
d_omega_by_pi = trans_bw*2/f_samp;
d_omega=d_omega_by_pi*pi;
N_min = ceil((A-8) / (2*2.285*(d_omega)));          

% Window length for Kaiser Window
addend = 8;  % N_min gives the minimum, add something to it to get exact
% addend is found such that the Stopband specification are met
n=N_min + addend

% Impulse response of Ideal Bandpass Filter for length n
w1=(ws1+wp1)/2;
w2=(ws2+wp2)/2;
bp_ideal = ideal_lp(w2,2*n+1) - ideal_lp(w1,2*n+1);

% Kaiser Window of length n using Shape Parameter Beta as calculated above the in-built function generates a function centred at n  
kaiser_win_shifted = (kaiser(2*n+1,beta))';    
num_arr = [-n:1:n];
kaiser_win = kaiser_win_shifted(num_arr+n+1);  % centered at 0
FIR_BandPass = bp_ideal .* kaiser_win;         
fvtool(FIR_BandPass,'Analysis','freq');         %frequency response
FIR_BandPass
[~,si]=size(FIR_BandPass);
n_axis = linspace((1-si)/2,(si-1)/2,si);
figure;
stem(n_axis,FIR_BandPass,'filled','MarkerSize',4);
hold on;
title("Time Domain Response")
ylabel("h[n]")
xlabel("n")
grid on;

%magnitude response
[H,f] = freqz(FIR_BandPass,1, 1024*1024,f_samp);
figure;
plot(f,abs(H),'LineWidth',1);
hold on;
title("Magnitude Response")
xlabel("Hz")
ylabel("|H(f)|")
xline(fs1,'--g');
xline(fp1,'--m');
yline(0.15,'r');
xline(fs2,'--g');
xline(fp2,'--m');
yline(1.15,'r');
yline(0.85,'r');
grid
legend('Magnitude Response','Passband edge','Stopband edge','Tolerances','location','northeast')