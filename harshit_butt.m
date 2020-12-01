%Butterworth Analog LPF parameters
Wc = 1.071;              %cut-off frequency
N = 7;                  %order 

%poles of Butterworth polynomial of degree 8 in the open CLHP 
p1 = Wc*cos(pi/2 + pi/14) + i*Wc*sin(pi/2 + pi/14);
p2 = Wc*cos(pi/2 + pi/14) - i*Wc*sin(pi/2 + pi/14);
p3 = Wc*cos(pi/2 + pi/14+pi/7) + i*Wc*sin(pi/2 + pi/14+pi/7);
p4 = Wc*cos(pi/2 + pi/14+pi/7) - i*Wc*sin(pi/2 + pi/14+pi/7);
p5 = Wc*cos(pi/2 + pi/14+2*pi/7) + i*Wc*sin(pi/2 + pi/14+2*pi/7);
p6 = Wc*cos(pi/2 + pi/14+2*pi/7) - i*Wc*sin(pi/2 + pi/14+2*pi/7);
p7 = Wc*cos(pi/2 + pi/14+3*pi/7) + i*Wc*sin(pi/2 + pi/14+3*pi/7);
%p8 = Wc*cos(pi/2 + pi/16+3*pi8) - i*Wc*sin(pi/2 + pi/16+3*pi/8);
% 16 -> 14 ; 8-> 7 

%Band Edge speifications
fp1 = 47.5;
fs1 = 51.5;
fs2 = 71.5;
fp2 = 75.5;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 260;         
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

%Parameters for Bandstop Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7],Wc^N);   %TF with poles p1-p8 and numerator Wc^N and no zeroes
                                                        %numerator chosen to make the DC Gain = 1

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bsf(s) = analog_lpf((B*s)/(s*s + W0*W0));        %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog bsf
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete bsf
[nz, dz] = numden(discrete_bsf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
kd = dz(1);                                             %normalisation factor
k = dz(1);    
dz = dz/k;
nz = nz/k;
fvtool(nz,dz,'Analysis','freq');                        %frequency response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,p,~]=tf2zp(nz,dz);
figure;
plot(real(p),imag(p),'rX');
title("Poles of the filter transfer function")
xlabel("Re(z)")
ylabel("Im(z)")
axis equal
grid on
t = linspace(0,2*pi,1000);
hold on
plot(cos(t),sin(t),'b-') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, f_samp);
figure;
plot(f,abs(H),'LineWidth',1);
hold on;
title("Magnitude Response")
xlabel("Hz")
ylabel("|H(f)|")
xline(fs1,'--m');
xline(fp1,'--m');
xline(fp2,'--m');
xline(fs2,'--m');
yline(0.85,'--m');
yline(0.15,'--m');
grid