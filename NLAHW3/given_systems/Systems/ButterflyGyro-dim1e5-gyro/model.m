[M, rows, cols, entries] = mmread('gyro.M');
[B, rows, cols, entries] = mmread('gyro.B'); %single input !!
[K, rows, cols, entries] = mmread('gyro.K');
[C, rows, cols, entries] = mmread('gyro.C'); %12 outputs



beta = 1e-6;
freq = 10.^[1:0.1:4];  
%bode_from_function(fun,b,c,s)
%FRF(j,:) = c'*(fun(s(j))\b);
%My TF= c^T* B / (Ms^2+beta*Ks+K)
H = @(s) s^2 * M + beta * K * s + K;

s = 2*pi*1i*freq;
numOutputs = size(C, 1);
%----DQPA---%
%[lambda, x, y] = QDPA(M, C, K, b, c, s0, epsilon)
%s0=1e2;
%epsilon=1e-10;
%[lambda, x, y] = QDPA(M, beta*K, K, B, C(1,:)', s0, epsilon);
%pole = lambda; 
%residue = (y' * B) * (B' * x);
%s_tf = tf('s');
%H_dqpa = residue / (s_tf - pole);
%figure;
%bode(H_dqpa);

%-------1I-1O------%
FRF = bode_from_function(H, B, C(1, :)', s)';
[pks, locs] = findpeaks(20*log10(abs(FRF)));
%-----SAQDPA----%
options = struct();
options.nwanted = 10; 
options.tol = 1e-6; 
options.displ = 1; 
options.selection_strategy = 'LM'; 
options.kmin = 1; 
options.kmax = 2*options.nwanted;
options.maxrestart = 100; 
options.use_lu = 1; 
options.dpa_bordered = 1; 


d = 0; 
%s0 = ones(options.nwanted, 1) + 1i * ones(options.nwanted, 1);
s0=10*ones(10,1);
%saqdpa(A, E, M, b, c, d, s0, options)
[poles, residues, leftev, rightev, nr_solves] = saqdpa(K, beta*K, M, B, C(1, :)', 0, s0, options);

% Construct the transfer function

H_saqdpa = zeros(size(s)); 
for i = 1:length(poles)
    H_saqdpa = H_saqdpa + residues(i) ./ (s-poles(i)); %!second order  
end
%-------QIRAK----%
tol=1e-6;
s0=10*ones(10,1);
%s0 = 2 * pi * 1i * freq(round(linspace(1, length(freq), 5)));
%[Mr, Dr, Kr, Br, Cr, sigma, Vr] = qirka(M, D, K, B, C, sigma0, tol)
[Mr, Dr, Kr, Br, Cr, sigma, Vr] = qirka(M, beta*K, K, B, C(1, :)', s0, tol);
H_r = @(s) s^2 * Mr + beta * Kr * s + Kr;
FRF_irka = bode_from_function(H_r, Br, Cr, s)';
figure;
subplot(2,1,1);
semilogx(freq, 20*log10(abs(FRF)));
title('Bode Plot - Magnitude');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
hold on;
subplot(2,1,1);
semilogx(freq, 20*log10(abs(H_saqdpa)), 'ro');
legend('Original', 'Reduced');
hold off;

subplot(2,1,2);
semilogx(freq, angle(FRF)*180/pi);
title('Bode Plot - Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (degrees)');
hold on;
semilogx(freq, angle(H_saqdpa)*180/pi, 'ro');
legend('Original', 'Reduced');
hold off;
