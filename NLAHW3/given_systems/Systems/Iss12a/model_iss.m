load('iss12a.mat', 'A', 'B', 'C', 'D');
freq = 10.^[-2:0.1:2];  % from 10^-2 to 10^2 
s = 2 * pi * 1i * freq; % s = 2*pi*i*omega
if ~exist('E', 'var')
    E = eye(size(A));
end
sigma0 = 10.^[-1:1:1];
FRF = bode_from_system(A, E, B(:,1), C(1,:)', s);
[pks, locs] = findpeaks(20*log10(abs(FRF)));
%only for one input output here
%-----DPA-----%
%works but not good estimate->logic
%tol=1e-10;
%sigma0=1;
%[lambda_hat, x, y] = DPA(A, E, B(:,1), C(1,:)', sigma0, tol);
%pole = lambda_hat; 
%residue = (y' * B) * (B' * x);
%s = tf('s');
%H = residue / (s - pole);
%figure;
%bode(H);
%title('Bode Plot of the Simplified System DPA');

%-----SADPA---%
tol=1e-10;
s0 =  pks;
options = struct();
options.nwanted = 20; 
options.tol = 1e-6; 
options.displ = 1; 
options.strategy = 'LM'; 
options.kmin = 0; 
options.kmax = 100;
options.maxrestarts = 100; 
options.use_lu=0;
options.use_bordered=0;
%options.yEx_scaling = 1;


[poles, residues, rightev, leftev, nr_solves, ress] = sadpa(A, E, B(:,1), C(1,:)', 0, s0, options);
H_reduced = zeros(size(s));
for i = 1:length(poles)
    H_reduced = H_reduced + (residues(i) ./ (s - poles(i)));
end
[pks, locs] = findpeaks(20*log10(abs(FRF)));
r=10;
s_0 = 10*ones(10,1);

%-----IRKA-----%
%[Ahat, Ehat, bhat, chat, V] = irka(A, E, B(:,1), C(1,:)', s_0, tol);
%FRF_reduced = bode_from_system(Ahat, Ehat, bhat, chat, s);

%------GRKA----%
sigma1=0.1;
smin=i*10^-1;
smax=i*10^1;
scount=50;
[Ahat, Ehat, bhat, chat, i] = grka(A, E, B(:,1), C(1,:)', sigma1, smin, smax, scount, tol);
FRF_reduced = bode_from_system(Ahat, Ehat, bhat, chat, s);
figure;
subplot(2,1,1);
semilogx(freq, 20*log10(abs(FRF)));
title('Bode Plot - Magnitude');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
hold on;
subplot(2,1,1);
semilogx(freq, 20*log10(abs(FRF_reduced)), 'ro');
legend('Original', 'Reduced');
hold off;

subplot(2,1,2);
semilogx(freq, angle(FRF)*180/pi);
title('Bode Plot - Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (degrees)');
hold on;
semilogx(freq, angle(FRF_reduced)*180/pi, 'ro');
legend('Original', 'Reduced');
hold off;