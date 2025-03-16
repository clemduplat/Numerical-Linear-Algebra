%[rows, cols, entries, rep, field, symm] = mminfo('spiral_inductor_peec.E');
[E, rows, cols, entries] = mmread('spiral_inductor_peec.E');

%[rows, cols, entries, rep, field, symm] = mminfo('spiral_inductor_peec.A');
[A, rows, cols, entries] = mmread('spiral_inductor_peec.A');

%[rows, cols, entries, rep, field, symm] = mminfo('spiral_inductor_peec.B');
[B, rows, cols, entries] = mmread('spiral_inductor_peec.B');

%---Bode From Function-----%
freq = 10.^[1:0.1:10];
s = 2 * pi * 1i * freq;
FRF = bode_from_system(A, E, B, B, s);
%tf_model=tf(FRF);

%----DPA-----%
s0 = 1e2; 
tol = 1e-10; 


[lambda_hat, x, y] = DPA(E, A, B,B, s0, tol);
pole = lambda_hat; 
residue = (y' * B) * (B' * x);
s_tf = tf('s');
H_dpa = residue / (s_tf - pole);

%-----SADPA-----%
s0 = [1e6 * 2 * pi * 1i, 1e8 * 2 * pi * 1i]; 


d = 0;
options=struct();
options.nwanted = 4;     
options.tol = 1e-10;     
options.displ = 1;      
options.strategy = 'LM';  
options.kmin = 1;        
options.kmax = 10;        
options.maxrestarts = 100;
options.use_lu=0;
options.use_bordered=0;
[poles, residues, rightev, leftev, nr_solves, ress] = sadpa(A, E, B, B, d, s0, options);

n = length(poles);

s_tf = tf('s'); 
H_sadpa = 0; 
%figure;
for i = 1:length(poles)
    H_sadpa = H_sadpa + residues(i) / (s_tf - poles(i));
    %bode(H_sadpa,'--');
    %legend(['Contribution of Pole ', num2str(i)]);
    %hold on;
end
%hold off;
%-----IRKA----%
s0=[1e2];
[Ahat, Ehat, bhat, chat, V] = irka(A, E, B, B, s0, tol);
FRF_irka = bode_from_system(Ahat, Ehat, bhat, chat, s);
%----GRKA----%
sigma1=0;
smin=i*10^2;
smax=i*10^8;
scount=2;
[Ahat, Ehat, bhat, chat, i] = grka(A, E, B, B, sigma1, smin, smax, scount, tol);
FRF_grka = bode_from_system(Ahat, Ehat, bhat, chat, s);
%----Bode plot ---%
  
figure;
subplot(2,1,1);
semilogx(freq, 20*log10(abs(FRF))); 
hold on;
semilogx(freq, 20*log10(abs(FRF_grka)),'ro');
%hold on;
%bodemag(H_sadpa,'ro');
hold off;
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
title('Bode Plot - Magnitude');
legend("original","Reduced")
subplot(2,1,2);
semilogx(freq, angle(FRF) * 180/pi); 
hold on;
semilogx(freq, angle(FRF_grka) * 180/pi,'ro');
%bodephase(H_sadpa,'ro');
hold off;

ylabel('Phase (degrees)');
xlabel('Frequency (Hz)');
title('Bode Plot - Phase');
legend("original","Reduced")

%------DPA reconstruction-----%

figure;
bodemag(H_dpa);
title('Bode Plot of the Simplified System DPA');

%-----SADPA reconstruction----%


figure;
bode(H_sadpa);
title('Bode Plot of the Reduced Transfer Function SADPA');

