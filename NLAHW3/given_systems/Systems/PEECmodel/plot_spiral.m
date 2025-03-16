%[rows, cols, entries, rep, field, symm] = mminfo('spiral_inductor_peec.E');
[spiral_inductor_peec.E, rows, cols, entries] = mmread('spiral_inductor_peec.E');

%[rows, cols, entries, rep, field, symm] = mminfo('spiral_inductor_peec.A');
[spiral_inductor_peec.A, rows, cols, entries] = mmread('spiral_inductor_peec.A');

%[rows, cols, entries, rep, field, symm] = mminfo('spiral_inductor_peec.B');
[spiral_inductor_peec.B, rows, cols, entries] = mmread('spiral_inductor_peec.B');

%inductance
L_spiral_peec = spiral_inductor_peec.E;
%resistance
R_spiral_peec = -spiral_inductor_peec.A;
%input
B_spiral_peec = spiral_inductor_peec.B;

freq = 10.^[1:0.1:10];

RinvB_spiral_peec = R_spiral_peec\B_spiral_peec;
RinvL_spiral_peec = R_spiral_peec\L_spiral_peec;

[V_spiral_peec,D_spiral_peec] = eig(RinvL_spiral_peec);

BtV_spiral_peec = transpose(B_spiral_peec)*V_spiral_peec;
VinvRinvB_spiral_peec = V_spiral_peec\RinvB_spiral_peec;

y_spiral_peec = zeros(size(freq));

%contribution of all the modes -> construction of N
for j = 1:size(freq,2)
  y_spiral_peec(1,j) = sum(transpose(BtV_spiral_peec)./(1+sqrt(-1)*2*pi*freq(j)*diag(D_spiral_peec)).*VinvRinvB_spiral_peec);
end

z_spiral_peec = 1./y_spiral_peec;

resis_spiral_peec = real(z_spiral_peec);
induc_spiral_peec = imag(z_spiral_peec)./(2*pi*freq);
%----------------------------------------------------------%
% m is PRIMA reduction order.  I found 50 is sufficient

m = 50;

fprintf('The PRIMA reduction order is chosen to be %i.\n',m);
%prima(C, G, B, m, verbose)
[Cr, Gr, X] = prima(L_spiral_peec, -R_spiral_peec, B_spiral_peec, m, 1);

L_prima = Cr(1:m,1:m);
R_prima = -Gr(1:m,1:m);
B_prima = transpose(X(:,1:m))*B_spiral_peec;

[V,D]=eig(Cr(1:m,1:m));
I = find(abs(diag(D))>=1e-20);
L_prima = D(I,I);
R_prima = -transpose(V(:,I))*Gr(1:m,1:m)*V(:,I);
B_prima = transpose(V(:,I))*transpose(X(:,1:m))*B_spiral_peec;

RinvB_prima = R_prima\B_prima;
RinvL_prima = R_prima\L_prima;

[V_prima,D_prima] = eig(RinvL_prima);

BtV_prima = transpose(B_prima)*V_prima;
VinvRinvB_prima = V_prima\RinvB_prima;

y_prima = zeros(size(freq));

for j = 1:size(freq,2)
  y_prima(1,j) = sum(transpose(BtV_prima)./(1+sqrt(-1)*2*pi*freq(j)*diag(D_prima)).*VinvRinvB_prima);
end

z_prima = 1./y_prima;

resis_prima = real(z_prima);
induc_prima = imag(z_prima)./(2*pi*freq);

Lhalf = sqrtm(L_prima);
Lhalfinv = Lhalf\eye(size(Lhalf));

A_symm_ss = -Lhalfinv*R_prima*Lhalfinv;
B_symm_ss = Lhalfinv*B_prima;
C_symm_ss = transpose(B_prima)*Lhalfinv;

[V_symm_ss,D_symm_ss]=eig(A_symm_ss);
VinvB_symm_ss = V_symm_ss\B_symm_ss;
CV_symm_ss = C_symm_ss*V_symm_ss;

y_symm_ss = zeros(size(freq));

for j = 1:size(freq,2)
  y_symm_ss(1,j) = sum(transpose(CV_symm_ss)./(sqrt(-1)*2*pi*freq(j)-diag(D_symm_ss)).*VinvB_symm_ss);
end

z_symm_ss = 1./y_symm_ss;

resis_symm_ss = real(z_symm_ss);
induc_symm_ss = imag(z_symm_ss)./(2*pi*freq);


%original peec models vs reduced models in terms or resistance and
%inducatnce across the frequency range
figure(1);
subplot(2,2,1); semilogx(freq,resis_spiral_peec,'-',freq,resis_prima,'o'); 
title('Resistance'); legend(cat(2,'spiral-peec-',num2str(entries)),cat(2,'PRIMA-reduc-',num2str(m)));
subplot(2,2,2); semilogx(freq,induc_spiral_peec,'-',freq,induc_prima,'o'); 
title('Inductance'); legend(cat(2,'spiral-peec-',num2str(entries)),cat(2,'PRIMA-reduc-',num2str(m)));
subplot(2,2,3); semilogx(freq,resis_spiral_peec,'-',freq,resis_symm_ss,'o'); 
title('Resistance'); legend(cat(2,'spiral-peec-',num2str(entries)),cat(2,'symm state-space-',num2str(m)));
subplot(2,2,4); semilogx(freq,induc_spiral_peec,'-',freq,induc_symm_ss,'o'); 
title('Inductance'); legend(cat(2,'spiral-peec-',num2str(entries)),cat(2,'symm state-space-',num2str(m)));

%---------------------------------------------------------------------------------------------------------%


