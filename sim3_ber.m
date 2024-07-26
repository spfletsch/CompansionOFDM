function sim3_ber(N, Nfr)
% Nfr = 1e5; N = 1024;

S = 2 * round( 0.4*N);    % ca. 80% of N, even number
ldM = 4; AcdB = 6;
disp(['ld(M) = ' int2str(ldM) ',  Ac = ' num2str(AcdB) ' dB,  N = ' int2str(N)])

% ++ fixed and derived system parameters
SNRdB = [10:40];
Ns = length(SNRdB);
A = pam_gray(ldM)';       % PAM alphabet as column vector
Ac = 10^(AcdB/20);

alpha = 1/sqrt(1-exp(-Ac^2));   % scaling factor to maintain signal power
clip = @(x)( alpha * min(abs(x), Ac) .* exp(1i*angle(x)) );

ber = zeros(Ns,1); tic
for ns = 1:Ns
    N0 = 10^(-SNRdB(ns)/10);
    for nfr = 1:Nfr
        u = randi([0 1], 2*ldM, S);    % two PAM symbols per column
        x1 = bit2int(u(1:ldM,:), ldM);  x2 = bit2int(u(ldM+1:2*ldM,:), ldM);
        x = A(x1+1) + 1i*A(x2+1);

        a = [0; x(1:S/2); zeros(N-S-1,1); x(S/2+1:S)];
        b = clip( N/sqrt(S) * ifft(a) );

        %++ channel
        w = sqrt(N0/2) * randn(S,2)*[1; 1i];
        w = [0; w(1:S/2); zeros(N-S-1,1); w(S/2+1:S)];
        w = N/sqrt(S) * ifft(w);
        y = b + w;

        %++ receiver
        y = sqrt(S)/N * fft(y);
        y = y([2:S/2+1 N-S/2+1:N]);     % subcarrier selection

        u2 = [pam_harddemap(real(y)', ldM); pam_harddemap(imag(y)', ldM)];
        ber(ns) = ber(ns) + sum(u(:) ~= u2(:));
    end
    disp(['SNR = ' int2str(SNRdB(ns)) ', ' int2str(ber(ns)) ' bit errors'])
end, toc

Pb = ber/(Nfr*2*ldM*S);
semilogy(SNRdB, Pb, 'Marker','+')
grid

filename = ['ber_hardclip_N' int2str(N) '.mat'];
if exist(filename) ~= 2
    res.ldM = ldM;  res.S = S;  res.Nfr = Nfr;  res.ber = ber;
else
    load(filename)
    res.Nfr = res.Nfr + Nfr;
    res.ber = res.ber + ber;
end
save(filename, 'res')






