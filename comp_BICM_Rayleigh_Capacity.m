clear, clc
% computes BICM capacity of PAM over fast Rayleigh fading channel

SNRdB = [0:0.5:40]; Nw = 1e5; % 1e4:8min,  1e5: 40min, 1e6: 
ldM = 6;

[A, A0, A1] = pam_gray(ldM);
Nsnr = length(SNRdB);  M = pow2(ldM);
Cbicm = zeros(Nsnr,1);
v = 1/sqrt(2) * randn(1,Nw);
h = abs(1/sqrt(2) * [1 1i]*randn(2,Nw)); 

tic
for ns = 1:Nsnr
    N0 = 10^(-SNRdB(ns)/10);
    w = sqrt(N0) * v;

    B = 0;
    for m = 1:ldM
        for j = 1:M/2
            B1 = 0;
            for i = 1:M
                B1 = B1 + exp( -((h*(A0(m,j) - A(i)) + w).^2)/N0 );
            end
            B2 = 0;
            for i = 1:M/2
                B2 = B2 + exp( -((h*(A0(m,j) - A0(m,i)) + w).^2)/N0 );
            end
            B = B + log2(B1./B2);

            B1 = 0;
            for i = 1:M
                B1 = B1 + exp( -((h*(A1(m,j) - A(i)) + w).^2)/N0 );
            end
            B2 = 0;
            for i = 1:M/2
                B2 = B2 + exp( -((h*(A1(m,j) - A1(m,i)) + w).^2)/N0 );
            end
            B = B + log2(B1./B2);
        end
    end
    Cbicm(ns) = ldM - sum(B)/(M*Nw);

    disp([datestr(now) ': SNR = ' num2str(SNRdB(ns)) ' dB,  Cbicm = ' num2str(Cbicm(ns))])
end, toc

snr = 10.^(SNRdB'/10);
Crayleigh = 1/(2*log(2)) * exp(1./snr) .* expint(1./snr);

plot(SNRdB, [Crayleigh Cbicm])

save(['capacity_pam_bicm_rayleigh' int2str(ldM) '.mat'], 'SNRdB', 'Cbicm', 'Nw')

