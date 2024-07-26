clear, clc
% computes BICM capacity of PAM over AWGN channel

SNRdB = [0:0.5:45]; Nw = 1e7; % 22h
%Nw = 1e6;
ldM = 6;

[A, A0, A1] = pam_gray(ldM);
Nsnr = length(SNRdB);  M = pow2(ldM);
Cbicm = zeros(Nsnr,1);
v = 1/sqrt(2) * randn(1,Nw);

tic
for ns = 1:Nsnr
    N0 = 10^(-SNRdB(ns)/10);
    a = A/sqrt(N0); a0 = A0/sqrt(N0); a1 = A1/sqrt(N0);
    
    B = 0;
    for m = 1:ldM
        for j = 1:M/2
            B1 = 0;
            for i = 1:M
                B1 = B1 + exp(- (a0(m,j) - a(i) + v).^2);
            end
            B2 = 0;
            for i = 1:M/2
                B2 = B2 + exp(- (a0(m,j) - a0(m,i) + v).^2);
            end
            B = B + log2(B1./B2);

            B1 = 0;
            for i = 1:M
                B1 = B1 + exp(- (a1(m,j) - a(i) + v).^2);
            end
            B2 = 0;
            for i = 1:M/2
                B2 = B2 + exp(- (a1(m,j) - a1(m,i) + v).^2);
            end
            B = B + log2(B1./B2);
        end
    end
    Cbicm(ns) = ldM - sum(B)/(M*Nw);

    disp([datestr(now) ': SNR = ' num2str(SNRdB(ns)) ' dB,  Cbicm = ' num2str(Cbicm(ns))])
end, toc

Cawgn = 0.5 * log2(1 + 10.^(SNRdB/10))';
plot(SNRdB, [Cawgn Cbicm])

save(['capacity_pam_bicm' int2str(ldM) '.mat'], 'SNRdB', 'Cbicm', 'Nw')