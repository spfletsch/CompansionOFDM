clear, clc, close
% PAM over AWGN channel: find mutual information after quantization with
% word width W bits

W = 10;         % word width of quantizer in bits
ldM = 6;       % bits per PAM symbol
Np = 1e7;      % number of PAM symbols
Nfr = 10;      % number of frames with Np PAM symbols

SNRdB = 0:45;
A = pam_gray(ldM);  Es = 1/2;  M = pow2(ldM);
Nsnr = length(SNRdB);    % number of SNR values
Iq = zeros(Nsnr,1);      % mutual information

tic
for ns = 1:Nsnr
    histc = zeros(pow2(W),2,ldM);    % histograms for joint probalities (z, cm)
    N0 = 2*Es * 10^(-SNRdB(ns)/10);
    quantize = @(y)( floor(1 + (pow2(W)-1)*qfunc(-y/sqrt(1+N0/2))) );
    
    for nfr = 1:Nfr
        c = randi([0 1], ldM, Np);
        u = bit2int(c, ldM);
        x = A(u+1);

        y = x + sqrt(N0/2)*randn(1,Np);
        z = quantize(y);

        for i = 1:Np
            for m = 1:ldM
                histc(z(i), c(m,i)+1, m) = histc(z(i), c(m,i)+1, m) + 1;
            end
        end
    end

    Iq(ns) = mutual_information(histc);
    disp(['SNR = ' int2str(SNRdB(ns)) ' dB, Iq = ' num2str(Iq(ns))])
end, toc

save(['pam_w' int2str(W) '.mat'], 'SNRdB', 'Iq', 'Np', 'Nfr')
