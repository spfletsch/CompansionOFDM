clear, close all
figure_settings
% plot simulated error rates for different compression methods

figure('Units','centimeters', 'Position',[80 30 18 10])

load('berwerIc.mat')

% +++ code B
load('werB_noclip.mat')
ldM = 6;  K = 486;
snr1 = T.SNRdB;

%semilogy(snr1, T.ber./(2*ldM*T.Nfr*K), 'Marker','+', 'Color','b'), hold on
semilogy(snr1, T.wer./(2*ldM*T.Nfr), 'Marker','+', 'Color','b'), hold on
load(['capacity_pam_bicm6.mat'])    % without clipping: BICM capacity
Ic = interp1(SNRdB, Cbicm, T.SNRdB)/ldM;    % mutual information per coded bit
%ber2 = 10.^interp1(IcB, log10(berB), Ic);
wer2 = 10.^interp1(IcB, log10(werB), Ic);
%semilogy(snr1, ber2, 'Linestyle',':', 'Color','b')
semilogy(snr1, wer2, 'Linestyle',':', 'Color','b')


% +++ DFT precoding
load('werB_dft.mat')
snr1 = T.SNRdB;
%semilogy(snr1, T.ber./(2*ldM*T.Nfr*K), 'Marker','>', 'Markersize',4, 'Color','r')
semilogy(snr1, T.wer./(2*ldM*T.Nfr), 'Marker','>', 'Markersize',4, 'Color','r')

load('dft_Ac5.mat')    % for mutual information with DFT precoding
Ic1 = zeros(46,1);
for i = 1:46
    hist = reshape(res.hist(:,i), 4096, 2, ldM);
    Ic1(i) = mutual_information(hist)/ldM;
end
Ic = interp1([0:45], Ic1, snr1);
%ber2 = 10.^interp1(IcB, log10(berB), Ic);
wer2 = 10.^interp1(IcB, log10(werB), Ic);
semilogy(snr1, wer2, 'Linestyle',':', 'Color','r')


% +++ Rapp compansion
load('werB_rapp.mat')
snr1 = T.SNRdB;
semilogy(snr1, T.wer./(2*ldM*T.Nfr), 'Marker','.', 'Color',[0 0.6 0])

load('rapp_Ac5.mat')   % for mutual information with Rapp compansion
Ic1 = zeros(46,1);
for i = 1:46
    hist = reshape(res.hist(:,i), 4096, 2, ldM);
    Ic1(i) = mutual_information(hist)/ldM;
end
Ic = interp1([0:45], Ic1, snr1);
wer2 = 10.^interp1(IcB, log10(werB), Ic);
semilogy(snr1, wer2, 'Linestyle',':', 'Color',[0 0.6 0])


% +++ code A +++++++++++++++++++++++++++++++
K = 38880;
load('werA_noclip.mat')
snr1 = T.SNRdB;
h1 = semilogy(snr1, T.wer./T.Nfr, 'Marker','+', 'Color','b'); hold on
load('capacity_pam_bicm6.mat')     % without clipping: BICM capacity
Ic = interp1(SNRdB, Cbicm, T.SNRdB)/ldM;    % mutual information per coded bit
wer2 = 10.^interp1(IcA, log10(werA), Ic);
semilogy(snr1, wer2, 'Linestyle',':', 'Color','b')


% +++ DFT precoding
load('werA_dft.mat')
snr1 = T.SNRdB;
h2 = semilogy(snr1, T.wer./T.Nfr, 'Marker','>', 'MarkerSize',4, 'Color','r');

load('dft_Ac5.mat')    % for mutual information with DFT precoding
Ic1 = zeros(46,1);
for i = 1:46
    hist = reshape(res.hist(:,i), 4096, 2, ldM);
    Ic1(i) = mutual_information(hist)/ldM;
end
Ic = interp1([0:45], Ic1, snr1);
wer2 = 10.^interp1(IcA, log10(werA), Ic);
semilogy(snr1, wer2, 'Linestyle',':', 'Color','r')


% +++ Rapp compansion
load('werA_rapp.mat')
snr1 = T.SNRdB;
h3 = semilogy(snr1, T.wer./T.Nfr, 'Marker','.', 'Color',[0 0.6 0]);

load('rapp_Ac5.mat')   % for mutual information with Rapp compansion
Ic1 = zeros(46,1);
for i = 1:46
    hist = reshape(res.hist(:,i), 4096, 2, ldM);
    Ic1(i) = mutual_information(hist)/ldM;
end
Ic = interp1([0:45], Ic1, snr1);
wer2 = 10.^interp1(IcA, log10(werA), Ic);
semilogy(snr1, wer2, 'Linestyle',':', 'Color',[0 0.6 0])

text(23, 4e-3, 'MCS 3', 'Backgroundcolor','w')
text(29.3, 0.06, 'MCS 4', 'Backgroundcolor','w')
annotation('ellipse', [0.4 0.5 0.08 0.04])
annotation('ellipse', [0.61 0.77 0.03 0.18])

xlabel('SNR [dB]'), ylabel('WER')
axis([18 38 2e-5 1])
legend([h1;h2;h3], 'no clipping', 'DFT precoding', 'Rapp compansion', 'Location','SW')
grid

exportgraphics(gcf, 'Figure-8-coded-ber.pdf',  'ContentType','vector')


