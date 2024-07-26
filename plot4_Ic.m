clear, close all
figure_settings
% plot mutual information over SNR for different clippling levels

Ac = [4 5 6 7 8];
ldM = 6;

figure('Units','centimeters', 'Position',[80 30 18 10])

SNRdB = linspace(0,45); snr = 10.^(SNRdB/10);
plot(SNRdB, 0.5 * log2(1 + snr), 'LineStyle','--', 'Color','k'); hold on

load('capacity_pam_bicm6.mat')
plot(SNRdB, Cbicm, 'Color',lightblue, 'Linewidth',1);

text(25, 4.8, '$$\frac12 \mathrm{ld}( 1 + \mathrm{SNR} )$$', 'rotation',35)
text(40,6.3, '$$C_\mathrm{BICM}$$', 'Color',lightblue)

co = [0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19; 0.64 0.08 0.18; 1 0.41 0.16; 0 1 0];
lw = [0.5 1 0.5 0.5 1 0.5];

SNRdB = [0:45];
for a = 1:length(Ac)
    load(['hardclip_Ac' int2str(Ac(a)) '.mat'])
    Iq = zeros(46,1);
    for i = 1:46
        hist = reshape(res.hist(:,i), 4096, 2, ldM);
        Iq(i) = mutual_information(hist);
    end

    plot(SNRdB, Iq, 'Color',co(a,:), 'Linewidth',lw(a))
    text(41, Iq(41)-0.1, ['$$' int2str(Ac(a)) '\,\mathrm{dB}$$'], 'Backgroundcolor','w', 'Color',co(a,:))
end

xlabel('$$10 \, \mathrm{lg}( \mathrm{SNR} )$$ in dB')
ylabel('Mutual information $I_\mathrm{q}$ [bpcu]')
grid
title(['$$\mathrm{ld} M = ' int2str(ldM) ', \; N = 1024, \; S = 648$$'])
axis([0 45 0 7])

exportgraphics(gcf, 'Figure-6-Iq-Ac.pdf',  'ContentType','vector')


%----------------------------------------------------------------
figure('Units','centimeters', 'Position',[80 12 18 10])

SNRdB = linspace(0,45); snr = 10.^(SNRdB/10);
plot(SNRdB, 0.5 * log2(1 + snr), 'LineStyle','--', 'Color','k'); hold on

load('capacity_pam_bicm6.mat')
plot(SNRdB, Cbicm, 'Color',lightblue);

text(25, 4.8, '$$\frac12 \mathrm{ld}( 1 + \mathrm{SNR} )$$', 'rotation',35)

filename{1} = 'rapp_Ac5.mat';
filename{2} = 'dft_Ac5.mat';
filename{3} = 'hardclip_Ac5.mat';
SNRdB = [0:45];

co = [0 0.6 0; 1 0 0; 0.49 0.18 0.56];
lw = [1 0.5 1];


for n = 1:3
  load(filename{n})

  Iq = zeros(46,1);
  for i = 1:46
    hist = reshape(res.hist(:,i), 4096, 2, ldM);
    Iq(i) = mutual_information(hist);
  end

  plot(SNRdB, Iq, 'Linewidth',lw(n), 'Color',co(n,:));
end

xlabel('$$10 \, \mathrm{lg}( \mathrm{SNR} )$$ in dB')
ylabel('Mutual information $I_\mathrm{q}$ [bpcu]')
grid
title(['$$\mathrm{ld} M = ' int2str(ldM) ', \; N = 1024, \; S = 648, \; A_\mathrm{dB} = 5 \,\mathrm{dB}$$'])
axis([0 45 0 7])
legend('AWGN capacity', '$$C_\mathrm{BICM}$$', 'Rapp compansion', 'DFT precoding', 'Hard clipping', ...
    'Location','SE')

exportgraphics(gcf, 'Figure-7-Iq-comp.pdf',  'ContentType','vector')


