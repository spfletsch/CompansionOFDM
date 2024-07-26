clear, close
figure_settings
figure('Units','centimeters', 'Position',[60 30 14 10])

SNRdB = [0:0.25:45]; Cawgn = 0.5 * log2(1 + 10.^(SNRdB/10));
plot(SNRdB, Cawgn, 'Color','m', 'Linewidth',1); hold on

load('capacity_pam_bicm6.mat')
plot(SNRdB, Cbicm, 'Color','b');

load('pam_w7.mat')
plot(SNRdB, Iq, 'Color',[0 0.6 0], 'Linestyle','--')

load('pam_w10.mat')
plot(SNRdB, Iq, 'Color','r', 'Linestyle','none', 'Marker','+', 'Markersize',3)

xlabel('$$10 \, \mathrm{lg}( \mathrm{SNR} )$$ in dB'); 
ylabel('$$C_\mathrm{BICM}, I_\mathrm{q}$$ [bit per channel use]')
grid, axis([0 45 0 7])
legend('$$C_\mathrm{AWGN}$$', '$$C_\mathrm{BICM}$$', '$$I_\mathrm{q}, W = 7$$', ...
    '$$I_\mathrm{q}, W = 10$$', 'Location','NW')

exportgraphics(gcf, 'Figure-1-Cbicm.pdf',  'ContentType','vector')



