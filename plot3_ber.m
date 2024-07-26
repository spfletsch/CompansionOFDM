clear, close all
figure_settings
% plot uncoded BER over SNR for different clippling levels

ma = {'+', '.', 'x', 's'};
ms = [4 4 4 2];
co = colororder;


figure('Units','centimeters', 'Position',[80 20 18 10])
N = [64 256 1024 8192];
SNRdB = [10:40];

for i = 1:length(N)
    load(['ber_hardclip_N' int2str(N(i)) '.mat'])
    Pb = res.ber/(res.Nfr*2*res.ldM*res.S);
    semilogy(SNRdB, Pb, 'Marker',ma{i}, 'MarkerFaceColor',co(i,:), 'MarkerSize',ms(i));
    if i==1, hold on, end
end

SNRdB = linspace(10,30);  snr = 10.^(SNRdB/10)';
M = pow2(res.ldM);
Pb = zeros(100,1);
for k = 1:res.ldM
  for i = 0:(1-pow2(-k))*M-1
    Pb = Pb + 2/(res.ldM*M) * (-1)^floor(i*pow2(k-1)/M) ...
             * (pow2(k-1) - floor(0i*pow2(k-1)/M + 0.5) ) ...
             * qfunc( (2*i+1)*sqrt(3*snr/(M^2-1)));
  end
end

semilogy(SNRdB, Pb, 'Color', orange, 'LineStyle','--');
xlabel('$$10 \, \mathrm{lg}( \mathrm{SNR} )$$ in dB'), ylabel('BER')
legend('$$N=64$$', '$$N=256$$', '$$N=1024$$', '$$N=8192$$', 'no clipping')
title('$\mathrm{ld} M = 4, A_\mathrm{dB} = 6\,\mathrm{dB}$')
grid on

exportgraphics(gcf, 'Figure-5-ber.pdf',  'ContentType','vector')


