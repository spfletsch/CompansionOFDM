clear, close all
figure_settings

ls = {'-', '--', '-', '--'};
ma = {'o', 's', '^', 'd'};

figure('Units','centimeters', 'Position',[80 40 18 10])

% +++ AWGN channel
for mcs = 1:4
    load(['mcs' int2str(mcs) '_awgn.mat'])
    semilogy(T.SNRdB, T.wer./T.Nw, 'Color','blue', 'LineStyle', ls{mcs}, 'Marker', ma{mcs}, 'MarkerFaceColor','b')
    if mcs == 1, hold on, end
end
% +++ Rayleigh fading channel
for mcs = 1:4
    load(['mcs' int2str(mcs) '_rayleigh.mat'])
    semilogy(T.SNRdB, T.wer./T.Nw, 'Color','r', 'LineStyle', ls{mcs}, 'Marker', ma{mcs})
end
xlabel('$$10 \, \mathrm{lg}( \mathrm{SNR} )$$ in dB'), ylabel('WER'), grid
axis([0 40 4e-6 1])
legend('MCS 1, AWGN', 'MCS 2, AWGN', 'MCS 3, AWGN', 'MCS 4, AWGN', ...
    'MCS 1, Rayleigh', 'MCS 2, Rayleigh', 'MCS 3, Rayleigh', 'MCS 4, Rayleigh', 'Location','best')
legend('boxoff')

exportgraphics(gcf, 'Figure-2-wer-snr.pdf',  'ContentType','vector')


% +++ plot over Mutual Information instead of SNR
figure('Units','centimeters', 'Position',[80 20 18 10])

% BI-AWGN channel
for mcs = 1:2
    load(['mcs' int2str(mcs) '_awgn.mat'])
    mi = Jfun(2*10.^(T.SNRdB/20));
     semilogy(mi, T.wer./T.Nw, 'Color','b', 'MarkerFaceColor','b', 'LineStyle',ls{mcs}, 'Marker',ma{mcs})
    if mcs == 1, hold on, end
end

% 64-PAM over AWGN
load('capacity_pam_bicm6.mat')
for mcs = 3:4
    load(['mcs' int2str(mcs) '_awgn.mat'])
    mi = interp1(SNRdB, Cbicm/6, T.SNRdB);
    semilogy(mi, T.wer./T.Nw, 'Color','b', 'MarkerFaceColor','b', 'LineStyle',ls{mcs}, 'Marker',ma{mcs})
end

% BI Rayleigh channel
beta = @(x)(0.5*(psi((x+1)/2) - psi(x/2)));
for mcs = 1:2
    load(['mcs' int2str(mcs) '_rayleigh.mat'])
    snr = 10.^(T.SNRdB/10);
    mi = 1./(log(2) * sqrt(1+2./snr)) .* beta( (1 + sqrt(1+2./snr))/2 );
    semilogy(mi, T.wer./T.Nw, 'Color','r', 'LineStyle',ls{mcs}, 'Marker',ma{mcs})
end

% 64-PAM over Rayleigh channel
load('capacity_pam_bicm_rayleigh6.mat')
for mcs = 3:4
    load(['mcs' int2str(mcs) '_rayleigh.mat'])
    mi = interp1(SNRdB, Cbicm/6, T.SNRdB);
    semilogy(mi, T.wer./T.Nw, 'Color','r', 'LineStyle',ls{mcs}, 'Marker',ma{mcs})
end
  
text(0.67, 5e-5, ['empty markers: Rayleigh channel   '; 'filled markers: AWGN channel      ';...
                  'solid lines: $R_\mathrm{c} = 3/5$ '; 'dashed lines: $R_\mathrm{c} = 3/4$'], ... 
    'BackgroundColor','w')
legend('MCS 1, AWGN', 'MCS 2, AWGN', 'MCS 3, AWGN', 'MCS 4, AWGN', ...
    'MCS 1, Rayleigh', 'MCS 2, Rayleigh', 'MCS 3, Rayleigh', 'MCS 4, Rayleigh', 'Location','SW')
xlabel('Mutual Information per Coded Bit $I_\mathrm{c}$')
ylabel('WER'), grid
ylim([4e-6 1])

exportgraphics(gcf, 'Figure-3-wer-Ic.pdf',  'ContentType','vector')


% +++ BER and WER curves over Ic for each code

% code A: (64800, 38880) code from DVB-S2
K = 38880;
Ic = cell(4,1);
Pw = cell(4,1);  Pb = cell(4,1);

load('mcs1_awgn.mat')
Ic{1} = Jfun(2*10.^(T.SNRdB/20)); Pw{1} = T.wer./T.Nw; Pb{1} = T.ber./(K*T.Nw);

load('mcs1_rayleigh.mat')
snr = 10.^(T.SNRdB/10);
Ic{2} = 1./(log(2) * sqrt(1+2./snr)) .* beta( (1 + sqrt(1+2./snr))/2 );
Pw{2} = T.wer./T.Nw;  Pb{2} = T.ber./(K*T.Nw);

load('mcs3_awgn.mat')
load('capacity_pam_bicm6.mat')
Ic{3} = interp1(SNRdB, Cbicm, T.SNRdB)/6;
Pw{3} = T.wer./T.Nw;  Pb{3} = T.ber./(K*T.Nw);

load('mcs3_rayleigh.mat')
load('capacity_pam_bicm_rayleigh6.mat')
Ic{4} = interp1(SNRdB, Cbicm, T.SNRdB)/6;
Pw{4} = T.wer./T.Nw;  Pb{4} = T.ber./(K*T.Nw);

Icmin = min([Ic{1}; Ic{2}; Ic{3}; Ic{4}]);
Icmax = max([Ic{1}; Ic{2}; Ic{3}; Ic{4}]);
IcA = linspace(Icmin, Icmax);
berA = zeros(100,1);  werA = zeros(100,1);
for i = 1:100
    berA(i) = 10^( (interp1(Ic{1}, log10(Pb{1}), IcA(i)) + interp1(Ic{2}, log10(Pb{2}), IcA(i)) ...
                  + interp1(Ic{3}, log10(Pb{3}), IcA(i)) + interp1(Ic{4}, log10(Pb{4}), IcA(i)) )/4 );
    werA(i) = 10^( (interp1(Ic{1}, log10(Pw{1}), IcA(i)) + interp1(Ic{2}, log10(Pw{2}), IcA(i)) ...
                  + interp1(Ic{3}, log10(Pw{3}), IcA(i)) + interp1(Ic{4}, log10(Pw{4}), IcA(i)) )/4 );
end
i = max(find(berA>0)); berA(i+1:end) = 0;
i = max(find(werA>0)); werA(i+1:end) = 0;

figure('Units','normalized', 'Position',[0.7 0.7 0.3 0.3])
for i = 1:4
    semilogy(Ic{i}, Pw{i}, 'Marker','.'); hold on
    semilogy(IcA, werA, 'Color','k', 'Linewidth',1)
end

% code B: (648, 486) code from WLAN
K = 486;
Ic = cell(4,1);
Pw = cell(4,1);  Pb = cell(4,1);

load('mcs2_awgn.mat')
snr = 10.^(T.SNRdB/10);
Ic{1} = Jfun(2*sqrt(snr)); Pw{1} = T.wer./T.Nw; Pb{1} = T.ber./(K*T.Nw);

load('mcs2_rayleigh.mat')
snr = 10.^(T.SNRdB/10);
Ic{2} = 1./(log(2) * sqrt(1+2./snr)) .* beta( (1 + sqrt(1+2./snr))/2 );
Pw{2} = T.wer./T.Nw;  Pb{2} = T.ber./(K*T.Nw);

load('mcs4_awgn.mat')
load('capacity_pam_bicm6.mat')
Ic{3} = interp1(SNRdB, Cbicm, T.SNRdB)/6;
Pw{3} = T.wer./T.Nw;  Pb{3} = T.ber./(K*T.Nw);

load('mcs4_rayleigh.mat')
load('capacity_pam_bicm_rayleigh6.mat')
Ic{4} = interp1(SNRdB, Cbicm, T.SNRdB)/6;
Pw{4} = T.wer./T.Nw;  Pb{4} = T.ber./(K*T.Nw);

Icmin = min([Ic{1}; Ic{2}; Ic{3}; Ic{4}]);
Icmax = max([Ic{1}; Ic{2}; Ic{3}; Ic{4}]);
IcB = linspace(Icmin, Icmax);
berB = zeros(100,1);  werB = zeros(100,1);
for i = 1:100
    berB(i) = 10^( (interp1(Ic{1}, log10(Pb{1}), IcB(i)) + interp1(Ic{2}, log10(Pb{2}), IcB(i)) ...
                  + interp1(Ic{3}, log10(Pb{3}), IcB(i)) + interp1(Ic{4}, log10(Pb{4}), IcB(i)) )/4 );
    werB(i) = 10^( (interp1(Ic{1}, log10(Pw{1}), IcB(i)) + interp1(Ic{2}, log10(Pw{2}), IcB(i)) ...
                  + interp1(Ic{3}, log10(Pw{3}), IcB(i)) + interp1(Ic{4}, log10(Pw{4}), IcB(i)) )/4 );
end
i = max(find(berB>0)); berB(i+1:end) = 0;
i = max(find(werB>0)); werB(i+1:end) = 0;

for i = 1:4
    semilogy(Ic{i}, Pw{i}, 'Marker','.'); hold on
    semilogy(IcB, werB, 'Color','k', 'Linewidth',1)
end
save('berwerIc.mat', 'IcA','berA','werA', 'IcB','berB','werB')




