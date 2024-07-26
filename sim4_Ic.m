function sim4_Ic(compmethod, AcdB, Nfr, SNRdB)
%clear, Nfr = 1e4;  SNRdB = 20; compmethod = 'dft'; AcdB = 5;

ldM = 6;
N = 1024; S = 648;
disp(['ld(M) = ' int2str(ldM) ',  Ac = ' num2str(AcdB) ' dB,  N = ' int2str(N) ', SNR = ' int2str(SNRdB)])

% ++ fixed and derived system parameters
A = pam_gray(ldM)';       % PAM alphabet as column vector
Ac = 10^(AcdB/20);
N0 = 10^(-SNRdB/10);

histb = zeros(4096,2,ldM);    % histogram for quantization with W = 12
quantize = @(y)( floor(1 + 4095*qfunc(-y/sqrt(1+N0/2))) );

if strcmp(compmethod,'hardclip') || strcmp(compmethod,'dft')
    alpha = 1/sqrt(1-exp(-Ac^2));   % scaling factor to maintain signal power
    compress = @(x)( alpha * min(abs(x), Ac) .* exp(1i*angle(x)) );
    expand = @(y)( y/alpha );
elseif strcmp(compmethod,'rapp')
    fun = @(x) (2 * x.^3 .* exp(-x.^2))./( (1+(x/Ac).^2));
    alpha1 = 1/sqrt(integral(fun,0,inf));
    compress = @(x)( alpha1 * x./sqrt(1+abs(x/Ac).^2) );
    expand = @(y)( exp(1i*angle(y)) ./ sqrt(abs(alpha1./y).^2 - 1/Ac^2) );
else
    error(['method ' compmethod ' not implemented'])
end


tic
for nfr = 1:Nfr
    u = randi([0 1], 2*ldM, S);    % two PAM symbols per column
    x1 = bit2int(u(1:ldM,:), ldM);  x2 = bit2int(u(ldM+1:2*ldM,:), ldM);
    x = A(x1+1) + 1i*A(x2+1);

    if strcmp(compmethod,'dft'), x = fft(x)/sqrt(S); end

    a = [0; x(1:S/2); zeros(N-S-1,1); x(S/2+1:S)];
    b = compress( N/sqrt(S) * ifft(a) );

    %++ channel
    w = sqrt(N0/2) * randn(S,2)*[1; 1i];
    w = [0; w(1:S/2); zeros(N-S-1,1); w(S/2+1:S)];
    w = N/sqrt(S) * ifft(w);
    y = b + w;

    %++ receiver
    y = expand(y);
    y = sqrt(S)/N * fft(y);
    y = y([2:S/2+1 N-S/2+1:N]);     % subcarrier selection

    if strcmp(compmethod,'dft'), y = sqrt(S)*ifft(y); end

    z1 = quantize(real(y));  z2 = quantize(imag(y));

    for s = 1:S
        for m = 1:ldM
            histb(z1(s), u(m,s)+1, m)     = histb(z1(s), u(m,s)+1, m) + 1;
            histb(z2(s), u(m+ldM,s)+1, m) = histb(z2(s), u(m+ldM,s)+1, m) + 1;
        end
    end
end, toc

filename = [compmethod '_Ac' int2str(AcdB) '.mat'];
ii = SNRdB + 1;
if exist(filename) ~= 2
    res.SNRdB = [0:45];
    res.Nfr = zeros(46,1); res.Nfr(ii) = Nfr;
    res.hist = zeros(4096*2*ldM, 46);  res.hist(:,ii) = histb(:);
else
    load(filename)
    res.Nfr(ii)    = res.Nfr(ii) + Nfr;
    res.hist(:,ii) = res.hist(:,ii) + histb(:);
end 
save(filename, 'res')
