function T = sim5B_wer(compmethod, Nfr, SNRdB)
%clear, Nfr = 1e2;  SNRdB = 35; compmethod = 'rapp'; 

ldM = 6;  N = 1024; S = 648;  AcdB = 5;

% ++ fixed and derived system parameters
Ac = 10^(AcdB/20);
N0 = 10^(-SNRdB/10);

K = 486; Nc = 648;      % number of information and coded bits
[A, A0, A1] = pam_gray(ldM);  
A = A';  M = pow2(ldM);
A0 = reshape(A0, 1, ldM, M/2);
A1 = reshape(A1, 1, ldM, M/2);
map = @(u)( A(bit2int(u(1:ldM,:),ldM)+1) + 1i*A(bit2int(u(ldM+1:2*ldM,:),ldM)+1) );
demap = @(y)( reshape( log( sum( exp( -(y - A1).^2/N0 ), 3)./sum( exp( -(y - A0).^2/N0 ), 3) )', S*ldM, 1) );
llr = zeros(2*ldM,S);

quantize = @(y)( floor(1 + 4095*qfunc(-y/sqrt(1+N0/2))) );

if strcmp(compmethod,'noclip')
    compress = @(x)(x);
    expand = @(y)(y);
elseif strcmp(compmethod,'hardclip') || strcmp(compmethod,'dft')
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

P = [16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
     25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
     25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
      9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
     24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
      2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0];
H = ldpcQuasiCyclicMatrix(27,P);
encoderObj = ldpcEncoderConfig(H);
decoderObj = ldpcDecoderConfig(H);

if ~strcmp(compmethod,'noclip')
    load([compmethod '_Ac' int2str(AcdB) '.mat'])
    [~,i] = min( abs(res.SNRdB - SNRdB));
    hist = reshape(res.hist(:,i), 4096, 2, ldM);
end


ber = 0; wer = 0; tic
for nfr = 1:Nfr
    u = randi([0 1], K, 2*ldM, 'int8');
    c = ldpcEncode(u, encoderObj);
    c = reshape(c, 2*ldM, S);
    x = map(c);

    if strcmp(compmethod,'dft'), x = fft(x)/sqrt(S); end

    a = [0; x(1:S/2); zeros(N-S-1,1); x(S/2+1:S)];
    b = compress( N/sqrt(S) * ifft(a) );

    %++ channel
    w = sqrt(N0/2) * randn(S,2)*[1; 1i];
    w = [0; w(1:S/2); zeros(N-S-1,1); w(S/2+1:S)];
    w = N/sqrt(S) * ifft(w);
    y = b + w;

    %++ receiver
    y = sqrt(S)/N * fft(expand(y));
    y = y([2:S/2+1 N-S/2+1:N]);     % subcarrier selection

    if strcmp(compmethod,'dft'), y = sqrt(S)*ifft(y); end

    if strcmp(compmethod,'noclip')
       llr(1:ldM,:) = reshape(demap(real(y)), ldM, S);
       llr(ldM+1:2*ldM,:) = reshape(demap(imag(y)), ldM, S);
    else
        z1 = quantize(real(y));  z2 = quantize(imag(y));
        for m = 1:ldM
            llr(m,:)     = log( hist(z1,2,m)./hist(z1,1,m) );
            llr(m+ldM,:) = log( hist(z2,2,m)./hist(z2,1,m) );
        end
    end
    llr2 = reshape(llr(:),S, 2*ldM);
    llr2(isnan(llr2)) = 0;
    llr2(llr2==Inf) = 50;
    llr2(llr2==-Inf) = -50;
    ud = ldpcDecode(-llr2, decoderObj, 200);

    ber = ber + sum(u(:) ~= ud(:));
    wer = wer + sum(any(u ~= ud));
end, toc

disp(['SNR = ' num2str(SNRdB) ' dB:  ber = ' int2str(ber) ', wer = ' int2str(wer)])

filename = ['werB_' compmethod '.mat'];
if exist(filename) ~= 2
    T = table(SNRdB, Nfr, ber, wer);
else
    load(filename)
    row  = find(abs(T.SNRdB - SNRdB) < 0.01);  % maximum resolution 0.1 dB
    if isempty(row)                            % new SNR value
        T = [T; table(SNRdB, Nfr, ber, wer)];
        T = sortrows(T);
    else
        T.Nfr(row) = T.Nfr(row) + Nfr;
        T.ber(row) = T.ber(row) + ber;       
        T.wer(row) = T.wer(row) + wer;
    end
end
save(filename, 'T');

