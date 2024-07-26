function sim2_pamldpc(channeltype, mcs, Nwstep)
% clear, channeltype = 'awgn'; mcs = 3; Nwstep = 1e3;

ldM = [1 1 6 6];  Rc = [3/5 3/4 3/5 3/4];
ldM = ldM(mcs);  Rc = Rc(mcs);  M = pow2(ldM);
targetWER = 200;

if strcmp(channeltype, 'awgn')
    SNRdB = floor(10*log10(pow2(2*ldM*Rc)-1));
    h = 1;
    ray = 0;
elseif strcmp(channeltype, 'rayleigh')
    fun = @(snr)( exp(1./snr).*expint(1./snr)/(2*log(2)) - ldM*Rc );
    x0 = pow2(2*ldM*Rc)-1;
    SNRdB = floor(10*log10( fzero(fun,x0) ));
    ray = 1;
end

if ismember(mcs, [1 3])
    H = dvbs2ldpc(3/5);
    [M1,N] = size(H); K = N-M1;
    SNRstep = 0.1;
elseif ismember(mcs, [2 4])
    K = 486; N = 648;      % number of information and coded bits
    P = [16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
         25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
         25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
          9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
         24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
          2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0];
    H = ldpcQuasiCyclicMatrix(27,P);
    SNRstep = 0.5;
end
encoderObj = ldpcEncoderConfig(H);
decoderObj = ldpcDecoderConfig(H);
Nq = N/ldM;

Nw = 0; ber = 0; wer = 0;
filename = ['mcs' int2str(mcs) '_' channeltype '.mat'];
if exist(filename) ~= 2    
    T = table(SNRdB, Nw, ber, wer);
    row = 1;
else
    load(filename)
    Nrows = size(T,1);
    for row = 1:Nrows
        if T.wer(row) < targetWER
            break
        end
    end
    if row == Nrows && T.wer(row) >= targetWER
        SNRdB = T.SNRdB(row) + SNRstep;
        T = [T; table(SNRdB, Nw, ber , wer)];
        row = row + 1;
    end
end
SNRdB = T.SNRdB(row); Nw = Nwstep;

N0 = 10^(-SNRdB/10);
if ismember(mcs, [1 2])
    map = @(c)(  (2*double(c)-1)/sqrt(2) );
    demap = @(y,h)( 2*sqrt(2)/N0 * h.*y );
elseif ismember(mcs, [3 4])
    [A, A0, A1] = pam_gray(ldM);
    A0 = reshape(A0, 1, ldM, M/2);
    A1 = reshape(A1, 1, ldM, M/2);
    map = @(b)( A(bit2int(b,ldM)+1)' );
    demap = @(y,h)( reshape( log( sum( exp( -(y - h.*A1).^2/N0 ), 3)./sum( exp( -(y - h.*A0).^2/N0 ), 3) )', N, 1) );
end



fprintf('SNR = %g dB: ', SNRdB);
ber = 0; wer = 0; tic
for nw = 1:Nw
    u = randi([0 1], K, 1, 'int8');
    c = ldpcEncode(u, encoderObj);
    x = map(c);
    
    % +++ channel
    if ray,  h = abs(1/sqrt(2) * randn(Nq,2) * [1; 1i]);  end
    y = h.*x + sqrt(N0/2)*randn(Nq,1);

    % +++ receiver
    llr = demap(y,h);
    ud = ldpcDecode(-llr, decoderObj, 200);
  
    ber = ber + sum(u ~= ud);
    wer = wer + any(u ~= ud);
end, t = toc;
fprintf('%d bit errors, %d word errors, %d seconds\n', ber, wer, round(t))

% +++ store results
T.Nw(row)  = T.Nw(row) + Nw;
T.ber(row) = T.ber(row) + ber;
T.wer(row) = T.wer(row) + wer;
save(filename, 'T')

