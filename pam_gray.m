function [A, A0, A1] = pam_gray(Rm)
% finds PAM constellation for Rm bits with binary reflected Gray labeling,
% scaled to Es = 1/2
% corresponds to A = sqrt(Es) * real(pammod(0:pow2(Rm)-1, pow2(Rm), 0, 'gray'))
    b = cell(Rm,1);

    b{1} = [0 1];

    for m = 2:Rm
        M = pow2(m);          % number of PAM symbols
        Es = (M^2-1)/3;       % average signal power
        b{m} = zeros(m, M);   % bit matrix: columns correspond to PAM symbols
        b{m}(1,1:M/2) = 0;  b{m}(1,M/2+1:end) = 1;
        b{m}(2:end, 1:M/2) = b{m-1};
        b{m}(2:end, M/2+1:end) = fliplr(b{m-1});

        bint = bit2int(b{m}, m);  % columns as integers, first row is MSB

        x = [-M+1:2:M-1];         % values of PAM symbols
        A(bint+1) = x/sqrt(2*Es);   % ordered and scaled PAM symbols

        % determine subconstellations
        A0 = zeros(m, M/2);  A1 = zeros(m, M/2);
        for i = 1:m
            A0(i,:) = x(b{m}(i,:) == 0)/sqrt(2*Es);
            A1(i,:) = x(b{m}(i,:) == 1)/sqrt(2*Es);
        end
    end
end