function Iq = mutual_information(hist)
% compute sum of mutual informations I(cm; z)
%
% the histogram contains in the first dimension (row) the received symbol
% y, in the second dimension (column) the transmitted bit x and in the
% third dimension the subchannel m = 1,...,Rm
% the mutual information is computed per subchannel 
% I = E[ ld( p(x,y)/(p(x) p(y)) )] and then added

[~,~,ldM] = size(hist);
Iq = 0;

for m = 1:ldM
    pxy = hist(:,:,m);
    pxy = pxy/sum(pxy(:));              % joint pmf p(x,y)

    px = sum(pxy);  py = sum(pxy,2);    % marginal probabilities p(x), p(y)
    pxpy = py*px;                       % p(x) p(y)

    pxy = pxy(:); pxpy = pxpy(:);       % column vectors
    ii = (pxy>0);                       % indices of non-zero probabilities

    pxy = pxy(ii);  pxpy = pxpy(ii);    % only positive probabilities
    
    Iq = Iq + sum(pxy .* log2(pxy./pxpy));
end
