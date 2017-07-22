function [YM,Y]=lmultigauss(x,mus,sigm,c)
% [lYM,lY]=lmultigauss(X,mu,sigm,c)
%
% computes multigaussian log-likelihood
%
% X   : (LxT) data (columnwise vectors)
% sigm: (LxM) variances vector  (diagonal of the covariance matrix)
% mu  : (LxM) means
% c   : (Mx1) the weights

[L,T]=size(x);
M=size(c,1);

% repeating, changing dimensions:
X=permute(repmat(x',[1,1,M]),[1,3,2]);      % (T,L) -> (T,M,L) one per mixture

Sigm=permute(repmat(sigm,[1,1,T]),[3,2,1]); % (L,M) -> (T,M,L)

Mu=permute(repmat(mus,[1,1,T]),[3,2,1]);     % (L,M) -> (T,M,L)

%Y=squeeze(exp( 0.5.*dot(X-Mu,(X-Mu)./Sigm))); % L dissapears: (L,T,M) -> (T,M)
lY=-0.5.*dot(X-Mu,(X-Mu)./Sigm,3);

% c,const -> (T,M) and then multiply by old Y
lcoi=log(2.*pi).*(L./2)+0.5.*sum(log(sigm),1); % c,const -> (T,M)
lcoef=repmat(log(c')-lcoi,[T,1]);

YM=lcoef+lY;            % ( T,M ) one mixture per column
Y=lsum(YM,2);           % add mixtures

