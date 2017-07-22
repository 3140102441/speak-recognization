function [mm,mv,mc]=gaussmixm(m,v,w,z)
% GAUSSMIXM estimate mean and variance of the magnitude of a GMM
%
%  Inputs:  M(K,P)   is the mean vectors (one row per mixture)
%           V(K,P)   are diagonal covariances (one row per mixture) [ones(K,P)]
%        or V(P,P,K) are full covariance matrices (one per mixture)
%           W(K)     are the mixture weights [ones(K,1)/K]
%           Z(T,P)   each row is an origin to measure magnitude from [zeros(1,P)]
%                    T must equal 1 if P>2 and the MC output is requested.
%
% Outputs:  MM(T,1)  mean of sqrt((x-z(t))'*(x-z(t))) where x is the GMM random variate
%           MV(T,1)  variance of sqrt((x-z(t))'*(x-z(t)))
%           MC(T,T)  covariance matrix of sqrt((x-z(t))'*(x-z(t)))
%
% We approximate the magnitude of each mixture as a Nakagami-m distribution and
% match the moments of its squared variate.
% Answers are exact for P=1 or when all M=0. Accuracy improves otherwise for
% large P or M. Note that v_chimv() if V is a multiple of I for each mixture.

%      Copyright (C) Mike Brookes 2014
%      Version: $Id: gaussmixm.m 4992 2014-08-08 08:51:12Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent cart

[k,p]=size(m);                          % k = # mixtures, p = vector dimension
if nargin<4 || ~numel(z)
    z=zeros(1,p);
end
if nargin<3 || ~numel(w)
    w=ones(k,1);
end% default to uniform weights
if nargin<2 || ~numel(v)
    v=ones(k,p);                    % default to unit variance
end
[t,p1]=size(z);
if p~=p1 || (p>2 && t>1 && nargout>2)
    error('unable to calculate a covariance matrix');
end
w=w(:)/sum(w);                          % make the weights sum to 1
if p==1                                 % treat 1D case specially
    s=sqrt(v(:));                       % normalize mixture means by the std dev
    mt=m(:,ones(1,t))-z(:,ones(1,k))';
    mts=mt./s(:,ones(1,t));
    ncdf=normcdf(-mts);
    npdf=normpdf(-mts);
    mm=(mts.*(1-2*ncdf)+2*npdf)'*(s.*w); % exact mean of |X|
    mv=(mts.^2+1)'*(v(:).*w); % expected value of |x|^2 (not the variance yet)
    mc=diag(mv);
    mv=mv-mm.^2;
    if nargout>2 && t>1 % we need to calculate a covariance matrix
        for it=1:t
            for jt=1:it-1
                mc(it,jt)=w.'*((v(:)+mt(:,it).*mt(:,jt)).*(1-2*abs(ncdf(:,it)-ncdf(:,jt)))+ ...
                    2*s.*sign(mt(:,jt)-mt(:,it)).*(npdf(:,it).*mt(:,jt)-npdf(:,jt).*mt(:,it)));
                mc(jt,it)=mc(it,jt);
            end
        end
        mc=mc-mm*mm';                       % convert to variance
    end
else
    if ndims(v)>2 || size(v,1)>k            % full covariance matrix is supplied
        uv=zeros(k,p);
        u=zeros(p,p,k);
        for i=1:k                           % loop for each mixtrue component
            [u(:,:,i),d]=eig(v(:,:,i));   	% decorrelate the vector components
            uv(i,:)=diag(d)';               % variances of transformed version
        end
        ums=zeros(k,p);
        ms=zeros(k,t);
        vs=zeros(k,t);
        for it=1:t
            for i=1:k                       % loop for each mixtrue component
                ums(i,:)=((m(i,:)-z(it,:))*u(:,:,i)).^2;   	% squares of transformed means
            end
            ms(:,it)=sum(ums,2)+sum(uv,2);          % mean of |x|^2 for each mixture
            vs(:,it)=sum(uv.*(2*uv+4*ums),2);       % variance of |x|^2 for each mixture
        end
    else                                % else diagonal covariance matrix supplied
        ms=repmat(sum(m.^2,2)+sum(v,2),1,t)-2*m*z'+repmat(sum(z.^2,2)',k,1);  % mean of |x|^2 for each mixture x origin
        vs=sum(v.*(2*v+4*m.^2),2)-8*v.*m*z'+4*v*z.^2';      % variance of |x|^2 for each mixture x origin
    end
    nm=ms.^2./vs;                       % Nakagami-m parameter per mixture x origin
    mmk=(exp(gammaln(nm+0.5)-gammaln(nm)).*sqrt(ms./nm)); % mean of Nakagami-m distrbution per mixture x origin
    mm=mmk'*w;  % mean per mixture
    mv=ms'*w-mm.^2; % variance  of |x|  per mixture
    if nargout>2 && t>1 % we need to calculate a covariance matrix
        if isempty(cart)
            pth=mfilename('fullpath');
            load([pth '_cart']);
        end
        msd=sqrt(ms-mmk.^2); % std dev per mixture x origin
        mc=zeros(t,t,k);
        for i=1:k
            if ndims(v)>2 || size(v,1)>k            % full covariance matrix is supplied
                vi=v(:,:,i);
            else                                % else diagonal covariance matrix supplied
                vi=diag(v(i,:));
            end
            mc(:,:,i)=diag(ms(i,:));
            sc=1/sqrt(vi(1)+vi(4)); % scale factor to make trace(v)=1
            for it=1:t
                zit=z(it,:);
                for jt=1:it-1
                    z0=(zit+z(jt,:))*0.5; % new origin
                    zt=(zit-z0)*sc;
                    mt=(m(i,:)-z0)*sc;
                    [s,c]=atan2sc(zt(2),zt(1));
                    rm=[c -s; s c];
                    rm=rm*diag(2*(mt*rm>0)-1);
                    zt=zt*rm;
                    mt=mt*rm; % mt is always in the first quadrant
                    vt=rm'*vi*rm*sc^2;
                    d=2*abs(zt(1));
                    rho=eval(cart,[d mt/d vt(2)/sqrt(vt(1)*vt(4)) vt(1)]);
                    mc(it,jt,i)=rho*msd(i,it)*msd(i,jt)+ mmk(i,it)*mmk(i,jt);
                    mc(jt,it,i)=mc(it,jt,i);
                end
            end
        end
        mc=reshape(reshape(mc,t^2,k)*w,t,t)-mm*mm';    % calculate the variance as <X^2> - <X>^2
    end

end