function [Alpha, Mu, Sigma] = GMM_EM(Data, Alpha0, Mu0, Sigma0)
%% EM 迭代停止条件
loglik_threshold = 1e-10;
%% 初始化参数
[dim, N] = size(Data);
M = size(Mu0,2);
loglik_old = -realmax;
nbStep = 0;
 
Mu = Mu0;
Sigma = Sigma0;
Alpha = Alpha0;
Epsilon = 0.0001;
while (nbStep < 1200)
  nbStep = nbStep+1;
  %% E-步骤 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for i=1:M
    % PDF of each point
    Pxi(:,i) = GaussPDF(Data, Mu(:,i), Sigma(:,:,i));         
  end
 
  % 计算后验概率 beta(i|x)
  Pix_tmp = repmat(Alpha,[N 1]).*Pxi;
  Pix = Pix_tmp ./ (repmat(sum(Pix_tmp,2),[1 M])+realmin);
  Beta = sum(Pix);
  %% M-步骤 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for i=1:M
    % 更新权值
    Alpha(i) = Beta(i) / N;
    % 更新均值
    Mu(:,i) = Data*Pix(:,i) / Beta(i);
    % 更新方差
    Data_tmp1 = Data - repmat(Mu(:,i),1,N);
    Sigma(:,:,i) = (repmat(Pix(:,i)',dim, 1) .* Data_tmp1*Data_tmp1') / Beta(i);
    %% Add a tiny variance to avoid numerical instability
    Sigma(:,:,i) = Sigma(:,:,i) + 1E-5.*diag(ones(dim,1));
  end
 
%  %% Stopping criterion 1 %%%%%%%%%%%%%%%%%%%%
%  for i=1:M
    %Compute the new probability p(x|i)
%    Pxi(:,i) = GaussPDF(Data, Mu(:,i), Sigma(i));
%  end
  %Compute the log likelihood
%  F = Pxi*Alpha';
%  F(find(F<realmin)) = realmin;
%  loglik = mean(log(F));
  %Stop the process depending on the increase of the log likelihood
%  if abs((loglik/loglik_old)-1) < loglik_threshold
%    break;
%  end
%  loglik_old = loglik;
 
  %% Stopping criterion 2 %%%%%%%%%%%%%%%%%%%%
  v = [sum(abs(Mu - Mu0)), abs(Alpha - Alpha0)];
  s = abs(Sigma-Sigma0);
  v2 = 0;
  for i=1:M
    v2 = v2 + det(s(:,:,i));
  end
 
  if ((sum(v) + v2) < Epsilon)
    break;
  end
  Mu0 = Mu;
  Sigma0 = Sigma;
  Alpha0 = Alpha;
end
nbStep
