%function [ x_Tikhonov,alpha,sigma2,mse_r,MSE_R,Alpha] = MMSE_Tikhonov_1( sigma2,A, y, P, R, x_prior)
function [ x_Tikhonov,QQ_Tikhonov,sigma2] = MMSE_Tik(A, y, x_prior, sigma2)
%=============================START PROGRAM===============================%
[m, n] = size(A);

    y1 =  y;
    A1 =  A ;
    %% SVD decomposition
    [U, S, V] = svd(A1);
    S = S(1:min(m,n), :);
    U = U(:, 1:min(m,n));
    
    %% Singular Value
    Lambda = diag(S);
    
    %% regularization parameter estimation based on MMSE
    alpha = 1;
    left_alpha = 0;
    right_alpha = 1e20;
    times = 1000;
    for iter = 1:times

        first_deri = sum(Lambda.^2.*(alpha*(V'*x_prior).^2.-sigma2)./(Lambda.^2.+alpha).^3);        
        
        if first_deri>=0
            right_alpha = alpha;
            alpha = (alpha+left_alpha)/2;            
        else
            left_alpha = alpha;
            alpha = (alpha+right_alpha)/2;
        end
    end
   %alpha=10.^(18);
   Alpha = alpha * ones(n,1);
    %% Tikhonov
    x_Tikhonov = V * diag(1./diag(S'*S+ diag(Alpha))) * S * U' * y1;

  Qalpha = V * diag(1./diag(S'*S+ diag(Alpha))) * V';

  v=y1-A1 *x_Tikhonov;
  sig1=v'*v-alpha.^2*x_Tikhonov'*(Qalpha-alpha.^2*Qalpha*Qalpha)*x_Tikhonov;
  sig2=m-n+alpha.^2*trace(Qalpha*Qalpha);
  sigma2=sqrt(sig1/sig2);

   N = A1'*A1;
   Dalpha = sigma2*Qalpha*N*Qalpha;
   Balpha = (Qalpha*A1'*A-eye(n))*x_Tikhonov;
   QQ_Tikhonov = Dalpha+Balpha*Balpha'; 

end
