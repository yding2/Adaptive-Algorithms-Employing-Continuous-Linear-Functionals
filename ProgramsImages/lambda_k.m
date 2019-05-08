function x=lambda_k(k,r)
    x=abs(k).^(-r);
    x(k==0)=1;