function [K,err,i] = mylqr(A,B,Q,R,varargin)
%% 
ni = nargin;
if(ni>=5)
    epslion=varargin{1};
    if (isempty(epslion)||epslion<=0)
        epslion=eps;
    end
    if(ni==6)
        iterations=varargin{2};
        if (isempty(iterations)||iterations<=0)
            iterations=inf;
        end
    else
        iterations=inf;
    end
else
    epslion=eps;
    iterations=inf;
end

%%求解黎卡提
[P,err,i]=RiccatiSolver(A,B,Q,R,epslion,iterations);

K=R\B'*P;
end

function [P,min_err,i]=RiccatiSolver(A,B,Q,R,epslion,iterations)

I=eye(size(A));
% iA=inv(I-A);
iA=(I-A)^(-1);
E=iA*(I+A);
G=2*iA^2*B;
H=R+B'*iA'*Q*iA*B;
W=Q*iA*B;

%迭代求解P
P0=zeros(size(A));
i=0;
err=10;
min_err=10;
min_P=zeros(size(A));
while(err>epslion&&i<iterations)
    i=i+1;
%     P=E'*P0*E-(E'*P0*G+W)*inv(G'*P0*G+H)*(E'*P0*G+W)'+Q;
    EP0=E'*P0;
    T=EP0*G+W;
    P=EP0*E-T/(G'*P0*G+H)*T'+Q; %%由于初值P0和求解精度问题，P可能发散最终为nan需不需要处理？如改变Q，R参数
    err=norm(P-P0);
    P0=P;
    if(min_err>=err)
        min_err=err;
        min_P=P;
    end
end
P=2*iA'*min_P*iA;
end




