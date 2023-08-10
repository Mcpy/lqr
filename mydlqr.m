function [K,err,i] = mydlqr(A,B,Q,R,varargin)
%% A B Q R epslion iterations
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

K=(R+B'*P*B)\B'*P*A;
end

function [P,min_err,i]=RiccatiSolver(A,B,Q,R,epslion,iterations)

%迭代求解P
% P0=zeros(size(A));
P0=Q;
i=0;
err=10;
min_err=10;
min_P=zeros(size(A));
while(err>epslion&&i<iterations)
    i=i+1;
    P=Q+A'*P0*A-(A'*P0*B)/(R+B'*P0*B)*(B'*P0*A);
    err=norm(P-P0);
    P0=P;
    if(min_err>=err)
        min_err=err;
        min_P=P;
    end
end
P=min_P;
end




