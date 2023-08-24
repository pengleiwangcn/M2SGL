function y=main(Avd,c,NITER,nonval,lambda)
[n,~]=size(Avd{1});
[v, d] = size(Avd);
r=lambda;
nonval = min(d*v, nonval);

S=zeros(n,n);
for i=1:v
    S = S + Avd{i,1};
end
S=S/v;
D = diag(sum(S));
L = D - S;
F=eig1(L, c+1, 0);
F = F(:,2:c+1);
F = F./repmat(sqrt(sum(F.^2,2)),1,c);
for i=1:NITER
    %solveP
    P = zeros(v, d);
    tr_SA = cellfun(@(x) trace(S' * x), Avd, 'UniformOutput', true);
    tr_num_vec = reshape(tr_SA, 1, v*d);
    [~,ind] = sort(tr_num_vec,'descend');
    P(ind(1:nonval)) = 1;
    %solve alpha
    alpha=2*tr_SA;
    %solveSj
    U=zeros(n,n);
    distf = L2_distance_1(F',F');
    for v_=1:v
        for d_=1:d
            U=U+P(v_,d_)*alpha(v_,d_)*Avd{v_,d_};
        end
    end
    for g=1:n
        P_i=distf(g,:)-U(:,g)'*2/r;
        S(g,:)=EProjSimplex_new(-P_i*r/(lambda*4));
    end
    A = (S+S')/2;
    L = diag(sum(A)) - A;
    %solveF
    F_old=F;
    [F, ~, ev]=eig1(L, c, 0,1);
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 > 0.000001
        r = 2*r;
    elseif fn2 < 0.0000001
        r = r/2;
        F = F_old;
    else
        break;
    end
end
% sprintf('Break:%d,  NITET:%d',i,NITER)
y=conncomp(graph(sparse(A)));
% sprintf('Correct: %d,result: %d', c,length(unique(y)))
end
