function yhat = LinPred(y,sp,ord,del,L)

st = L+ord+del;
yhat = y;
for i = st:length(y)
    H = zeros(L,ord);
    for j = 1:ord
        id = i-del-L-j+1:i-del-j;
        pc = y(id);
        H(:,j) = pc(:);
    end
    id = i-L+1:i;
    yout = y(id);
    param = inv(H'*H)*H'*yout(:);
    yhat(i) = H(end,:)*param;
end
    