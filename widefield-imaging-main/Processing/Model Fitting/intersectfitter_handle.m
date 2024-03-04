function errout = intersectfitter_handle(ro)

global S D err iter;
 
iter = iter + 1;

ro = ones(length(S(:,1)),1)*ro;

Sn = S-ro;
En = sqrt(sum(Sn.*Sn,2));
En = En*ones(1,length(Sn(1,:)));

Sn = Sn./En;


errdum = sum(Sn(:).*D(:),2);

err(iter) = -sum(abs(errdum));

errout = err(iter);