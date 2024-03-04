function err = RadonTrxfitter_handle(param)

global RF odomglob pdomglob Gsfit;

ffitDori = Gsfit'*ones(1,length(pdomglob));

[xpg thetag] = meshgrid(pdomglob,odomglob);

prefpos = param(1)*cos(thetag*pi/180-param(2));
ffitDpos = exp( -(xpg-prefpos).^2/(2*param(3).^2) );

ffit = param(4)*ffitDpos.*ffitDori + param(5);

err = sum((ffit(:)-RF(:)).^2);