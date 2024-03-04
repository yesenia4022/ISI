function intersect = getintersect(satel,dir,ro)

global S D

En = sqrt(sum(dir.*dir,2));
En = En*ones(1,length(dir(1,:)));
dir = dir./En;

%%%search%%%

S = satel;
D = dir;
intersect = intersectfitter(ro);
%%%%%%%%%%%







