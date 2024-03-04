function [mapmu mapsig vertDomain] = retVScolorDist(ret,cmap)

[mat xdom ydom] = smoothscatter(ret,cmap,1,.05,[-45 45],[-1.5 1.5]);

mat = mat-min(mat(:));
mat = mat/max(mat(:));
imagesc(xdom,ydom,1-mat), colormap gray
xlim([-45 45]), ylim([-1.1 1.1])
axis xy

%vertBinEdges = [-45:15:45];

vertBinEdges = [ ];
prcDom = 0:10:100;
for i = 1:length(prcDom)
    vertBinEdges = [vertBinEdges prctile(ret,prcDom(i))];
end

clear mapmu mapsig vertDomain
for j = 1:length(vertBinEdges)-1
    
    id = find(ret>=vertBinEdges(j) & ret<vertBinEdges(j+1));
    cmapdum = cmap(id);
    vertdum = ret(id);
    
    mi = prctile(cmapdum,5);
    ma = prctile(cmapdum,95);
    idOutlier = find(cmapdum<mi | cmapdum>ma);
    vertdum(idOutlier) = [];
    cmapdum(idOutlier) = [];
    
    vertDomain(j) = nanmedian(vertdum);
    mapmu(j) = nanmedian(cmapdum);
    mapsig(j) = nanstd(cmapdum);
end

hold on,
   
errorbar(vertDomain,mapmu,mapsig,'r')

