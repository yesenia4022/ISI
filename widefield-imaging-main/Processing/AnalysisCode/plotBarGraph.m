function plotBarGraph(Y,logflag,varargin)

if logflag
    [mu sig] = getGeoStats(2.^(Y)); 
else
    mu = mean(Y); sig = std(Y);
end

if ~isempty(varargin)
    hdom = varargin{1};
    [h hdom] = hist(Y,hdom);
else
    [h hdom] = hist(Y);
end
    
bar(hdom,h)

mu = round(mu*100)/100; sig = round(sig*100)/100;

if logflag
    title(['geomean/sig = ' num2str(mu) '/' num2str(sig) '; N = ' num2str(length(Y))])
else
    title(['mean/sig = ' num2str(mu) '/' num2str(sig) '; N = ' num2str(length(Y))])
end
% a = 10.^str2num(get(gca,'XtickLabel'));
% a = round(a*100)/100;
%set(gca,'XtickLabel',num2str(a));  %This is buggy

set(gca,'Xtick',hdom(1:3:end))

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.5 .5 .5])
set(gca,'TickDir','out')
