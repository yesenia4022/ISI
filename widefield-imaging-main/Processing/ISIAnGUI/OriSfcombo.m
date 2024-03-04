%function OriSfcombo

global bw f0m f0m_var funcmap bcond ACQinfo maskS TCWin symbolInfo Analyzer G_handles


set(G_handles.HPflag,'Value',0);
set(G_handles.LPflag,'Value',1);
set(G_handles.Lwidth,'string','5');
hh = makeMapFilter;

anatomyflag = 0;
bwCellPlot = ones(size(funcmap));


%% get orimap
for i = 1:length(Analyzer.loops.conds{1}.symbol)
    if strcmp(Analyzer.loops.conds{1}.symbol{i},'ori');
        idsym = i;
        break
    end
end

symbolInfo.ID(1) = idsym;
set(G_handles.primSymbol,'value',1); 
set(G_handles.secSymbol,'value',2); 

setsymbolstruct

funcmap = GprocessAxis(f0m,hh);  %output is a vector image

oriang = angle(funcmap);
orimag = abs(funcmap);
oriang = (oriang+pi*(1-sign(oriang)))/2*180/pi;  %Put domain as [0 180].


%% get sfmap

for i = 1:length(Analyzer.loops.conds{1}.symbol)
    if strcmp(Analyzer.loops.conds{1}.symbol{i},'s_freq');
        idsym = i;
        break
    end
end

symbolInfo.ID(1) = idsym;
set(G_handles.primSymbol,'value',2); 
set(G_handles.secSymbol,'value',1); 

setsymbolstruct

funcmap = GprocessLog(f0m,bwCellPlot,hh);   %output is complex
sfmag = real(funcmap);
sfpref = imag(funcmap);


%% Get mask for data selection

dim = size(f0m{1});
Tens = zeros(dim(1),dim(2),length(f0m),'single'); %preallocate
Tens_var = Tens;
for k = 1:length(f0m)
    Tens(:,:,k) = ifft2(fft2(f0m{k}).*abs(fft2(hh)));
    Tens_var(:,:,k) = ifft2(fft2(f0m_var{k}).*abs(fft2(hh)));
    
%     Tens(:,:,k) = f0m{k};
%     Tens_var(:,:,k) = f0m_var{k};
end

[pk_dF idpk] = max(Tens(:,:,1:end-1),[],3);
pkSE_dF = zeros(size(idpk));
for i = 1:size(Tens_var,1)
    for j = 1:size(Tens_var,2)
        pkSE_dF(i,j) = sqrt(Tens_var(i,j,idpk(i,j)))/sqrt(getnorepeats(1));  %standard error of best response at each pixel
    end
end

base_dF = Tens(:,:,end);
baseSE_dF = sqrt(Tens_var(:,:,end))/sqrt(getnorepeats(getnoconditions));

dprime = (pk_dF-base_dF)./(pkSE_dF+baseSE_dF);

%% plot

logdom = logspace(log10(.5),log10(8),30);
oridom = [0:5:175];

dpthresh = 1;
dprimemask = zeros(size(dprime));
dprimemask(find(dprime>dpthresh)) = 1;

figure
Gplotaxismap(orimag.*dprimemask.*bw,oriang,anatomyflag), title(symbolInfo.str{1},'FontWeight','bold','FontSize',15);
hold on
sfdum = sfpref;
sfdum = sfdum.*dprimemask;
contour(sfdum,logdom,'k')
hold off
axis square
 
figure
Gplotlogmap(sfmag.*dprimemask.*bw,sfpref,anatomyflag)
hold on
oriangdum = oriang*64/180/2;  
oriangdum = oriangdum.*dprimemask;
contour(oriangdum,[oridom]*64/180/2,'k')
axis square


figure
contour(oriangdum,[oridom]*64/180/2,'k')
hold on
contour(sfdum,logdom,'r')
axis ij

%%
dim = size(oriang);
%dorix = oriang(:,3:end) - oriang(:,1:end-2);
%doriy = oriang(3:end,:) - oriang(1:end-2,:);

dorix = oridiff(oriang(:,3:end)*pi/180,oriang(:,1:end-2)*pi/180);
doriy = oridiff(oriang(3:end,:)*pi/180,oriang(1:end-2,:)*pi/180);

dorix = [zeros(dim(1),1) dorix zeros(dim(1),1)];
doriy = [zeros(1,dim(2)); doriy; zeros(1,dim(2))];

magorigrad = sqrt(dorix.^2 + doriy.^2);
angorigrad = atan(doriy./dorix);
id = find(angorigrad<0);
angorigrad(id) = angorigrad(id) + pi;
angorigrad = angle(exp(1i*angorigrad)); %make it an axis


dsfx = log2(sfpref(:,3:end)./sfpref(:,1:end-2));
dsfy = log2(sfpref(3:end,:)./sfpref(1:end-2,:));

dsfx = [zeros(dim(1),1) dsfx zeros(dim(1),1)];
dsfy = [zeros(1,dim(2)); dsfy; zeros(1,dim(2))];

magsfgrad = sqrt(dsfx.^2 + dsfy.^2);
angsfgrad = atan(dsfy./dsfx);
id = find(angsfgrad<0);
angsfgrad(id) = angsfgrad(id) + pi;
angsfgrad = angle(exp(1i*angsfgrad)); %make it an axis

%id = find(orimag>prctile(orimag(:),60));
%id = find(sfmag>prctile(sfmag(:),20));
id = find(dprime>dpthresh & bw==1);

dax = oridiff(angorigrad,angsfgrad);
figure,hist(abs(dax(id))*180/pi)

%figure,scatter(magsfgrad(id),magorigrad(id),'.')
[r p] = corrcoef(magsfgrad(id),magorigrad(id))
