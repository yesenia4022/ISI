function [axismap] = GprocessAxis(f0,hh)

global bsflag Analyzer symbolInfo G_handles
%Each element of the cell array 'f0dum' is the average image for the
%corresponding condition

idsym = symbolInfo.ID(1);
idsym2 = symbolInfo.ID(2);

dom2index = get(G_handles.secCollapse,'value') - 1;  %Use zero to indicate the mean

bflag = stimblank(getnoconditions); %if a blank exists in this experiment
if bflag
    f0blank = f0{end};
    f0(end) = [];
end

for i = 1:length(f0)
    axisdomcond(i) = Analyzer.loops.conds{i}.val{idsym};
    dom2cond(i) = Analyzer.loops.conds{i}.val{idsym2};
end

idslice = 1:length(f0);
if dom2index>0
    dom2 = eval(Analyzer.L.param{idsym2}{2});
    idslice = find(dom2(dom2index) == dom2cond);
end

dim = size(f0{1});
Tens = zeros(dim(1),dim(2),length(idslice),'single'); %preallocate
for k = 1:length(idslice)
    Tens(:,:,k) = f0{idslice(k)};
end

%axisdomcond = axisdomcond*180/max(axisdomcond);

%I used to normalize by the blank if baseline subtraction was not checked.
%This is now done on the Tensor beforehand in CondTensor3

mi = min(Tens,[],3);
for k = 1:length(Tens(1,1,:))
    if bflag == 1  %if baseline subtraction button is set
        %Tens(:,:,k) = Tens(:,:,k)-f0blank;
    else
        Tens(:,:,k) = Tens(:,:,k)-mi;
    end
end
% % 

if ~isempty(hh)
    for k = 1:length(Tens(1,1,:))        
        Tens(:,:,k) = ifft2(abs(fft2(hh)).*fft2(Tens(:,:,k)));
    end
end

id = find(Tens(:)<0);
Tens(id) = 0;

% su = (sum(Tens,3));
% for k = 1:length(Tens(1,1,:))
%     Tens(:,:,k) = Tens(:,:,k)./su;
% end



axismap = zeros(size(f0{1}));
for k = 1:length(idslice)
    %axismap = axismap + f0{k}*exp(1i*2*axisdomcond(k)*pi/180);    %Linear combination
    axismap = axismap + Tens(:,:,k)*exp(1i*2*axisdomcond(idslice(k))*pi/180);    %Linear combination
end



%axismap = mean(Tens.^2,3).*exp(1i*angle(axismap));

%if a filter exists, use it...
% id = find(isnan(axismap));
% axismap(id) = 0;
% if ~isempty(hh)
%     axismap = ifft2(abs(fft2(hh)).*fft2(axismap));
% end

% CBnorm = 1;
% if CBnorm
%     su = sum((Tens),3);
%     if ~isempty(hh)
%         su = ifft2(abs(fft2(hh)).*fft2(su));
%     end
%     id = find(su(:)<0);
%     su(id) = NaN;
%     axismap = axismap./su;
% end


