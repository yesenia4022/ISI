function [y y_var] = CondTensor3(b,slowMo)

%Compute the tensor for each condition

%The new version builds the entire Tensor, as opposed to truncating within a time interval 
%
%b is a 2D vector corresponding the the beginning and end of
%the baseline subtraction images, in milliseconds. e.g. varargin = {[0 500]} sums
%the images from 0 to .5 seconds for each repetition and then subtracts it
%from the mean response in the repeat.
%
%slowMo performs movement correction 
%Rflag fits a line to the red/green scatter plot and subtracts this trend
%from the data in the green channel

global ACQinfo GUIhandles Analyzer bsflag repDom G_handles maskS cellS mbestall nbestall datadir

varflag = 0;


alignCh = 1;
if slowMo
    
    dum = GetTrialData([-inf inf],1);
    temp = median(dum,3);    
    
else
    temp = [];
end

B_Flim = getframeidx(b,1);

%Get acquired bit depth by pulling in a raw frame
ue = datadir(end-8:end-1);
fname = [datadir ue  '_' sprintf('%03d',0)];
fnamedum = [fname '_f1'];
load(fnamedum)
k = whos('im');
bitDepth = str2num(k.class(5:end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nt = getnotrials;
nc = getnoconditions;
y = cell(1,nc);
y2 = cell(1,nc);
ybase = cell(1,nc);

Px = 0;
Py = 0;

ExptMu = 0;

cellS = struct; %clear it
cellS.cellMat = cell(1,getnoconditions);

mbestall = [];
nbestall = [];

mbest = 0;
nbest = 0;
k = 0;
for t = 1:nt
    
    percentdone = round(t/nt*100);
    set(G_handles.status,'string',[num2str(percentdone) '%']), drawnow
    
    [c r] = getcondrep(t);

    if stimblank(c);      %Identify blank condition since it will have a unique number of repeats
        repeatDomain = 1:getnorepeats(nc);        
        %The following is to take only the blank repeats that are contained
        %within the selected repeats 
        if repDom(end) < getnorepeats(c)
            bperr = getnorepeats(nc)/(getnorepeats(1));    %blanks per repeat
            repeatDomain = 1:floor(repDom(end)*bperr);
        end
    else
        repeatDomain = repDom;
    end
    
    if ~isempty(find(r == repeatDomain))
        
        k = k+1;
        
        if repeatDomain(1) == r
            y{c} = 0;     
            y2{c} = 0; 
        end        
        
        CH_raw = GetTrialData([-inf inf],t);          
        if get(G_handles.negativeSignalFlag,'Value')            
            CH_raw = 2^bitDepth-CH_raw;            
        end

        CH = CH_raw;

        
        %ImLast = CH;
        
        ExptMu = ExptMu+CH;
        
        %Apply movement correction
        
        if slowMo
    
            imdum = mean(CH_raw{alignCh}(:,:,2:end-2),3); 

            
            [mbest nbest] = getShiftVals(imdum.^2,temp.^2,[mbest nbest])  %squaring seems to really help sometimes
            
            mbestall = [mbestall mbest];
            nbestall = [nbestall nbest];

            CH = circshift(CH,[-round(mbest) -round(nbest) 0]);

        end
        
        %I do this before the baseline normalization (below), because
        %averaging in the cell ROIs before the
        %divisive normalization can reduce noise. 'getCellStats' now
        %applies the baseline norm
        
        %Baseline normalization
        if bsflag == 1
                 
            bimg1 = mean(CH(:,:,B_Flim(1):B_Flim(2)),3);
            
            for z = 1:length(CH(1,1,:))
                CH(:,:,z) = CH(:,:,z) - bimg1;   %% baseline subtraction
            end
            
            id = find(bimg1(:) == 0);
            bimg1(id) = NaN;
            for z = 1:length(CH(1,1,:))
                CH(:,:,z) = CH(:,:,z)./bimg1;   %% baseline division
            end
            
        end
        
        
        %Average repeats:  /nr is important for when blanks have different
        %number of reps
        
%         for k = 1:length(CH(1,1,:))
%             hh = makeMapFilter;
%             CH(:,:,k) = ifft2(fft2(CH(:,:,k)).*abs(fft2(hh)));
%         end
        
        y{c} = y{c} + CH/length(repeatDomain);
        
        if varflag
            y2{c} = y2{c} + (CH.^2)/length(repeatDomain);  %to compute variance later
        end
               
  
        clear CHs
        
    end
    
end


y_var = cell(1,length(y));

if varflag
    for i = 1:length(y)
       y_var{i} = y2{i} - y{i}.^2;
    end
end


% bflag = stimblank(getnoconditions); %if a blank exists in this experiment
% if ~bsflag && bflag  %if baseline subtraction not checked and blanks are in this experiment    
%     dum = mean(y{end},3);
%     f0blank = zeros(size(y{1}));
%     for z = 1:length(y{1}(1,1,:))
%         f0blank(:,:,z) = dum;
%     end
%     
%     for c = 1:length(y)
%         y{c} = (y{c} - f0blank)./f0blank;
%     end
%     
% end



