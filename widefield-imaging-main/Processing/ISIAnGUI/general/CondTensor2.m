function [y y_var] = CondTensor2(b,shiftflag,normflag,varflag)

%Compute the tensor for each condition

%The new version builds the entire Tensor, as opposed to truncating within a time interval 
%
%b is a 2D vector corresponding the the beginning and end of
%the baseline subtraction images, in milliseconds. e.g. varargin = {[0 500]} sums
%the images from 0 to .5 seconds for each repetition and then subtracts it
%from the mean response in the repeat.
%
%shiftflag performs movement correction 
%Rflag fits a line to the red/green scatter plot and subtracts this trend
%from the data in the green channel

global ACQinfo bsflag repDom

nc = getnoconditions;

alignCh = 2;
if shiftflag
    chvec = [0 0 0 0];
    chvec(alignCh) = 1;
    temp = getExptMean(chvec,1); %Use first channel for alignment (use first trial as the template)
    temp = temp{1};
else
    temp = [];
end

if normflag
    [RGslope RGbase mbest nbest] = RtoG_trx2(shiftflag,temp,alignCh);
    sigIm = 0;
    nIm = 0;
end

framePer = ACQinfo.linesPerFrame*ACQinfo.msPerLine;  %frame period in ms
b = b+getparam('predelay')*1000;  %user input is relative to stimulus onset, not trial beginning

bframe1 = floor(b(1)/framePer) + 1;
bframe2 = ceil(b(2)/framePer) + 1;


y = cell(1,nc);

for c = 1:nc

    y{c} = 0;

    if stimblank(c);      %Identify blank condition since it will have a unique number of repeats
        repeatDomain = 1:getnorepeats(nc);        
        %The following is to take only the blank repeats that are contained
        %within the selected repeats 
        if repDom(end) < getnorepeats(1)
            bperr = getnorepeats(nc)/(getnorepeats(1));    %blanks per repeat
            repeatDomain = 1:floor(repDom(end)*bperr);
        end
    else
        repeatDomain = repDom;
    end

    for rid = 1:length(repeatDomain)
        r = repeatDomain(rid);
        
        if ~normflag
            CHs = GetTrialData([1 0 0 0],[c r]);  %don't use syncs here
        else
            %requires second channel
            CHs = GetTrialData([1 1 0 0],[c r]);  %don't use syncs here
        end

        %Apply movement correction
        if shiftflag
            imdum = mean(CHs{1}(:,:,2:end-2),3);
            [mbest nbest] = getShiftVals(imdum,temp);  %get the transformation for this trial
            for z = 1:length(CHs{1}(1,1,:))
                if ~normflag  %shift values were already computed above if normflag is set
                    CHs{1}(:,:,z) = circshift(CHs{1}(:,:,z),[-mbest -nbest]); %transform
                else
                    CHs{1}(:,:,z) = circshift(CHs{1}(:,:,z),[-mbest(c,r,z) -nbest(c,r,z)]); %transform
                end
            end
        end
    
        if normflag
            for z = 1:length(CHs{1}(1,1,:))
%                 dumG = CHs{1}(:,:,z);
%                 dumG = (dumG-mean(dumG(id)))/std(dumG(id));
%                 CHs{1}(:,:,z) = dumG;

                %                     gain = CHs{1}(:,:,z)./CHs{2}(:,:,z);
                %                     gain = trimmean(gain(:),40);
                %                     normer = CHs{2}(:,:,z)*gain;
                %                     CHs{1}(:,:,z) = (CHs{1}(:,:,z) - normer);

                %CHs{1}(:,:,z) = (dumG - normer)/normer;
                %CHs{1}(:,:,z) = (dumG - min(dumG(:)))/min(dumG(:));

                CHs{1}(:,:,z) = CHs{1}(:,:,z) - (CHs{2}(:,:,z).*RGslope + RGbase);

                sigIm = sigIm + CHs{1}(:,:,z).^2;
                nIm = nIm + 1;
                
            end
        end


        if bsflag == 1

            bimg1 = mean(CHs{1}(:,:,bframe1:bframe2),3);

            for z = 1:length(CHs{1}(1,1,:))
                CHs{1}(:,:,z) = CHs{1}(:,:,z) - bimg1;   %% baseline subtraction
            end

            for z = 1:length(CHs{1}(1,1,:))
                CHs{1}(:,:,z) = CHs{1}(:,:,z)./(bimg1+eps);   %% baseline division (add 'eps', becauase...
                %there were wierd cases with pixel values of zero...
                %Scanimage bug)
            end

        end
        
        %Average repeats:  /nr is important for when blanks have different
        %number of reps
        y{c} = y{c} + CHs{1}/length(repeatDomain); 

        clear CHs

    end

end


if normflag    
    sigIm = sqrt(sigIm/(nIm-1));
    for i = 1:length(y)
        if ~isempty(y{i})
            for j = 1:length(y{i}(1,1,:))
                y{i}(:,:,j) = y{i}(:,:,j)./sigIm;
            end
        end
    end
end

if varflag
    y_var = CondTensor_sigma(b,shiftflag,normflag,y);
else
    y_var = [];
end
