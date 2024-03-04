function [Tens pardom f0blank] = GetF0Tensor(f0dum)

%This is a modification of "Gprocessori"

global pepANA popState
%Each element of the cell array 'f0dum' is the average image for the
%corresponding condition


bflag = 0;
k = 1;
for(i=0:length(f0dum)-1)
    pepsetcondition(i)
    if(~pepblank)       %This loop filters out the blanks  
        for z = 1:length(pepANA.listOfResults{i+1}.values)  %loop through each loop parameter
            if strcmp(pepANA.listOfResults{i+1}.symbols(z),popState.funcSymbol)
                paramID = z;
            end
        end
        v = pepgetvalues;
        par(k) = v(paramID);
        
        f0{k} = f0dum{i+1};
        k = k+1;
    else
        f0blank = f0dum{i+1};
        bflag = 1;
    end
end

%If functionality is orientation, then combine across opposite directions:
if strcmp(popState.funcSymbol,'ori')
    par = angle(exp(1i*2*par*pi/180))*180/pi;
    par = par/2;
    par = par + 90*(1-sign(par+.01));
end

par = round(par*100)/100;
pardom = unique(par);

Tens = zeros(length(f0{1}(:,1)),length(f0{1}(1,:)),length(pardom));
for k = 1:length(pardom)
    id = find(par == pardom(k));
    for i = 1:length(id)
        Tens(:,:,k) = Tens(:,:,k) + f0{id(i)}/length(id);
    end
end

%Normalize by the blank response
mi = min(Tens,[],3);
for k = 1:length(Tens(1,1,:))
    if bflag == 1
        Tens(:,:,k) = (Tens(:,:,k)-f0blank)./f0blank;
    else
        Tens(:,:,k) = Tens(:,:,k)-mi;
    end
end



