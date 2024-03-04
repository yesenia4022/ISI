function [Tens oridom ang mag] = GetOriTensor(f0dum,hh)

%This is a modification of "Gprocessori"

global pepANA
%Each element of the cell array 'f0dum' is the average image for the
%corresponding condition


bflag = 0;
k = 1;
for(i=0:length(f0dum)-1)
    pepsetcondition(i)
    if(~pepblank)       %This loop filters out the blanks  
        for z = 1:length(pepANA.listOfResults{i+1}.values)  %loop through each loop parameter
            if strcmp(pepANA.listOfResults{i+1}.symbols(z),'ori')
                paramID = z;
            end
        end
        v = pepgetvalues;
        ori(k) = v(paramID);
        
        f0{k} = f0dum{i+1};
        k = k+1;
    else
        f0blank = f0dum{i+1};
        bflag = 1;
    end
end

%Make direction into orientation:

ori = angle(exp(1i*2*ori*pi/180))*180/pi;
ori = ori/2;
ori = ori + 90*(1-sign(ori+.01));

ori = round(ori);
oridom = unique(ori);

Tens = zeros(length(f0{1}(:,1)),length(f0{1}(1,:)),length(oridom));
for k = 1:length(oridom)
    id = find(ori == oridom(k));
    for i = 1:length(id)
        Tens(:,:,k) = Tens(:,:,k) + f0{id(i)}/length(id);
    end
end


mi = min(Tens,[],3);
for k = 1:length(Tens(1,1,:))
    if bflag == 1
        Tens(:,:,k) = (Tens(:,:,k)-f0blank)./f0blank;
    else
        Tens(:,:,k) = Tens(:,:,k)-mi;
    end
end

% su = sum(abs(Tens),3);
% for k = 1:length(Tens(1,1,:))
%     Tens(:,:,k) = Tens(:,:,k)./su;
% end

% id = find(Tens(:)<0);
% Tens(id) = 0;

% orimap = zeros(size(f0{1}));
% for k = 1:length(oridom)
%     orimap = orimap + Tens(:,:,k)*exp(1i*2*oridom(k)*pi/180);    %Linear combination
% end
% 
% ang = angle(orimap);
% ang = (ang+pi*(1-sign(ang)))/2*180/pi;
% 
% mag = abs(funcmap1);

%if a filter exists, use it...
if ~isempty(hh)
    id = find(isnan(orimap));
    orimap(id) = 0;
    orimap = ifft2(abs(fft2(hh)).*fft2(orimap));    
end

