function hh = makeMapFilter

global G_handles maskS ACQinfo

togstateHP = get(G_handles.HPflag,'Value');
togstateLP = get(G_handles.LPflag,'Value');

sizeIM = [ACQinfo.linesPerFrame ACQinfo.pixelsPerLine];
if togstateHP == 1
    Hwidth = str2double(get(G_handles.Hwidth,'string'));
    ind = get(G_handles.HPWind,'value');

    switch ind
        case 1
            H = -fspecial('gaussian',sizeIM,Hwidth);
            H(round(sizeIM/2),round(sizeIM/2)) = 1+H(round(sizeIM/2),round(sizeIM/2));
        case 2
            H = zeros(sizeIM);
            Hd = hann(Hwidth)*hann(Hwidth)';
            Hd = -Hd./sum(Hd(:));
            Hd(round(Hwidth/2),round(Hwidth/2)) = 1+Hd(round(Hwidth/2),round(Hwidth/2));
            H(1:Hwidth,1:Hwidth) = Hd;
        case 3
            H = zeros(sizeIM);
            Hd = -fspecial('disk',round(Hwidth/2));
            Hsize = length(Hd(1,:));  %~=Hwidth
            Hd(round(Hsize/2),round(Hsize/2)) = 1+Hd(round(Hsize/2),round(Hsize/2));
            H(1:Hsize,1:Hsize) = Hd;
    end
    if togstateLP == 0
        hh = H;   %Eliminate phase information
    end
end

if togstateLP == 1
    Lwidth = str2double(get(G_handles.Lwidth,'string'));
    ind = get(G_handles.LPWind,'value');

    switch ind
        case 1
            L = fspecial('gaussian',sizeIM,Lwidth);
        case 2
            L = zeros(sizeIM);
            Ld = hann(Lwidth)*hann(Lwidth)';
            Ld = Ld./sum(Ld(:));
            L(1:Lwidth,1:Lwidth) = Ld;
        case 3
            L = zeros(sizeIM);
            Ld = fspecial('disk',round(Lwidth/2));
            Lsize = length(Ld(1,:));
            L(1:Lsize,1:Lsize) = Ld;
    end
    if togstateHP == 0
        hh = L;   %Eliminate phase information
    else
        hh = ifft2(fft2(L).*fft2(H));   %Take mag because phase gives a slight shift.
    end
end

if ~or(togstateLP,togstateHP)
    hh = [];
end
