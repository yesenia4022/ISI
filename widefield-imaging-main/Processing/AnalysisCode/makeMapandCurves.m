function makeMapandCurves(locs)

%You must plot and click on the map before running this code

global ACQinfo Tens Tens_var Flim TCWin Fsymbol oppCollapse

Nlocs = length(locs(:,1));

for i = 1:Nlocs
    hold on
    plot(locs(i,2),locs(i,1),'ok','LineWidth',1.5,'MarkerSize',18)
    
    text(locs(i,2)+8,locs(i,1)-8,num2str(i),'FontWeight','bold')
    
end
hold off

figure
for i = 1:Nlocs
    
    W = TCWin;

    varflag = 0;
    if ~isempty(Tens_var)
        varflag = 1;
    end

    tdom = 0:length(Tens{1}(1,1,:))-1;
    tdom = tdom*ACQinfo.msPerLine/1000*ACQinfo.linesPerFrame;
    if isfield(ACQinfo,'stimPredelay')
        predelay = ACQinfo.stimPredelay;
        trialtime = ACQinfo.stimTrialtime;
        tdom = tdom-predelay;
    end

    nr = pepgetnorepeats;

    SEn = sqrt(length(Flim(1):Flim(2))*nr);  %standard error normalizer for tuning curve
    %
    %pos = round(get(event_obj,'Position')); %pos(1) is column dimension
    pos = [locs(i,2) locs(i,1)];

    %%%
    [tc tcourseHi tcourseLo axisdom blank legStr] = getpixeldata(pos,W);  %This does the work
    %%%

    tau = pos(2)*ACQinfo.msPerLine/1000;
    tdom = tdom + tau;

    subplot(2,Nlocs,i)
    if ~isempty(blank)
        plot([axisdom(1) axisdom(end)],[blank blank],'k'), hold on
    else
        %Even if no blank was shown we put a line at zero.
        plot([axisdom(1) axisdom(end)],[0 0],'k'), hold on
    end

    if ~varflag
        plot(axisdom,tc,'-o'), hold off
        if i == Nlocs
            legend(legStr)
        end
    else
        errorbar(axisdom,tc(id),sqrt(tc_var(id))/SEn,'b'), hold off
    end
    xlabel(Fsymbol)
    xlim([axisdom(1)-5 axisdom(end)+5])
    
    if i == 1
        ylabel('dF/F')
    end

    %Get 'orientation selectivity index' and put into the title
    if ~isempty(blank)
        tc = tc-blank;
    end


    if oppCollapse == 3
        [y x] = find(tc == max(tc(:)));
        tcdum = tc(:,x);
    else
        tcdum = tc;
    end

    TCSel = abs(sum(tcdum'.*exp(1i*2*axisdom*pi/180)));
    TCSel = TCSel/sum(tcdum);
    TCSel =  round(TCSel*100)/100;
    title(['Cell ' num2str(i) '    OSI = ' num2str(TCSel)])

    Fi = 3;

    subplot(2,Nlocs,i+Nlocs)
    if varflag
        dum_var = squeeze(sum(sum(Tens_var{idma}(yran,xran,:),1),2))/nopix/nr;
        errorbar(tdom(1:end-Fi),tcourseHi(1:end-Fi),sqrt(dum_var(1:end-Fi))), hold on
    else
        plot(tdom(1:end-Fi),tcourseHi(1:end-Fi),'.-'), hold on
    end

    if varflag
        dum_var = squeeze(sum(sum(Tens_var{idmi}(yran,xran,:),1),2))/nopix/nr;
        errorbar(tdom(1:end-Fi),tcourseLo(1:end-Fi),sqrt(dum_var(1:end-Fi)),'r')
    else
        plot(tdom(1:end-Fi),tcourseLo(1:end-Fi),'.-r')
    end

    if isfield(ACQinfo,'stimPredelay')
        ylimits = get(gca,'Ylim');
        plot([0 trialtime],[ylimits(1) ylimits(1)]+(ylimits(2)-ylimits(1))/10,'k')
    end
    hold off
    xlabel('sec')
    xlim([tdom(1) tdom(end-Fi)])
    
    if i == 1
        ylabel('dF/F')
    end

end
       