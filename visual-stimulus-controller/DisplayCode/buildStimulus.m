function buildStimulus(cond,trial)

global DcomState

%Sends loop information and buffers

global looperInfo Mstate Pstate

mod = getmoduleID;

msg = ['B;' mod ';' num2str(trial)];

bflag = strcmp(looperInfo.conds{cond}.symbol{1},'blank');

%This is done just in case there are dependencies in the 'formula' on
%Mstate.
Mf = fields(Mstate);
for i = 1:length(fields(Mstate))
    eval([Mf{i} '= Mstate.'  Mf{i} ';' ])
end


if bflag  %if it is a blank condition
    
    msg = sprintf('%s;%s=%.4f',msg,'y_size',.5); %shrink it.
    msg = sprintf('%s;%s=%.4f',msg,'y_pos',-90); %push it off the screen
    %msg = sprintf('%s;%s=%.4f',msg,'contrast',0); %push it off the screen
    
    %Used to set constrast to 0, but that leaves a gray screen between
    %trials.  I decided that I wanted it set to the 'background' value for blanks. 
    
else %if it is not a blank condition
    
    %%%Send the size and position in Pstate in case last trial was a blank%%%
 
    pval = getParamVal('y_size',1);
    msg = sprintf('%s;%s=%.4f',msg,'y_size',pval); 
    pval = getParamVal('y_pos',1);
    msg = sprintf('%s;%s=%.4d',msg,'y_pos',pval); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Nparams = length(looperInfo.conds{cond}.symbol);
    for i = 1:Nparams
        pval = looperInfo.conds{cond}.val{i};
        psymbol = looperInfo.conds{cond}.symbol{i};
        msg = updateMsg(pval,psymbol,msg);
        eval([psymbol '=' num2str(pval) ';'])  %May be used to evaluate formula below (dependencies);
        
        eyefunc(psymbol,pval)  %This moves the eye shutters if its the right symbol
    end
    
    %Append the message with the 'formula' information
    fmla = looperInfo.formula
    id = find(fmla == ' ');
    if ~isempty(fmla)
        fmla = [';' fmla ';'];
        ide = find(fmla == '=');
        ids = find(fmla == ';');

        for e = 1:length(ide);

            delim1 = max(find(ids<ide(e)));
            delim1 = ids(delim1)+1;
            delim2 = min(find(ids>ide(e)));
            delim2 = ids(delim2)-1;

            Estate = 1;
            while Estate == 1
                try
                    eval([fmla(delim1:delim2) ';'])  %any dependencies should have been established above
                    Estate = 0;
                catch ME
                    
                    if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
                        varname = ME.message(1:end-2);
                        varname = varname(find(varname == '''')+1 : end);
                        pval = getParamVal(varname,1);  %get the value from Pstate
                        eval([varname '=' num2str(pval) ';'])
                        %eval([fmla(delim1:delim2) ';'])  %try again
                        
                        Estate = 1;
                    else
                        Estate = 0;
                    end
                end
            end
            
            psymbol_Fmla = fmla(delim1:ide(e)-1);
            pval_Fmla = eval(psymbol_Fmla);
            
            eyefunc(psymbol_Fmla, pval_Fmla)  %This moves the eye shutters if its the right symbol
            
            msg = updateMsg(pval_Fmla,psymbol_Fmla,msg);
        end
    end
    
end

msg = [msg ';~'];  %add the "Terminator"

fwrite(DcomState.serialPortHandle,msg);


function eyefunc(sym,bit)
id = find(sym == ' ');
sym(id) = [];

if strcmp(sym,'Leye_bit')
    moveShutter(1,bit);
    waitforDisplayResp;
elseif strcmp(sym,'Reye_bit')
    moveShutter(2,bit);
    waitforDisplayResp;
elseif strcmp(sym,'eye_bit')
    switch bit
        case 0
            moveShutter(1,1);
            waitforDisplayResp
            moveShutter(2,0);
            waitforDisplayResp
        case 1
            moveShutter(1,0);
            waitforDisplayResp
            moveShutter(2,1);
            waitforDisplayResp
        case 2
            moveShutter(1,1);
            waitforDisplayResp
            moveShutter(2,1);
            waitforDisplayResp
        otherwise
    end
end
    
    


function msg = updateMsg(pval,psymbol,msg)

global Pstate

id = find(psymbol == ' ');
psymbol(id) = []; %In case the user put in spaces with the entry

%Find parameter in Pstruct
idx = [];
for j = 1:length(Pstate.param)
    if strcmp(psymbol,Pstate.param{j}{1})
        idx = j;
        break;
    end
end

%change value based on looper
if ~isempty(idx)  %its possible that looper variable is not a grating parameter
    prec = Pstate.param{idx}{2};  %Get precision
    switch prec
        case 'float'
            msg = sprintf('%s;%s=%.4f',msg,psymbol,pval);
        case 'int'
            msg = sprintf('%s;%s=%d',msg,psymbol,round(double(pval)));
        case 'string'
            msg = sprintf('%s;%s=%s',msg,psymbol,pval);
    end
end
        
