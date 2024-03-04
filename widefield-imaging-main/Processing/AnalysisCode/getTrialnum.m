function trial = getTrialnum(cond,rep)

%The input condition value goes from 1 to n (not 0 to n-1)

global pepANA

trial = pepANA.listOfResults{cond}.repeat{rep}.timeTag + 1;