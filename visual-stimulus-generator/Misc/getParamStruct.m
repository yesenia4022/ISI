function [PStruct] = getParamStruct

global Pstate

for i = 1:length(Pstate.param)
    
    feeld = Pstate.param{i}{1};
    
    eval(['PStruct.' feeld '= Pstate.param{i}{3} ;'])
    
%     if strcmp(Pstate.type,'PG') %only do this for periodic grater, just in case
%         %create a second structure that replaces the default grating variables with the
%         %ones in the plaid parameters
%         eval(['PStruct2.' feeld '= Pstate.param{i}{3} ;'])
%         if strcmp(feeld(end),'2') %if a plaid parameter
%             eval(['PStruct2.' feeld(1:end-1) '= Pstate.param{i}{3} ;']) %replace the grating parameter
%         end
%     end
        
end

