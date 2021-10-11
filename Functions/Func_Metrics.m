function [NDpc, RtimeA, RateA, RtimeB, RateB,...
     NDpcIn, RtimeAIn, RateAIn,...
     RtimeBIn, RateBIn] = Func_Metrics(N, tPost, yrDist, lDist, InBreak,...
                                       Interpt, recov_tol)
    
    sNPD = N; 
    
    % 1) total % knocked-down
    [TD,~] = min(sNPD(yrDist:end));
    NDpc = TD/sNPD(yrDist)*100;
    
    % 2A) time since disturbance end to yr closest pre-disturbance level
    % note mi = time since the disturbance ENDED until the point the
    % populaiton is considered recovered
    [~, mi] = min(abs(sNPD(yrDist+lDist+1:end) - sNPD(yrDist)));
    RtimeA = max([0,mi-1]);
    % 3A) recovery rate ([recovered abundance - KD]/time to recover)
    RateA = (sNPD(yrDist+lDist+mi-1)-TD)/RtimeA;
    
    % 2B) time since disturbance end to x% within pre-disturbance level
    [~, mi] = min(abs(sNPD(yrDist+lDist+1:end) - recov_tol.*sNPD(yrDist)));
    RtimeB = max([0,mi-1]);
    % 3B) recovery rate ([recovered abundance - KD]/time to recover)
    RateB = (sNPD(yrDist+lDist+mi-1)-TD)/RtimeB;
    
    % Intertpolate to try for better estimates
    % Interpolate abundances based on this
    sNPDinterp = interp1(1:tPost, sNPD, Interpt);
    InterpYrDist = find(Interpt == yrDist);
    
    % 1) total % knocked-down
    [TD,~] = min(sNPDinterp(InterpYrDist:end));
    NDpcIn = TD/sNPDinterp(InterpYrDist)*100;
            
    % 2A) time since disturbance end to yr closest pre-knockdown level
    [~, mi] = min(abs(sNPDinterp((InterpYrDist+lDist/InBreak+1):end) - sNPDinterp(InterpYrDist)));
    RtimeAIn = max([0,mi-1]).*InBreak;
    % 3A) recovery rate ([recovered abundance - KD]/time to recover)
    RateAIn = (sNPDinterp(InterpYrDist+(lDist/InBreak)+mi-1)-TD)/RtimeAIn;
    
	% 2B) time since disturbance end to x% within pre-disturbance level
    [~, mi] = min(abs(sNPDinterp((InterpYrDist+lDist/InBreak+1):end) - recov_tol.*sNPDinterp(InterpYrDist)));
    RtimeBIn = max([0,mi-1]).*InBreak;
    % 3B) recovery rate ([recovered abundance - KD]/time to recover)
    RateBIn = (sNPDinterp(InterpYrDist+(lDist/InBreak)+mi-1)-TD)/RtimeBIn;
    
    
end 