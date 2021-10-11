function [StgD] = Func_pop_rest(rest_pop, StgD, iDist, yrDist, rDist)

% Assigning the disturbancce to the correct populations & times
            % dictate which pop is restored
            switch rest_pop
                case 'rest_res'
                    % first of the two reserve pops (pop1) restored
                    StgD(yrDist:yrDist+rDist-1,1) = iDist;
                    StgD(yrDist:end,[2,3,4]) = iDist;
                case 'rest_fish'
                    % first of the two fished pops (pop3) restored
                    StgD(yrDist:yrDist+rDist-1,3) = iDist;
                    StgD(yrDist:end,[1,2,4]) = iDist;
                case 'rest_both'
                    % first of the two fished pops (pop3) restored
                    StgD(yrDist:yrDist+rDist-1,[1,3]) = iDist;
                    StgD(yrDist:end,[2,4]) = iDist;
            end
end

