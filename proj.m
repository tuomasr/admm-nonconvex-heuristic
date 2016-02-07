% See Takapoui et al. (2015)
% The projection stage
function x = proj(x0,no_cont,lb,ub)
    ints = x0((no_cont+1):end);
    temp = [x0(1:no_cont);round(ints)];
    x = min(max(temp,lb),ub);
end

