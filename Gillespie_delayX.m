function [Xt, tspan, Xbirth, Xdeath] = Gillespie_delayX(theta, T_max, Iter_max)
X = 0;
XList = nan(Iter_max,1);
Xbirth = zeros(T_max,1);
Xdeath = zeros(T_max,1);
TList = nan(Iter_max,1);
currentTime = 0;
stackTimeX = [];
for ii = 1:Iter_max
    lambda = theta(1) + theta(2) * X;
    if lambda == 0
        break;
    else
        r1 = rand();
        Tbin = - 1/lambda * log(r1);
        currentTime = currentTime + Tbin;
        stackTimeX = sort(stackTimeX);
        if currentTime>T_max
            break;
        else
            if ~isempty(stackTimeX)
                minStack = stackTimeX(1);
            else
                minStack = Inf;
            end
            if currentTime < minStack
                r2 = rand();
                if lambda*r2 < theta(1)
                % X birth
                    XList(ii) = X;
                    TList(ii) = currentTime;
                    stackTimeX = [stackTimeX, currentTime + gamrnd(theta(3), theta(4))];
                    %stackTimeX = [stackTimeX, currentTime];
                else
                    % X death
                    X = X-1;
                    XList(ii) = X;
                    TList(ii) = currentTime;
                    Xdeath(ceil(currentTime)) = Xdeath(ceil(currentTime)) + 1;
                end
            else
                X = X+1;
                XList(ii) = X;
                TList(ii) = minStack;
                currentTime = minStack;
                if length(stackTimeX) == 1
                    stackTimeX = [];
                else
                    stackTimeX = stackTimeX(2:end);
                end
              Xbirth(ceil(currentTime)) = Xbirth(ceil(currentTime)) + 1;
            end
        end
    end
end
tspan = TList(~isnan(TList));
Xt = XList(~isnan(XList));
end