function xDot = networkedModel(~, x, N, v, k, r, a, u, C, D)

xDot = zeros(2 * N, 1);
    
    for i = 1:N
        x1iIndex = 2 * i - 1;
        x2iIndex = 2 * i;
        
        x1Sum = 0;
        x2Sum = 0;

        for j = 1:N
            if j ~= i
                x1jIndex = 2 * j - 1;
                x2jIndex = 2 * j;
                
                % Summing over differences between regions
                x1Sum = x1Sum + C(i,j) * (x(x1iIndex) - x(x1jIndex));
                x2Sum = x2Sum + D(i,j) * (x(x2iIndex) - x(x2jIndex));
            end
        end
        
        % Susceptible population equation (x1i)
        x1Dot = -1 * (v * x(x1iIndex) * x(x2iIndex))/(k(1) + x(x2iIndex)) ...
                + a * x(x2iIndex) + u(i, 1) - x1Sum;
        
        % Infected population equation (x2i)
        x2Dot = (v * x(x1iIndex) * x(x2iIndex))/(k(1) + x(x2iIndex)) ...
                - (r * x(x2iIndex)) / (x(x2iIndex) + k(2)) ...
                - a * x(x2iIndex) + u(i, 2) - x2Sum;
        
        % Store for output
        xDot(x1iIndex) = x1Dot;
        xDot(x2iIndex) = x2Dot;
    end
end