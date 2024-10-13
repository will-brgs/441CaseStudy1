function xDot = networkedModel(~, x, N, v, k, r, a, u)

xDot = zeros(2 * N, 1);
    
    for i = 1:N
        x1Index = 2 * i - 1;
        x2Index = 2 * i;
        
        % x1
        x1Dot = (v * x(x1Index) * x(x2Index)) / (k(1) + x(x2Index)) + a * x(x2Index) + u(i, 1);
        % x2
        x2Dot = (v * x(x1Index) * x(x2Index)) / (k(1) + x(x2Index))- (r * x(x2Index)) / (x(x2Index) + k(2)) ...
                - a * x(x2Index) + u(i, 2);
        
        % Store for output
        xDot(x1Index) = x1Dot;
        xDot(x2Index) = x2Dot;
    end
end