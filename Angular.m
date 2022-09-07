function output = Angular(input)
    dy = input.dy;
    yMin = 0;
    yMax = 1;

    M = input.M;
    a = input.a; 
    mu = input.mu;

    n   = input.n;
    SCm = input.m;
    SCl = input.l;
    
    if isfield(input,'dLambda') == 0
        input.dLambda = 0;
        SCLambda = SCl * (SCl + 1);
    else
        SCLambda = input.Lambda + input.dLambda;
    end
    if isfield(input,'domegan') == 0
        input.domegan = 0;
        omegan = mu * (1 - M^2 * mu^2 / (2 * (n + SCl + 1)^2));
    else
        omegan = input.omegan + input.domegan;
    end
    
    y = (yMin : dy : yMax)';
    yLeft  = y + dy;
    yRight = y - dy;
    
    theta = pi * y;
    thetaLeft = pi * yLeft;
    thetaRight = pi * yRight;

    N = length(y);

    JElementsLeft    = pi * cot(thetaLeft);
    JElementsRight   = pi * cot(thetaRight);
    Wt0Elements      = pi^2 * ( - (SCm./sin(theta)).^2 + (omegan^2 - mu^2) * a^2 * cos(theta).^2);
    Wt1Elements      = pi^2 + 0 * theta;

    D2Left   = (1/(dy^2) - JElementsLeft /(2 * dy));
    D2Middle = (Wt0Elements + SCLambda * Wt1Elements - 2 /(dy^2));
    D2Right  = (1/(dy^2) + JElementsRight/(2 * dy));

    DLeft   = D2Left;
    DMiddle = D2Middle;
    DRight  = D2Right;
    LHS = spdiags([DLeft DMiddle DRight] ,-1:1,N,N);
    LHS(1,1) = 1;
    LHS(1,2) = -1 * (abs(SCm) == 0);
    LHS(end,end) = 1;
    LHS(end,end-1) = -1 * (abs(SCm) == 0);
    RHS = spdiags([Wt1Elements],[0],N,N);
    RHS(1,1) = 0;
    RHS(end,end) = 0;

    [output.Theta,output.dLambda] = eigs(LHS,-RHS,1,'smallestabs');
    output.Lambda = SCLambda;
    output.theta = theta;

end