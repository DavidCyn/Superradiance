function output = Radial(input)
    dx = input.dx;
    xMin = input.xMin;
    xMax = input.xMax;

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

    rp = M + sqrt(M^2 - a^2);
    rm = M - sqrt(M^2 - a^2);

    omegac = a * SCm / (2 * M * rp);

    x = (xMin : dx : xMax)';
    xLeft  = x + dx;
    xRight = x - dx;

    r = rp * (exp(x) + 1);
    rLeft  = rp * (exp(xLeft ) + 1);
    rRight = rp * (exp(xRight) + 1);

    N = length(x);

    HElements        =  ((r      - rp)./(r      - rm));
    HElementsLeft    =  ((rLeft  - rp)./(rLeft  - rm));
    HElementsRight   =  ((rRight - rp)./(rRight - rm));
    Vt0Elements      =  (SCm^2 * a^2)./((r - rm).^2) - HElements .*(mu^2 * r.^2 + SCLambda);
    Vt1Elements      = -(4 * M * a * SCm * r)./((r - rm).^2);
    Vt2Elements      =  ((r.^2 +  a^2).^2)./((r - rm).^2) - HElements * a^2;
    B0               = 1 + dx * omegac * 2 * 1i * (rp/(rp - rm));
    B1               =            - dx * 2 * 1i * (rp/(rp - rm));

    D2Left   = (1/(dx^2) - HElementsLeft /(2 * dx));
    D2Middle = (Vt0Elements - 2 /(dx^2));
    D2Right  = (1/(dx^2) + HElementsRight/(2 * dx));

    DLeft   = D2Left;
    DMiddle = D2Middle;
    DRight  = D2Right;
    D2 = spdiags([DLeft DMiddle DRight] ,-1:1,N,N);
    D2(1,1) = 0;
    D2(1,2) = 0;
    Vt1 = spdiags(Vt1Elements,0,N,N);
    Vt1(1,1) = 0;
    Vt2 = spdiags(Vt2Elements,0,N,N);
    Vt2(1,1) = 0;
    Id = speye(N,N);
    Z = sparse(N,N);
    S11 = Z;
    S12 = Z;
    S11(1,1) = 1;
    S12(1,2) = 1;


    A11 = S11 * B0 - S12 + D2;
    A12 = S11 * B1 + Vt1 + omegan * Vt2;
    A21 = Id * omegan;
    A22 = Id * (-1);

    A = [A11 A12;A21 A22];

    B11 = Z;
    B12 = Vt2 * (-1);
    B21 = Id  * (-1);
    B22 = Z;

    B = [B11 B12;B21 B22];

    [output.R,output.domegan] = eigs(A,B,1,'smallestabs');
    output.omegan = omegan;
    output.r = r;
    output.x = x;

end