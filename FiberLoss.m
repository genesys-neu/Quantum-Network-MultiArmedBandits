function [rho_fibre,p_fibre] = FiberLoss(p_init,f_attenuation,l_elem,rho_init_EP)
    %%%% EPinRightLeftMemory signifies whether EP is in left (L) or right (R) memory of the node
    %%%% Fibre Loss Noisy 
    p_fibre = ( 1-(1-p_init) * (10^(-f_attenuation * (l_elem) / 10)) );   %% Fibre Loss Probability due to quantum channel
    
    %%%% Operators
    I = [1,0;0,1]; %%%% I Operator
    X = [0,1;1,0]; %%%% X Operatr (Bit flip)
    Y = [0,-1i;1i,0]; %%%% Y Operator
    Z = [1,0;0,-1];   %%%% Z Operator (Phase flip)
    
    %%%% Fibre Loss Possibility on Different Operator
    p_fibreX1 = (1/3) * p_fibre;
    p_fibreY1 = (1/3) * p_fibre;
    p_fibreZ1 = (1/3) * p_fibre;
   
    p_fibreX2 = (1/3) * p_fibre;
    p_fibreY2 = (1/3) * p_fibre;
    p_fibreZ2 = (1/3) * p_fibre;

    %%%% Qubits Density Matrix Change through Depolarization noise
    rho_fibre = ((1 - p_fibreX2 - p_fibreY2 - p_fibreZ2) * kron(I,I) * rho_init_EP * kron(I,I)') +...
                (p_fibreX2 * kron(I,X) * rho_init_EP * kron(I,X)') +...
                (p_fibreY2 * kron(I,Y) * rho_init_EP * kron(I,Y)') +...
                (p_fibreZ2 * kron(I,Z) * rho_init_EP * kron(I,Z)');
                
    rho_fibre = ((1 - p_fibreX1 - p_fibreY1 - p_fibreZ1) * kron(I,I) * rho_fibre * kron(I,I)') +...
                (p_fibreX1 * kron(X,I) * rho_fibre * kron(X,I)') +...
                (p_fibreY1 * kron(Y,I) * rho_fibre * kron(Y,I)') +...
                (p_fibreZ1 * kron(Z,I) * rho_fibre * kron(Z,I)');

end