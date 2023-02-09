function [F_elem] = Link_EntangledPair_Fidelity(p_init,f_attenuation,l_elem)

%%%% In this Code, we assume all qubit are initially pure sate 
%%%% Fibre Loss 
    p_fibre = ( 1-(1-p_init) * (10^(-f_attenuation * (l_elem) / 10)) );   %% Fibre Loss Probability
    
    ket0 = [1; 0];
    ket1 = [0; 1];
    CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];
    Hadamard = (1/sqrt(2))*[1 1; 1 -1];
    Xgate = [0 1;1 0];
    
    maxEntangledState = (kron(ket0,ket0)+kron(ket1,ket1))/sqrt(2);
    
    %%%% Operators
    I = [1,0;0,1]; %%%% I Operator
    X = [0,1;1,0]; %%%% X Operatr
    Y = [0,-1i;1i,0]; %%%% Y Operator
    Z = [1,0;0,-1];   %%%% Z Operator

    %%%% Density Matrix of entangled qubits
    rho_init = maxEntangledState*maxEntangledState'; %%%% Initial Density Matrix of one qubit
    
    p_fibreX1 = (1/3) * p_fibre;
    p_fibreY1 = (1/3) * p_fibre;
    p_fibreZ1 = (1/3) * p_fibre;
   
    p_fibreX2 = (1/3) * p_fibre;
    p_fibreY2 = (1/3) * p_fibre;
    p_fibreZ2 = (1/3) * p_fibre;

    rho_fibre = ((1 - p_fibreX2 - p_fibreY2 - p_fibreZ2) * kron(I,I) * rho_init * kron(I,I)') +...
                (p_fibreX2 * kron(I,X) * rho_init * kron(I,X)') +...
                (p_fibreY2 * kron(I,Y) * rho_init * kron(I,Y)') +...
                (p_fibreZ2 * kron(I,Z) * rho_init * kron(I,Z)');
                
    rho_fibre = ((1 - p_fibreX1 - p_fibreY1 - p_fibreZ1) * kron(I,I) * rho_fibre * kron(I,I)') +...
                (p_fibreX1 * kron(X,I) * rho_fibre * kron(X,I)') +...
                (p_fibreY1 * kron(Y,I) * rho_fibre * kron(Y,I)') +...
                (p_fibreZ1 * kron(Z,I) * rho_fibre * kron(Z,I)');

    %%%% Fidelity calculation of one qubit one elementary link
    F_elem = trace(rho_init * rho_fibre);
    
end
