function [reward,total_time_ES] = Environment_QCNetwork(Path,f_attenuation,p_init,refractive_index,r_dephase,c_light,t_BSM,t_d,Distance_Adjacency_Matrix, p_BSM, p_GateErrors)
    
    total_time_ES = 0;
    Nodes_inPath = length(Path);
    ket0 = [1; 0];
    ket1 = [0; 1];
    EP_used = (kron(ket0,ket0)+kron(ket1,ket1))/sqrt(2);
    rho_EP_used = EP_used*EP_used';
    Node_Psi_MemoryTime(1:Nodes_inPath) = struct('NodeNo',0,'PsiLeft_DM_EP',rho_EP_used,'PsiRight_DM_EP',rho_EP_used,'MemoryLeftWaitTime',0,'MemoryRightWaitTime',0,'P_fibreLeft_EP',0,'P_fibreRight_EP',0,'P_bitflipLeft',0,'P_bitflipRight',0); 
    %PsiLeft and PsiRight of all nodes in the path are same initially.
    %Ideally, PsiLeft of Node 'a' is same as PsiRight of Node 'a+1'. Nodes 'a' and 'a+1' share same link. 
    %Ideally, PsiRight of Node 'a' is same as PsiLeft of Node 'a+1'. Nodes 'a' and 'a+1' share same link. 
    %Also, PsiLeft of first node (Tx) does not play any role. It is not required.
    for nn =1:Nodes_inPath
        Node_Psi_MemoryTime(nn).NodeNo = Path(nn);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Fiber Loss during ES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp_dist = [];
    t_q = [];                                   %Time taken by entangled qubits to reach consecutive nodes from entangled pair generator 
% %     Passed in the function call
% %     c_light = 3.0e8;                            % Speed of light
% %     t_BSM = 10e-9;                              %Time taken for BSM operation
% %     t_d = 10e-9;                                %Time taken for manipulation (gate) operations at end of link-level teleportation and entanglement swapping
    for ll = 1:length(Path)-1                   %or equivalent to entangled pair %doing link/channel-wise

        l_elem = Distance_Adjacency_Matrix(Path(ll),Path(ll+1));
        temp_dist = [temp_dist l_elem];
        t_q = [t_q (l_elem*1e3/2)/(c_light/refractive_index)];          %(l_elem*1e3) for writing 'l_elem Km'

        [DM_EP, P_fibre_EP] = FiberLoss(p_init,f_attenuation,l_elem/2,rho_EP_used);        

        %Left node of this channel/link ...... Right side memory of this left node
        Node_Psi_MemoryTime(ll).PsiRight_DM_EP = DM_EP;
        Node_Psi_MemoryTime(ll).P_fibreRight_EP = P_fibre_EP;
        
        %Right node of this channel/link...... Left side memory of this right node
        Node_Psi_MemoryTime(ll+1).PsiLeft_DM_EP = DM_EP;
        Node_Psi_MemoryTime(ll+1).P_fibreLeft_EP = P_fibre_EP; 
    end

    [L_max, L_max_Indx] = max(temp_dist);
    t_q_max = t_q(L_max_Indx);                        %%%% Time delay for qubit transfer between entangled pair generator and repeater
    t_q_diff = t_q_max.*ones(1,length(t_q)) - t_q;
    total_time_ES = total_time_ES + t_q_max;

    %%Update wait time due to difference in TOA of qubits from entangled pair generator at different nodes (This differs because of different link length between two nodes)
    for ll = 1:length(Path)-1
        Node_Psi_MemoryTime(ll).MemoryRightWaitTime = Node_Psi_MemoryTime(ll).MemoryRightWaitTime + t_q_diff(ll);
        Node_Psi_MemoryTime(ll+1).MemoryLeftWaitTime = Node_Psi_MemoryTime(ll+1).MemoryLeftWaitTime + t_q_diff(ll);
    end
    
    %%%%%%%% Memory decoherence (dephasing), imperfect BSM and gate operations-related noises/errors during ES %%%%%%%%%
    [Fidelity,total_time_ES] = EntanglementSwap_NoiseProbabilities(Node_Psi_MemoryTime,r_dephase,Distance_Adjacency_Matrix,t_BSM,c_light,t_d,rho_EP_used, p_BSM, p_GateErrors,total_time_ES);
  
    %%%%%%%%%%%%%%%%%% Calc. reward based on Fidelity %%%%%%%%%%%%%%%%%%%%
    if (Fidelity >= 0) && (Fidelity < 0.25) 
        reward = -2*( 0.5 - Fidelity ); 
    elseif (Fidelity >= 0.25) && (Fidelity < 0.5) 
        reward = -1*( 0.5 - Fidelity ); 
    elseif (Fidelity >= 0.5) && (Fidelity < 0.6) 
        reward = 1*Fidelity; 
    elseif (Fidelity >= 0.6) && (Fidelity < 0.7) 
        reward = 2*Fidelity; 
    elseif (Fidelity >= 0.7) && (Fidelity < 0.8) 
        reward = 3*Fidelity; 
    elseif (Fidelity >= 0.8) && (Fidelity < 0.9) 
        reward = 4*Fidelity; 
    elseif (Fidelity >= 0.9) && (Fidelity <= 1) 
        reward = 5*Fidelity; 
    end

end
