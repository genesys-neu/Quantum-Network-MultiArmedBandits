function [Fidelity,total_time_ES] = EntanglementSwap_NoiseProbabilities(Node_Psi_MemoryTime,r_dephase,Distance_Adjacency_Matrix,t_BSM,c_light,t_d,rho_EP_used, p_BSM, p_GateErrors,total_time_ES)
    
    % % %% Two Entangled pairs (q1-q2, q3-q4) shared b/w three nodes. 
    % % Node1 has q1, node2 has q2 and q3, node3 has q4. Representation 1-2    2-3
    % % 
    % % %%Entanglement Swapping
    % %                                             Tx 1, Rx 7, Odd no. of nodes in path
    % % 1-2   2-3   3-4   4-5   5-6   6-7           
    % %                                             BSM on 2-2, 4-4, 6-6.........ES: 1-2-2-3, 3-4-4-5, and 5-6-6-7
    % %    1-3        3-5         5-7
    % %                                             BSM on 3-3.................ES: 1-3-3-5 
    % %         1-5               5-7
    % %                                             BSM on 5-5.................ES: 1-5-5-7
    % %                1-7
    % % 
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % % 
    % %                                            Tx 1, Rx 7, Even no. of nodes in path
    % % 1-2   2-3   3-4   4-5   5-6   6-7   7-8
    % %                                            BSM on 2-2,4-4 and 6-6.........ES: 1-2-2-3, 3-4-4-5, 5-6-6-7  
    % %    1-3         3-5        5-7       7-8
    % %                                            BSM on 3-3, 7-7.................ES: 1-3-3-5 and 5-7-7-8
    % %          1-5                  5-8
    % %                                            BSM on 5-5.................ES: 1-5-5-8
    % %                    1-8

    % Structure of 'Node_Psi_MemoryTime(1:Nodes_inPath)' for reference.
    % Node_Psi_MemoryTime(1:Nodes_inPath) = struct('NodeNo',0,'PsiLeft_DM_EP',rho_EP_used,'PsiRight_DM_EP',rho_EP_used,'MemoryLeftWaitTime',0,'MemoryRightWaitTime',0,'P_fibreLeft_EP',0,'P_fibreRight_EP',0,'P_bitflipLeft',0,'P_bitflipRight',0); 

    %%%%Entanglement Swapping and Noise generation
    TotalIter_Err = 1;  
    Node_Psi_MemoryTime_Backup = Node_Psi_MemoryTime;
    Fidelity_Path = zeros(TotalIter_Err,1);
    for iter1 = 1:TotalIter_Err
        
        Node_Psi_MemoryTime = Node_Psi_MemoryTime_Backup;
        NodesIdx_unused = 1:length(Node_Psi_MemoryTime);
        while length(NodesIdx_unused) > 2      
    
            tc_ES_LevelWise_Left = zeros(length(Node_Psi_MemoryTime),1);
            tc_ES_LevelWise_Right = zeros(length(Node_Psi_MemoryTime),1);
            Nodes_unused_updated = NodesIdx_unused;
            
            for nn = 2:2:length(NodesIdx_unused)-1       
                
                %For q2 in function 'EntanglementSwap_Noises'
                Current_NodeIndx = NodesIdx_unused(nn);
                Prev_NodeIndx = NodesIdx_unused(nn-1); 
                Next_NodeIndx = NodesIdx_unused(nn+1);
                t_q_left = Node_Psi_MemoryTime(Current_NodeIndx).MemoryLeftWaitTime;
                Node_Psi_MemoryTime(Current_NodeIndx).P_bitflipLeft = (1-exp(-t_q_left * r_dephase));            %p_depolarize = (1-exp(-t_depolarize * r_depolarize));  %%% Depolarize Probability, unit is [1];
                EP_ToTeleport_DM = Node_Psi_MemoryTime(Current_NodeIndx).PsiLeft_DM_EP;
                
                %For q3 in function 'EntanglementSwap_Noises'
                t_q_right = Node_Psi_MemoryTime(Current_NodeIndx).MemoryRightWaitTime;
                Node_Psi_MemoryTime(Current_NodeIndx).P_bitflipRight = (1-exp(-t_q_right * r_dephase));
                EP_Channel_DM = Node_Psi_MemoryTime(Current_NodeIndx).PsiRight_DM_EP;
    
                %For q4, t wait = t_q + t_c + t_BSM
                classical_dis = Distance_Adjacency_Matrix(Node_Psi_MemoryTime(Current_NodeIndx).NodeNo,Node_Psi_MemoryTime(Next_NodeIndx).NodeNo);
                t_c = (classical_dis)*(1e3)/c_light;                            %%%% Time delay for classical information transmission
                tc_ES_LevelWise_Left(Next_NodeIndx) = t_c;
                Node_Psi_MemoryTime(Next_NodeIndx).MemoryLeftWaitTime = Node_Psi_MemoryTime(Next_NodeIndx).MemoryLeftWaitTime + t_c + t_BSM;
                t_wait = Node_Psi_MemoryTime(Next_NodeIndx).MemoryLeftWaitTime;
                Node_Psi_MemoryTime(Next_NodeIndx).P_bitflipLeft = (1-exp(-t_wait * r_dephase));      %prob error for q4 in fn ES_ProbErr
                
                %Call function 'EntanglementSwap_Noises'
                [ finalQBitDensityMatrix, avg_fidelity ] = EntanglementSwap_Noises(EP_ToTeleport_DM,EP_Channel_DM,Node_Psi_MemoryTime(Current_NodeIndx).P_bitflipLeft,Node_Psi_MemoryTime(Current_NodeIndx).P_bitflipRight,Node_Psi_MemoryTime(Next_NodeIndx).P_bitflipLeft, p_BSM, p_GateErrors);    
                %Node_Psi_MemoryTime(nn).P_bitflipRight will be for q3, while Node_Psi_MemoryTime(nn).P_bitflipLeft for q2, Node_Psi_MemoryTime(nn+1).P_bitflipLeft for q4
                
                %Update density matrix for EP (obtained from Eswapping) on link (nn-1,nn+1)
                %Left node of this channel/link ...... Right side memory of this left node
                Node_Psi_MemoryTime(Prev_NodeIndx).PsiRight_DM_EP = finalQBitDensityMatrix;
                %Right node of this channel/link...... Left side memory of this right node
                Node_Psi_MemoryTime(Next_NodeIndx).PsiLeft_DM_EP = finalQBitDensityMatrix;
    
                %For q2, q3, and q4, update/reset memory wait times 
                Node_Psi_MemoryTime(Current_NodeIndx).MemoryLeftWaitTime = 0;
                Node_Psi_MemoryTime(Current_NodeIndx).MemoryRightWaitTime = 0;
                Node_Psi_MemoryTime(Next_NodeIndx).MemoryLeftWaitTime = 0;
                
                %For q1, update t_wait = 
                tc_ES_LevelWise_Right(Prev_NodeIndx) = t_c;
                Node_Psi_MemoryTime(Prev_NodeIndx).MemoryRightWaitTime = Node_Psi_MemoryTime(Prev_NodeIndx).MemoryRightWaitTime + t_c + t_BSM + t_d;
    
                %Update status of used/unused nodes
                Nodes_usedFlag(Current_NodeIndx) = 1;
                Nodes_unused_updated = setdiff(Nodes_unused_updated,Current_NodeIndx);
            end
            
            %%Update wait times of unused nodes (i.e., the nodes on which no operation is performed)
            if mod(length(NodesIdx_unused),2) == 0
                %%Pick the last node because that is the only unused one
                mm = NodesIdx_unused(end);
                Node_Psi_MemoryTime(mm-1).MemoryRightWaitTime = Node_Psi_MemoryTime(mm-1).MemoryRightWaitTime + max(tc_ES_LevelWise_Right) + t_BSM + t_d;
                tc_ES_LevelWise_Right(mm-1) = max(tc_ES_LevelWise_Right);
                Node_Psi_MemoryTime(mm).MemoryLeftWaitTime = Node_Psi_MemoryTime(mm).MemoryLeftWaitTime + max(tc_ES_LevelWise_Left) + t_BSM + t_d;
                tc_ES_LevelWise_Left(mm-1) = max(tc_ES_LevelWise_Left);
            end
    
            %%Update wait times of all nodes based on the one that took max time in classical signalling
            for nn = 1:1:length(NodesIdx_unused)
                max_tc = max([tc_ES_LevelWise_Left;tc_ES_LevelWise_Right]);
                NodeIndx_to_be_chkd = NodesIdx_unused(nn);
                if isempty(intersect(Nodes_unused_updated,NodeIndx_to_be_chkd)) ~= 1
                    if NodeIndx_to_be_chkd ~= 1
                    Node_Psi_MemoryTime(NodeIndx_to_be_chkd).MemoryLeftWaitTime = Node_Psi_MemoryTime(NodeIndx_to_be_chkd).MemoryLeftWaitTime + (max_tc - tc_ES_LevelWise_Left(NodeIndx_to_be_chkd));
                    end
                    if NodeIndx_to_be_chkd ~= length(Node_Psi_MemoryTime)
                    Node_Psi_MemoryTime(NodeIndx_to_be_chkd).MemoryRightWaitTime = Node_Psi_MemoryTime(NodeIndx_to_be_chkd).MemoryRightWaitTime + (max_tc - tc_ES_LevelWise_Right(NodeIndx_to_be_chkd));
                    end
                end
            end
            
            total_time_ES = total_time_ES + max_tc + t_BSM + t_d;

            %Update unused nodes for next level of ES
            NodesIdx_unused = Nodes_unused_updated;
    
        end %end of entanglement swapping while loop

%         finalQBitDensityMatrix      %%This will represent system Density matrix of enatngled qubits shared b/w TX and Rx.
        Fidelity_Path(iter1) = trace(rho_EP_used * finalQBitDensityMatrix);
  
    end %end of iteration (used due to error creation)

    %Calculate reward
    Fidelity = mean(Fidelity_Path);      %Reward of this path

end %end of function