% Quantum Entanglement Swapping with noises/errors encountered by qubits in memory, during BSM, and gate operations.
% Scenario: Three quantum nodes, two entangled pairs (q1-q2, q3-q4) shared b/w the three nodes- Node1 has q1, node2 has q2 and q3, node3 has q4. 
% BSM is carried out on q2 and q3 at node2 and measurement results are sent to node3 to perform gate operations on q4 so that q1-q4 becomes entangled.
% Types of noise/error considered: Dephasing noise in memory, bit-flip error during BSM, bit-phase flip errorin gate operations.

function [ finalQBitDensityMatrix, avg_fidelity ] = EntanglementSwap_Noises(qBitToTeleportDensityMatrix,maxEntangledState_DM,Perr_q2,Perr_q3,Perr_q4, p_BSM, p_GateErrors)

%     qlib; %Qlib dependecy

%     Constants
    Id = [1 0; 0 1];
    ket0 = [1; 0];
    ket1 = [0; 1];
    CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];
    Hadamard = (1/sqrt(2))*[1 1; 1 -1];
    Xgate = [0 1;1 0];
    Zgate = [1,0;0,-1];

%    Density Matrix of the system    
    systemDensityMatrix = kron(qBitToTeleportDensityMatrix,maxEntangledState_DM);
        
    NoiseFlag = 2;

%    Applying noise
    systemDensityMatrix_orig = systemDensityMatrix;
    fidelity = [];
    Iterations = 1;     %Generation of noise instances and not the average behavior
    
    for nn = 1:Iterations
        
        systemDensityMatrix = systemDensityMatrix_orig;
        %Noise corresponding to decoherence in memory of qubits 2 and 3 (Dephasing noise)
        if (NoiseFlag == 2)    
        
            prob_error_q2 = Perr_q2;
            q2_error_status = (rand(1) <= prob_error_q2);   %1/0 = phase-flip/No phase-flip
            NoiseGate_q2 = ( (1-q2_error_status).*Id ) + ( (q2_error_status).*Zgate );

            prob_error_q3 = Perr_q3;
            q3_error_status = (rand(1) <= prob_error_q3);   %1/0 = phase-flip/No phase-flip
            NoiseGate_q3 = ( (1-q3_error_status).*Id ) + ( (q3_error_status).*Zgate );
            
            % Any operation on q1 and q4 are not required so far, that's why not considered noise/error in them. Considered memory decoherence of q1 and q4 in later part of code. 
            NoiseGate_q1 = ( (1-0).*Id ) + ( (0).*Zgate );
            NoiseGate_q4 = ( (1-0).*Id ) + ( (0).*Zgate );
                        
            CombinedNoise_q1q2q3q4 = kron(kron(kron(NoiseGate_q1,NoiseGate_q2),NoiseGate_q3),NoiseGate_q4);
            systemDensityMatrix = CombinedNoise_q1q2q3q4 * systemDensityMatrix * CombinedNoise_q1q2q3q4';
        end
        
        %BSM process starts
        % Applying CNOT Gate
        CNOT_q2q3 = kron(Id,kron(CNOT, Id));
        systemDensityMatrix = CNOT_q2q3 * systemDensityMatrix * CNOT_q2q3';
        
        BSM_error_status = (rand(1) <= p_BSM);   %1/0 = bit-flip/No bit-flip
        Hadamard_q2 = ( (1-BSM_error_status).*Hadamard ) + ( (BSM_error_status).*(Xgate*Hadamard) );            %BSM errors are bit-flip errors     

        H = kron(Id,kron(kron(Hadamard_q2,Id),Id));
        systemDensityMatrix = H * systemDensityMatrix * H';
                
        % Measuring the second qubit q2
        systemDensityMatrix = measureSingleQBit(systemDensityMatrix, [0 1 0 0]);
        
        % Measuring the third qubit q3
        systemDensityMatrix = measureSingleQBit(systemDensityMatrix, [0 0 1 0]);
        
        qBitToTeleportAfterMeasurementDensityMatrix = partial_trace(systemDensityMatrix, [0 1 0 0]);
        qBitToTeleportAfterMeasurement = dm2pure(qBitToTeleportAfterMeasurementDensityMatrix);
        
        channelQBitAfterMeasurementDensityMatrix = partial_trace(systemDensityMatrix, [0 0 1 0]);
        channelQBitAfterMeasurement = dm2pure(channelQBitAfterMeasurementDensityMatrix);
        
        teleportedQBitDensityMatrix = partial_trace(systemDensityMatrix, [1 0 0 1]);
        teleportedQBit = dm2pure(teleportedQBitDensityMatrix);
        
        % Noise on q4 before manipulation operations/gate operations post q2 qnd q3 measurements
        if (NoiseFlag == 2)
        
            % Decoherence of q4 in memory (dephasing) before applying gate operations 
            prob_error_q4 = Perr_q4;
            q4_error_status = (rand(1) <= prob_error_q4);   %1/0 = phase-Flip/No phase-flip
            NoiseGate_q4 = ( (1-q4_error_status).*Id ) + ( (q4_error_status).*Zgate );
        
            %Any operation on q1 is not required in time duration from BSM of q2,q3 to manipulation of q4, that's why not considered error in q1. 
            %Memory decoherence of q1 is calculated when it will be involved in operations in the end-to-end entanglement distribution over a path.  
            %Since q2 and q3 are already collapsed due to mmt, that's why not considered here.
            NoiseGate_q1 = ( (1-0).*Id ) + ( (0).*Zgate );
            
            CombinedNoise_q1q4 = kron(NoiseGate_q1,NoiseGate_q4);
            teleportedQBitDensityMatrix = CombinedNoise_q1q4 * teleportedQBitDensityMatrix * CombinedNoise_q1q4';
        end
        
        % Density matrix of entangled q1-q4 after swapping entanglement of q1-q2 and q3-q4. 
        finalQBitDensityMatrix = operationAfterMeasure_ES_GateErr(qBitToTeleportAfterMeasurement, channelQBitAfterMeasurement, teleportedQBitDensityMatrix, p_GateErrors);
        
        if (NoiseFlag == 2)
            fidelity = [fidelity trace(qBitToTeleportDensityMatrix * finalQBitDensityMatrix)];
        end
    
    end
    avg_fidelity = mean(fidelity);      %equals to fidelity because iterations = 1.

end    