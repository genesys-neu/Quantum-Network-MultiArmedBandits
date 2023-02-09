function [ qBit ] = operationAfterMeasure_ES_GateErr( firstQBit, secondQBit, teleportedQBitDensityMatrix, p_GateErrors )
%
% Function that applies an operation on the teleported/entangle swapped QBits according to the two measured qBits
% 
% @returns array when success, -1 when failure
%

    % Constants
    ket0 = [1; 0];
    ket1 = [0; 1];
    sigmaX = [0 1; 1 0];
    sigmaY = [0,-1i;1i,0];
    sigmaZ = [1 0;0 -1];
    Id = [1 0; 0 1];
    
    qBit = -1;

    Gate_error_status = (rand(1) <= p_GateErrors);   %1/0 = Flip/No flip
    if rand(1) < 0.5
%         del_x = (-0.5*pi/180) + (pi/180)*rand(1);         %U(-0.5*pi/180,0.5*pi/180)
%         theta_x = pi + del_x;
%         Rx_theta_x = [cos(theta_x/2)    -1i*sin(theta_x/2);
%                      -1i*sin(theta_x/2)  cos(theta_x/2)];
%         sigmaX_General = ( (1-Gate_error_status).*sigmaX ) + ( (Gate_error_status).*Rx_theta_x );
        sigmaX_General = ( (1-Gate_error_status).*sigmaX ) + ( (Gate_error_status).*(sigmaY*sigmaX) );          %Gate errors are bitphaseflip errors
        sigmaZ_General = sigmaZ;
    else
%         del_z = (-0.5*pi/180) + (pi/180)*rand(1);         %U(-0.5*pi/180,0.5*pi/180)
%         theta_z = pi + del_z;
%         Rz_theta_z = [cos(theta_z/2)-(1i*sin(theta_z/2)) 0;
%                       0                                 cos(theta_z/2)+(1i*sin(theta_z/2))];
%         sigmaZ_General = ( (1-Gate_error_status).*sigmaZ ) + ( (Gate_error_status).*Rz_theta_z );
        sigmaZ_General = ( (1-Gate_error_status).*sigmaZ ) + ( (Gate_error_status).*(sigmaY*sigmaZ) );
        sigmaX_General = sigmaX;
    end
        
    if (firstQBit == ket0)
        if (secondQBit == ket0)
           qBit = teleportedQBitDensityMatrix;
        elseif (secondQBit == ket1)
%             qBit = kron(Id,sigmaX) * teleportedQBitDensityMatrix * (kron(Id,sigmaX))';
            qBit = kron(Id,sigmaX_General) * teleportedQBitDensityMatrix * (kron(Id,sigmaX_General))';
        end
    elseif (firstQBit == ket1)
        if (secondQBit == ket0)
%            qBit = kron(Id,sigmaZ) * teleportedQBitDensityMatrix * (kron(Id,sigmaZ))';
            qBit = kron(Id,sigmaZ_General) * teleportedQBitDensityMatrix * (kron(Id,sigmaZ_General))';
        elseif (secondQBit == ket1)
%             qBit = (kron(Id,sigmaX)) * ( (kron(Id,sigmaZ)) * teleportedQBitDensityMatrix * (kron(Id,sigmaZ))' ) * (kron(Id,sigmaX))';
            qBit = (kron(Id,sigmaX_General)) * ( (kron(Id,sigmaZ_General)) * teleportedQBitDensityMatrix * (kron(Id,sigmaZ_General))' ) * (kron(Id,sigmaX_General))';
        end
    end

end

