%%%%%%%%%%%%% Cascaded Fidelity %%%%%%%%%%%%%%%%%%%%%%
function [F12] = CascadedFidelity(F1,F2)
    if (F1~=0)&&(F2~=0)
        F12 = (F1*F2) + ((1-F1)*(1-F2)/3);
    else    %case when both/either fidelities are zero (as if no path between nodes)
        F12 = 0;
    end
end