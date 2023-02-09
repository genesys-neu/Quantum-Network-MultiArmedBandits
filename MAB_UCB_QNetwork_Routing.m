function [Best_Path] = MAB_UCB_QNetwork_Routing(rounds,experiments,Feasible_paths_count,Feasible_paths,f_attenuation,p_init,refractive_index,r_dephase,c_light,t_BSM,t_d,Distance_Adjacency_Matrix, p_BSM, p_GateErrors)

    qlib; %Qlib dependecy
    
    No_of_arms = Feasible_paths_count;
    Arms_route = Feasible_paths;
%     rewards_Initial = Feasible_paths_Fidelity;

    Reward_Round_Exp = zeros(rounds,experiments);
    Reward_Arm_Round_Exp = zeros(No_of_arms,rounds,experiments);

    Pulls_Arm_Exp = zeros(No_of_arms,experiments);    
    Pulls_Arm_Round_Exp = zeros(No_of_arms,rounds,experiments);
    
    AverageReward_Arm_Exp = zeros(No_of_arms,experiments);
    AverageReward_Arm_Round_Exp = zeros(No_of_arms,rounds,experiments);
    
    total_time_ES_Round_Exp = zeros(rounds,experiments);
    
    for iter = 1:experiments
       iter
       
        reward_Round = [];
        pulls_Arm = []; 
        average_Arm = []; 
        Reward_Arm_Round = []; 
        Pulls_Arm_Round = []; 
        AverageReward_Arm_Round = [];
        total_time_ES_Round = [];
        [reward_Round, pulls_Arm, average_Arm, Reward_Arm_Round, Pulls_Arm_Round, AverageReward_Arm_Round, total_time_ES_Round] = UCB(No_of_arms,Arms_route,rounds,f_attenuation,p_init,refractive_index,r_dephase,c_light,t_BSM,t_d,Distance_Adjacency_Matrix, p_BSM, p_GateErrors);

        Reward_Round_Exp(:,iter) = reward_Round;
        Reward_Arm_Round_Exp(:,:,iter) = Reward_Arm_Round;
    
        Pulls_Arm_Exp(:,iter) = pulls_Arm;    
        Pulls_Arm_Round_Exp(:,:,iter) = Pulls_Arm_Round;
        
        AverageReward_Arm_Exp(:,iter) = average_Arm;
        AverageReward_Arm_Round_Exp(:,:,iter) = AverageReward_Arm_Round;
    
        total_time_ES_Round_Exp(:,iter) = total_time_ES_Round;
    end
    
%     filename = 'D:\Quantum Comput n Comm\Multi-Arm-Bandit-Simulation-master\Sept9_Temp\S2D15_UCBE600R1400I1_MemDeph_BSMX_ImperfGateY_FiberDepolModel_c40.mat';
%     save(filename);
% 
%     figure
%     plot(mean(Reward_Round_Exp,2),'b');
%     ylabel('Reward averaged over Experiments');
%     xlabel('Rounds (or Time Steps)');
%     title('UCB');
% 
%     figure
%     stem(mean(Pulls_Arm_Exp,2),'b');
%     ylabel('Pull averaged over Experiments');
%     xlabel('Arms');
%     title('UCB');
% 
%     figure
%     Pulls_Arm_Round_AvgExp = mean(Pulls_Arm_Round_Exp,3);       %Arms X Rounds
%     for aa = 1:No_of_arms
%         plot(Pulls_Arm_Round_AvgExp(aa,:)./experiments);
%         if aa ~= No_of_arms
%             hold on
%         end
%     end
%     ylabel('Percentage of Pulls averaged over Experiments');
%     xlabel('Rounds');
%     
%     figure
%     for aa = 1:No_of_arms
%         subplot(ceil(No_of_arms/2),2,aa)
%         histogram(AverageReward_Arm_Exp(aa,:))
%         ylabel('Prob. Distr. average rewards in diff. Exp.');
%         xlabel('Average Reward');
%     end

    %%%%%%%%%%%%%%%% Obtain the least noisy path %%%%%%%%%%%%%%%%%%%%%%% 
    Pathwise_Pull = mean(Pulls_Arm_Exp,2);
    [Times_BestPathSelected Best_Path_Indx] = max(Pathwise_Pull);
    Best_Path = Arms_route{Best_Path_Indx};
end

function [reward, pulls, average, Reward_Arm_Round, Pulls_Arm_Round, AverageReward_Arm_Round, total_time_ES_Round]= UCB(No_of_arms,Arms_route,rounds,f_attenuation,p_init,refractive_index,r_dephase,c_light,t_BSM,t_d,Distance_Adjacency_Matrix, p_BSM, p_GateErrors)
    %     Input :
    %         No_of_arms: No of feasible paths
    %         Arms_route: Feasible paths
    %         rounds: Iterations (used for exploration and exploitation)
    %     Output: 
    %         reward: reward in each round (or iteration)
    %         pulls: count pulls for each arm after last iteration
    %         average: average reward armwise after last iteration
    %         Reward_Arm_Round: Stores rewars of each arm in each round
    %         Pulls_Arm_Round: Stores count of pulls of each arm in each round
    %         AverageReward_Arm_Round: average reward armwise and roundwise
    %         total_time_ES_Round: Stores total time taken in each round
    

    pulls = zeros(No_of_arms,1);          % initialize
    reward = zeros(1,rounds);     
    average = zeros(No_of_arms,1);
    AverageReward_Arm_Round = zeros(No_of_arms,rounds);
    Reward_Arm_Round = zeros(No_of_arms,rounds);
    Pulls_Arm_Round = zeros(No_of_arms,rounds);
    total_time_ES_Round = zeros(1,rounds);
    
    % Play each machine once
    for iter = 1:No_of_arms
        index = iter;
        [result, total_time_ES] = Environment_QCNetwork(Arms_route{index},f_attenuation,p_init,refractive_index,r_dephase,c_light,t_BSM,t_d,Distance_Adjacency_Matrix, p_BSM, p_GateErrors);       %normrnd(arms(index,1),std);
        average(index,1) = (average(index,1)*pulls(index,1)+ result)/(pulls(index,1)+1);                        %This is Qvalue calculation only
%             average(index,1) = average(index,1) + ( ( 1/( pulls(index,1)+1 ) )*(result - average(index,1) ) );      %Q value implementation (another form)
        reward(1,iter) = result;        
        pulls(index,1) = pulls(index,1) + 1;
        total_time_ES_Round(1,iter) = total_time_ES;
        for aa = 1:No_of_arms
            AverageReward_Arm_Round(aa,iter) = average(aa,1);
            Pulls_Arm_Round(aa,iter) = pulls(aa,1);
            if aa == index
                Reward_Arm_Round(aa,iter) = reward(1,iter);
            else
                Reward_Arm_Round(aa,iter) = 0;
            end
        end        
    end

    % play machine j that maximize UCB rules
    for iter = (No_of_arms+1):rounds
        % select the arm to pull
        temp = zeros(No_of_arms,1);
        idx=1;
        m=0;
        for j=1:No_of_arms
            temp(j) = average(j,1) + (4)*sqrt(2*log(iter)/pulls(j,1));
            if m < temp(j)
                m = temp(j);
                idx = j;
            end
        end
        
        % pull the arms with idx
        [res, total_time_ES] = Environment_QCNetwork(Arms_route{idx},f_attenuation,p_init,refractive_index,r_dephase,c_light,t_BSM,t_d,Distance_Adjacency_Matrix, p_BSM, p_GateErrors);          %normrnd(arms(idx,1),std);
        average(idx,1) = (average(idx,1)*pulls(idx,1)+ res)/(pulls(idx,1)+1);
        reward(1,iter) =  res;
        pulls(idx,1) = pulls(idx,1) + 1;
        total_time_ES_Round(1,iter) = total_time_ES;
        for aa = 1:No_of_arms
            AverageReward_Arm_Round(aa,iter) = average(aa,1);
            Pulls_Arm_Round(aa,iter) = pulls(aa,1);
            if aa == idx
                Reward_Arm_Round(aa,iter) = reward(1,iter);
            else
                Reward_Arm_Round(aa,iter) = 0;
            end
        end

    end             %end of for (rounds)

end                 %UCB function end
