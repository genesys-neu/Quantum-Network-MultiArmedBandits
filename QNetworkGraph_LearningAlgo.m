clear all
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%% Quantum network Deployment parameters %%%%%%%%%%%%%%%%%%%%%%%     
%%Network topology creation and edge fidelity calculation
%%Distance parameters
No_of_nodes = 15;
x_min = -10;                        %% wrt origin (0,0)
x_max = 10;                         %% wrt origin (0,0)
Link_distance_Threshold = 10;       %%In Kms
     
%%Assigning x,y coordinates to nodes
x_loc = zeros(No_of_nodes,1);
y_loc = zeros(No_of_nodes,1);
ConnectedGraph_Flag = 0;
loop_count = 0;

%%%% Generation of a fully connected graph
%%%% Considered source and destination nodes (configurable)
sources = 2;
destination = 15; 
while (ConnectedGraph_Flag == 0)
    %%Random deployment of nodes (Comment lines 29 and 30 if want to generate a new connected network graph)
%     x_loc = ( (x_max-x_min).*rand(No_of_nodes,1) ) + ( x_min.*ones(No_of_nodes,1) );
%     y_loc = ( (x_max-x_min).*rand(No_of_nodes,1) ) + ( x_min.*ones(No_of_nodes,1) );

%     %%Use Sample connected graph (Comment lines 25 and 26)
    x_loc = [-9.20000000000000, -6.50000000000000, -7.50000000000000, -4, -5.10000000000000, -0.200000000000000, 0.100000000000000, -1.20000000000000, -1, 4, 3, 5, 6, 7.50000000000000, 8.20000000000000];
    y_loc = [-3.20000000000000, -7, 4.50000000000000, 8, -0.200000000000000, -8.60000000000000, 7.60000000000000, 3.20000000000000, -3, 4.30000000000000, -6, -2, 7, -6, 2];

    %%Link establishment
    Distance_Matrix = zeros(No_of_nodes,No_of_nodes);
    Distance_Adjacency_Matrix = zeros(No_of_nodes,No_of_nodes);
    Adjacency_Matrix = zeros(No_of_nodes,No_of_nodes);
    for nn = 1:No_of_nodes
        for mm = 1:No_of_nodes
            if (nn ~= mm)
                Distance_Matrix(nn,mm) = sqrt( ((x_loc(nn)-x_loc(mm))^2) + ((y_loc(nn)-y_loc(mm))^2) );
                if Distance_Matrix(nn,mm) < Link_distance_Threshold                                       %%Link will exist b/w two adjacent nodes if they are less than 2 Km apart.
                    %%Distance Adjacency Matrix
                    Distance_Adjacency_Matrix(nn,mm) = Distance_Matrix(nn,mm);
                    Adjacency_Matrix(nn,mm) = 1;                %%Will beused in finding if there is atleast 1 route b/w each pair of nodes
                else
                    %%Distance Adjacency Matrix
                    Distance_Adjacency_Matrix(nn,mm) = inf;
                    Adjacency_Matrix(nn,mm) = 0;                %%Will beused in finding if there is atleast 1 route b/w each pair of nodes
                end
            else
                Distance_Adjacency_Matrix(nn,mm) = inf;         %either no link b/w 2 nodes or no self loop
                Adjacency_Matrix(nn,mm) = 0;                %%Will beused in finding if there is atleast 1 route b/w each pair of nodes
            end
        end
    end
    
    %%Create graph with default numbering/naming for vertices
    NetworkGraph = graph(Adjacency_Matrix); 
    Node_Pairs = nchoosek(1:No_of_nodes,2);
    No_of_Paths_Node_Pairs = zeros(size(Node_Pairs,1),1);
    ConnectedGraph_Flag = 1;
    for np = 1:size(Node_Pairs,1)
        paths = allpaths(NetworkGraph,Node_Pairs(np,1),Node_Pairs(np,2));      
        No_of_Paths_Node_Pairs(np) = length(paths);       %%No of paths for each node pair....used to iteratively call our algorithm until Fidelity threshold is violated
        if isempty(paths) == 1
            ConnectedGraph_Flag = 0;
            break;
        elseif ( (Node_Pairs(np,1)==sources) && (Node_Pairs(np,2)==destination) ) || ( (Node_Pairs(np,2)==sources) && (Node_Pairs(np,1)==destination) ) 
            AllPaths_S_D = paths;      %save all paths b/w source and destination
            S_D_Indx_in_NodePairs = np;
        end
    end
    loop_count = loop_count+1;
end     %end while

figure
plot(x_loc,y_loc,'*');
grid on

%%%% Quantum Network Settings (Customizable)
p_init = 0.00001;                           % Probability of loss of entangled pair after generation
p_BSM = 0.2;                                % Probability of imperfect (or error in) BSM operation
p_GateErrors = 0.2;                         % Probability of occurrence of gate errors during swapping operations
r_dephase = 10000;                          % Dephasing rate of quantum memories
c_light = 3.0e8;                            % Speed of light
t_BSM = 10e-9;                              % Time taken for BSM operation
t_d = 10e-9;                                % Time taken for manipulation (gate) operations at end of link-level teleportation and entanglement swapping
f_attenuation = 0.05;                       % Fibre Loss attenuation db/km
refractive_index = 1.5;                     % speed of light in fiber = c_light/refractive_index
 
%%%% Find initial (feasible) paths b/w source and node based on Fiber loss only. These will be set as arms in Bandit learning algorithm
%%%% One can build their own logic to choose paths for setting arms in Bandit algorithm
%%Link level Fidelity calculation
EdgeFidelityMatrix = zeros(No_of_nodes,No_of_nodes);
for nn = 1:No_of_nodes
    for mm = 1:No_of_nodes
        if (nn ~= mm) && (Distance_Matrix(nn,mm) < Link_distance_Threshold)                                       %%Link will exist b/w two adjacent nodes if they are less than 2 Km apart.
                EdgeFidelityMatrix(nn,mm) = Link_EntangledPair_Fidelity(p_init,f_attenuation,Distance_Matrix(nn,mm)/2);      
        end
    end
end
netCostMatrix = 1./EdgeFidelityMatrix;
Fidelity_threshold_EachFlow = 0.5; 

source = sources(1);
No_of_arms = 8;     
max_paths = No_of_arms*1;
[Feasible_paths_count,Feasible_paths,Feasible_paths_Fidelity] = Fidelity_FeasiblePaths(netCostMatrix, source, destination, Fidelity_threshold_EachFlow,max_paths);

%%%%% call to MAB Learning Algorithm for Route Selection
rounds = 1200;      %Customizable                  
experiments = 100;  %Customizable
[Best_Path] = MAB_UCB_QNetwork_Routing(rounds,experiments,Feasible_paths_count,Feasible_paths,f_attenuation,p_init,refractive_index,r_dephase,c_light,t_BSM,t_d,Distance_Adjacency_Matrix, p_BSM, p_GateErrors);       %r_dephase for memory decoherence, r_depolarize for fiber loss, p_BSM and p_GateErrors for imperfect gates
