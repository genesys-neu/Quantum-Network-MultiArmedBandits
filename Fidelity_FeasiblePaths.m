function [Feasible_paths_count,Feasible_paths,Feasible_paths_Fidelity] = Fidelity_FeasiblePaths(netCostMatrix, source, destination, Fidelity_Threshold, max_Paths)
    
    Fidelity_paths = [];
    k_paths_loop = max_Paths;
    min_achieved_Fidelity = 1;
    [shortestPaths, totalCosts] = kShortestPath_QNetwork(netCostMatrix, source, destination, k_paths_loop);
    Fidelity_paths = 1./totalCosts;
    [min_achieved_Fidelity] = min(Fidelity_paths);
    while (min_achieved_Fidelity < Fidelity_Threshold)       
        k_paths_loop = k_paths_loop - 1;
        [min_achieved_Fidelity] = min(Fidelity_paths(1:k_paths_loop));
    end
    Feasible_paths_count = k_paths_loop;                                
    Feasible_paths = shortestPaths(1:Feasible_paths_count);             
    Feasible_paths_Fidelity = Fidelity_paths(1:Feasible_paths_count);
    
end