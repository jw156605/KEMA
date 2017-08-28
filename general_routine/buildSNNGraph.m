%-------------------------------------------------------------------
% Build a shared nearest neighbor graph
%
% [G,ICA] = buildSNNGraph(data,num_ICs,k_param,k_scale,prune_snn);
%
% Input:
%    - data     :      the input data matrix
%    - k_param  :      number of neighbors
%    - k_scale  :      number of neighbors
%    - prune_snn:	 pruning parameter between 0 and 1
%                                     
%
% Output:
%    - G    :      the graph
%    - ICA  :      ICA projection used to compute the SNN graph
%
% Joshua Welch
% 
% jwelch@cs.unc.edu
%  
%-------------------------------------------------------------------

function G = buildSNNGraph(data,num_ICs,k_param,k_scale,prune_snn)

csvwrite('temp.csv',data);
[res] = system(sprintf('Rscript /broad/macosko/jwelch/CellIntegration/KEMA/general_routine/build_SNN_graph.R temp.csv graph.mm ica.mm %i %i %i %f',num_ICs,k_param,k_scale,prune_snn),'-echo');
[G] = mmread('graph.mm');
[ICA] = mmread('ica.mm');
[res] = system('rm temp.csv graph.mm ica.mm');

end