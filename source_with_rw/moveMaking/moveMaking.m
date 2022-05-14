function [ ncut_energies, labels ] = moveMaking( A, k, clustering, maxiter, energy_type )
%MOVEMAKIN 此处显示有关此函数的摘要
%   此处显示详细说明

ncut_energies = [];
itrnum = 0;
while 1
    ncut_energies = [ncut_energies NormalizedCutEnergy( A, clustering )];
    disp(['current ' energy_type ' energy: ' num2str(ncut_energies(end))]);
    
    [prob, unaries] = doMove(A, k, clustering, energy_type);
    [~, new_clustering] = min(unaries, [], 2); 
    rw_energies = sum(prob(new_clustering));
    
    %showFigs(img, reshape(new_clustering, [X, Y]));
    % check convergence
    if 0 == sum(new_clustering~=clustering)
        disp('converged!');
        break;
    else
        clustering = new_clustering;
    end
    itrnum = itrnum + 1;
    if itrnum > maxiter
        disp('max iteration');
        break;
    end
end

labels = clustering;



