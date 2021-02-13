 function DeRain = OnlineMSCSC_SS(X_Fold, param, model)  
    if (~isfield(param,'tol'))
        tol = 1.0e-7;
    else
        tol = param.tol;
    end
    if (~isfield(param,'lambda'))
        param.lambda = 5;
    end     
    if (~isfield(param,'sigma'))
        param.sigma = [];
    end       
    if (~isfield(param,'weight'))
        param.weight = 1;
    end
    if (~isfield(model,'ro'))
       model.ro=0.99;
    end
    
    [height,wide] = size(X_Fold);
    B_Omega = model.U*model.v'; 

    X_pre = model.lastfrmae;
    X = X_Fold(:);
    X_Omega = [X_pre(:), X];
    T = zeros(size(X));rho=1;
        
    % graph cuts initialization
    % GCO toolbox is called
    Omega = param.Omega; 
    ObjArea = sum(~Omega(:)); minObjArea = numel(X)/1e4; % minimum number of outliers
    Psigma = param.sigma; Plambda = param.lambda; Prho = param.rho;
    Pbeta = 0.5*(std(X))^2;                              % Start from a big value
    minbeta = 0.5*(3*std(X)/20)^2;                       % lower bound: suppose SNR <= 20
    converged = false;
    hMRF = GCO_Create(numel(X_Omega),2);
    GCO_SetSmoothCost( hMRF, [0 1;1 0] );                % Smooth Cost(3DTV)
    AdjMatrix = getAdj(size(X_Omega),param.weight);
    amplify = 10 * Plambda;
    GCO_SetNeighbors( hMRF, amplify * AdjMatrix );
    energy_cut = 0; energy_old = inf;
     
     iter = 1; maxloopiter = 2;  
     while ~converged && iter <= model.MaxIter   
        %% update Omega paramcters
        loopiter = 1;
        while loopiter <= maxloopiter
            E = X_Omega-B_Omega;
            if model.Mask==0
               Omega = ones(size(X_Omega));
            else
               if isempty(Psigma)
                    sigma_old = Psigma;
                    residue = sort(E(Omega(:)));
                    truncate = 0.005;
                    idx1 = round(truncate*length(residue))+1;  
                    idx2 = round((1-truncate)*length(residue));
                    Psigma = std(residue(idx1:idx2));
                    if abs(sigma_old-Psigma)/abs(sigma_old) < 0.01
                        Psigma = sigma_old;
                    end
                end
    
                if ObjArea < minObjArea
                    Pbeta = Pbeta/2;
                else
                    Pbeta =min(max([Pbeta/2,0.5*(3*Prho*Psigma)^2,minbeta]),Pbeta);
                end
                alpha = Plambda * Pbeta; 
        
                % comment these part if there is no moving object
                if Plambda > 0
                    % call GCO to run graph cuts  
                    GCO_SetDataCost( hMRF, int32((amplify/alpha)*[ 0.5*(E(:)).^2, ones(numel(X_Omega),1)*Pbeta ]'));  
                    GCO_Expansion(hMRF);
                    Omega = ( GCO_GetLabeling(hMRF) == 1 )';
                    energy_cut = energy_cut + double( GCO_ComputeEnergy(hMRF) ); 
                    ObjArea = sum(Omega(:)==0);
                    energy_cut = (alpha/amplify) * energy_cut;
                else
                    % direct hard thresholding if no smoothness
                    Omega = 0.5*E.^2 < Pbeta;
                    ObjArea = sum(Omega(:)==0);   
                    energy_cut =  0.5*norm(X_Omega-B_Omega-E)^2+Pbeta*ObjArea;
                end  
                Omega = reshape(Omega,size(X_Omega));
            end 
            h_omega = Omega(:,end);
 
          %% initial Q
            if iter == 1 
                Q = h_omega.*E(:,end);
            end            
          %% Upadata V
            model.v(end,:)=(model.U'.* repmat(h_omega',model.r,1)*model.U+0.0001*eye(model.r))^-1*model.U'*(h_omega.*(X-Q));
            B_Omega = model.U*model.v';
            loopiter = loopiter+1;
        end
        param.Omega = Omega;
        B=B_Omega(:,end);
            
       %% updata U
        model=update_subspace(X-Q,h_omega,model.v(end,:),model);% updata subspace          
            
       %% Update Foreground F
        %Temp = reshape(~h_omega.*(X-Q),size(X_Fold));
        % F_Fold = FAD(Temp, model.Flam, 10, 50, [1e-4, 1e-4],0);
%         F_Fold = TVL1denoise(Temp, model.Flam, 50);
%         F = F_Fold(:);
        F = ~h_omega.*(X-Q);
              
       %% updata model param
        iter = iter+1;
        DeRain = reshape(h_omega.* B + ~h_omega.*F,[height,wide]);
    end
    model.v(1,:) = model.v(2,:);
end


%% function to get the adjacent matirx of the graph
function W = getAdj(sizeData,weight)
numSites = prod(sizeData);
id1 = [1:numSites, 1:numSites, 1:numSites];
id2 = [ 1+1:numSites+1,...
1+sizeData(1):numSites+sizeData(1),...
1+sizeData(1)*sizeData(2):numSites+sizeData(1)*sizeData(2)];
%value = ones(1,3*numSites);
value = [weight*ones(1,2*numSites),1*ones(1,numSites)];
W = sparse(id1,id2,value);
W = W(1:numSites,1:numSites);
end
    
function [ M, Rain ] = SynthesisSC( X, Filters, par, FInd )
% UNTITLED2 Summary of this function goes here
% Detailed explanation goes here
[h,w] = size(X); mu=par.b; f_size=par.f_size;
for k=1:size(Filters,3)
    Fk = Filters(:,:,k);Temp = reshape(Fk(FInd(:,:,k)),[f_size(k),f_size(k)]);
    FFilters(:,:,k) = psf2otf(rot90(Temp,2),[h,w]);       % [h+2*CutEdge,w+2*CutEdge,K]
end
[M, Rain]= CSC_ADMM_CPU( fft2(X), single(FFilters), mu, 200);
end

