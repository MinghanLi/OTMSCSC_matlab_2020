 function [model,DeRain, Rain_Fold, Rains_Fold, h_omega_Fold] = OnlineMSCSC(X_Fold, param, model)  
    if (~isfield(param,'tol'))
        tol = 1.0e-7;
    else
        tol = param.tol;
    end
    if (~isfield(param,'lambda'))
        param.lambda = 5;
    end     
    if (~isfield(param,'weight'))
        param.weight = 1;
    end
    if (~isfield(model,'ro'))
       model.ro=0.99;
    end
    if (~isfield(param,'sigma'))
        param.sigma = [];
    end 
    
    [height,wide] = size(X_Fold);
    f_size = model.f_size;
    B_Omega = model.U*model.v';
    Filters = reshape(model.D,[max(f_size),max(f_size),length(f_size)]);
    FInd = FilterInd(f_size); 

    X_pre = model.lastframe; X = X_Fold(:);
    X_Omega = [X_pre(:), X];
    T = zeros(size(X));rho=1;
        
    % graph cuts initialization
    % GCO toolbox is called
    Omega = param.Mask; 
    ObjArea = sum(~Omega(:)); minObjArea = numel(X)/1e4; % minimum number of outliers
    Palpha = param.alpha; 
    converged = false;
    hMRF = GCO_Create(numel(X_Omega),2);
%     GCO_SetLabeling(hMRF,Omega(:)')
    GCO_SetSmoothCost( hMRF, [0 1;1 0] );                % Smooth Cost(3DTV)
    AdjMatrix = getAdj([size(X_Fold),size(X_Omega,2)],param.weight);
    GCO_SetNeighbors( hMRF, 10*Palpha * AdjMatrix );
    energy_cut = 0; energy_old = inf;
     
     iter = 1; maxloopiter = 2;  
     while ~converged && iter <= model.MaxIter   
        %% update Omega paramcters
        loopiter = 1;
        while loopiter <= maxloopiter 
            E = B_Omega-X_Omega;
            if model.Mask==0
               Omega = ones(size(X_Omega));
            else
                % comment these part if there is no moving object
                if Palpha > 0
                    % call GCO to run graph cuts  
                    GCO_SetDataCost( hMRF, int32(10*[ 5*abs(E(:))/std(X), ones(numel(X_Omega),1)]'));
                    GCO_Expansion(hMRF);
                    Omega = ( GCO_GetLabeling(hMRF) == 1 )';
                    energy_cut = energy_cut + double( GCO_ComputeEnergy(hMRF) ); 
                    ObjArea = sum(Omega(:)==0);
                else
                    % direct hard thresholding if no smoothness
                    Omega = E.^2 > (std(E(:)))^2;
                    ObjArea = sum(Omega(:)==0);   
                    energy_cut =  0.5*norm(X_Omega-B_Omega-E)^2+(std(E(:)))^2*ObjArea;
                end  
                Omega = reshape(Omega,size(X_Omega));
            end 
            h_omega = Omega(:,end);
            h_omega_Fold = reshape(~h_omega, size(X_Fold));
 
          %% initial Q
            if iter == 1; Q = h_omega.*E(:,end); end            
          %% Upadata V
            model.v(end,:)=(model.U'.* repmat(h_omega',model.r,1)*model.U+0.0001*eye(model.r))^-1*model.U'*(h_omega.*(X-Q));
            B_Omega = model.U*model.v';
            loopiter = loopiter+1;
        end
        param.Omega = Omega;
        B=B_Omega(:,end);
%         B_Fold=reshape(B,size(X_Fold));
            
       %% updata U
        model=update_subspace(X-Q,h_omega,model.v(end,:),model);% updata subspace  
            
       %% Update Foreground F
        Temp = reshape(~h_omega.*(X-Q),size(X_Fold));
        % F_Fold = FAD(Temp, model.Flam, 10, 50, [1e-4, 1e-4],0);
        F_Fold = TVL1denoise(Temp, model.Flam, 50);
        F = F_Fold(:);
%         F = ~h_omega.*(X-Q);
              
       %% upadata feathur map 
        if param.use_gpu
            [Map, Rain, Rains]= CSC_ADMM_GPU(reshape(Q-T,[height,wide]), Filters, model.b, 100);
        else
            [Map, Rain]= CSC_ADMM_CPU(reshape(Q-T,[height,wide]), Filters, model.b, 100);
        end
        Rain = Rain(:);
         
       %% updata filters
        if mod(model.frame,10) == 1
            [Filters,model] = UpFilter_online(reshape(Q-T,[height,wide]), model, FInd);
        end
           
            
       %% updata Q
        Temp =h_omega.*B+~h_omega.*F+rho*model.Sigma*(Q+T);
        Q = (X-Temp)/(1-rho*model.Sigma);
        rho = rho*1.05;
            
       %% updata T
        T = T+Rain-Q;
    
       %% updata model param
        frame = model.frame;
        Error = h_omega.*(X-B-Rain);
        model.Sigma = 1/frame* (sum(Error.^2)/sum(h_omega)+(frame-1)*model.Sigma);
        model.b = 1/frame*(mean(abs(Unfold(Map,size(Map),3)'))+(frame-1)*model.b);
        energy = energy_cut + sum(Error(:).^2);
        if ObjArea > minObjArea && abs(energy_old-energy)/energy < tol
            break;
        end
        energy_old = energy; iter = iter+1;
        DeRain = reshape(h_omega.* B + ~h_omega.*F,[height,wide]);
        Rain_Fold = reshape(Rain,[height, wide]);  %% for experiments
        Rains_Fold = reshape(Rains,[height, wide, length(model.f_size)]);  %% for experiments
    end
    model.v(1,:) = model.v(2,:);
    model.lastframe = X_Fold;
end


%% function to get the adjacent matirx of the graph
function W = getAdj(sizeData)
numSites = prod(sizeData);
id1 = [1:numSites, 1:numSites, 1:numSites];
id2 = [ 1+1:numSites+1,...
1+sizeData(1):numSites+sizeData(1),...
1+sizeData(1)*sizeData(2):numSites+sizeData(1)*sizeData(2)];
%value = ones(1,3*numSites);
value = [ones(1,2*numSites),ones(1,numSites)];
W = sparse(id1,id2,value);
W = W(1:numSites,1:numSites);
end
    
