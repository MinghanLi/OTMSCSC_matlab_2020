 function [model,DeRain, Rain_Fold, Rains_Fold, h_omega_Fold] = OnlineMSCSC_new(X_Fold, param, model)  
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

    X_pre = model.lastframe; 
    X = X_Fold(:);
    X_Omega = [X_pre(:), X];
    T = zeros(size(X));rho=1;
    R = zeros(size(X_Omega));
        
    % graph cuts initialization
    % GCO toolbox is called
    Omega = param.Mask; 
    ObjArea = sum(~Omega(:)); minObjArea = numel(X)/1e4; % minimum number of outliers
    Palpha = param.alpha; 
    Psigma = param.sigma; 
    Pbeta = 0.5*(std(X))^2;                              % Start from a big value
    minbeta = 0.5*(3*std(X)/20)^2; 
    converged = false;
    hMRF = GCO_Create(numel(X_Omega),2);
    GCO_SetLabeling(hMRF,Omega(:)')                    % initial
    GCO_SetSmoothCost( hMRF, [0 1;1 0] );                % Smooth Cost
    AdjMatrix = getAdj([size(X_Fold),size(X_Omega,2)],param.weight);
    % AdjMatrix = getAdj_grad(reshape(X_Omega, [size(X_Fold),size(X_Omega,2)]));
    GCO_SetNeighbors( hMRF, 10*Palpha * AdjMatrix );
    energy_cut = 0; energy_old = inf;
     
     iter = 1; LoopIter=1; LoopMaxIter=2;
     while ~converged && iter <= model.MaxIter   
        %% update Omega paramcters
        while LoopIter <= LoopMaxIter
            E = X_Omega-B_Omega;
            if model.isMask==0
                Omega = ones(size(X_Omega));
            else
                % comment these part if there is no moving object
                if isempty(Psigma)
                    sigma_old = Psigma;
                    residue = sort(abs(E(Omega(:))));
                    truncate = 0.05;
                    idx1 = round(truncate*length(residue))+1;  
                    idx2 = round((1-truncate)*length(residue));
                    Psigma = std(residue(idx1:idx2));
                    if abs(sigma_old-Psigma)/abs(sigma_old) < 0.01
                        Psigma = sigma_old;
                    end
                end
    
                if ObjArea <  minObjArea
                    Pbeta = Pbeta/2;
                else
                    Pbeta =min(max([Pbeta/2,0.5*(3*Psigma)^2,minbeta]),Pbeta);
                end
                
                if Palpha > 0
                    sigma = 1;
                    window = double(uint8(3*sigma)*2+1);
                    H_gauss = fspecial('gaussian', window, sigma);
                    E_gauss = imfilter(reshape(E, [360,540, 2]), H_gauss, 'replicate');
                    FG_cost = E_gauss;
                    % call GCO to run graph cuts
                    % GCO_SetDataCost( hMRF, int32(5/Pbeta*[ 5*abs(E(:))./X_Omega(:), ones(numel(X_Omega),1) ]'));
                    GCO_SetDataCost( hMRF, int32(10*[ 0.5*(FG_cost(:)).^2/Pbeta, ones(numel(X_Omega),1) ]'));  
                    GCO_Expansion(hMRF);
                    Omega = ( GCO_GetLabeling(hMRF) == 1 )';
                    energy_cut = energy_cut + double( GCO_ComputeEnergy(hMRF) );
                    ObjArea = sum(Omega(:)==0);
                    energy_cut = Pbeta/10 * energy_cut;
                else
                    % direct hard thresholding if no smoothness
                    Omega = E.^2 > Pbeta; % (std(E(:)))^2;
                    ObjArea = sum(Omega(:)==0);
                    energy_cut =  0.5*norm(X_Omega-B_Omega-E)^2+Pbeta*ObjArea;
                end
                Omega = reshape(Omega,size(X_Omega));
            end
            h_omega = Omega(:,end);
            
            %% initial R
            if iter == 1
                R = h_omega.*E(:,end); 
            end
            
            %% Upadata V
            model.v(end,:)=(model.U'.* repmat(h_omega',model.r,1)*model.U+0.0001*eye(model.r))^-1*model.U'*(h_omega.*(X-R));
            B_Omega = model.U*model.v';
            LoopIter = LoopIter+1;
        end
        B=B_Omega(:,end);
        param.Omega = Omega;
        h_omega_Fold = reshape(~h_omega, size(X_Fold));
        
        %% Update Foreground F
        Temp = reshape(~h_omega.*(X-R),size(X_Fold));
        F_Fold = FAD(Temp, model.Flamda, 10, 50, [1e-4, 1e-4],0);
%         F_Fold = TVL1denoise(Temp, model.Flam, 50);
            
        %% updata U
        model=update_subspace(X-R,h_omega,model.v(end,:),model);% updata subspace  
                 
        %% updata filters
%         if mod(model.frame,10) == 1
%             [Filters,model] = UpFilter_online(reshape(R-T,[height,wide]), model, FInd);
%         end

        %% upadata feathur map 
        sign = (R-T) < 0;
        if model.use_gpu
            [Map, R_conv, Rains]= CSC_ADMM_GPU(reshape(abs(R-T),[height,wide]), Filters, model.b, 100);
        else
            [Map, R_conv]= CSC_ADMM_CPU(reshape(R-T,[height,wide]), Filters, model.b, 100);
        end
        R_conv(sign == 1) = -1 * R_conv(sign == 1);
        R_conv = R_conv(:);
        
         %% updata R
        Temp =h_omega.*B+~h_omega.*F_Fold(:)-rho*model.Sigma*(R_conv+T);
        R = (X-Temp)/(1+rho*model.Sigma);
        rho = rho*1.05;
    
       %% updata T
        T = T+R_conv-R;
    
       %% updata model param
        frame = model.frame;
        Error = h_omega.*(X-B-R_conv);
        model.Sigma = 1/frame* (sum(Error.^2)/sum(h_omega)+(frame-1)*model.Sigma);
        model.b = 1/frame*(mean(abs(Unfold(Map,size(Map),3)'))+(frame-1)*model.b);
        energy = energy_cut + sum(Error(:).^2);
        if ObjArea > minObjArea && abs(energy_old-energy)/energy < tol
            break;
        end
        energy_old = energy; iter = iter+1;
        % imshow(reshape(F_Fold, [height, wide]))
        if model.derain == 1
            DeRain = reshape(h_omega.* B + ~h_omega.*F_Fold(:),[height,wide]);
        elseif model.derain == 2
            DeRain = reshape(X-R_conv, [height,wide]);
        end
        Rain_Fold = reshape(R_conv,[height, wide]);  %% for experiments
        Rains_Fold = reshape(Rains,[height, wide, length(model.f_size)]);  %% for experiments
    end
    model.v(1,:) = model.v(2,:);
    model.lastframe = X_Fold;
end


%% function to get the adjacent matirx of the graph
function W = getAdj(sizeData,weight)
numSites = prod(sizeData);
id1 = [1:numSites, 1:numSites, 1:numSites];
id2 = [ 
    1+1:numSites+1,...
    1+sizeData(1):numSites+sizeData(1),...
    1+sizeData(1)*sizeData(2):numSites+sizeData(1)*sizeData(2)];
value = weight * ones(1,3*numSites);
W = sparse(id1,id2,value);
W = W(1:numSites,1:numSites);
end

%% function to get the adjacent matirx of the graph
function W = getAdj_grad(X)
sizeData = size(X);
numSites = prod(sizeData);
[Fx,Fy, Fz] = gradient(X);
Fx = abs(Fx)*100;
Fy = abs(Fy)*100;
Fz = abs(Fz)*100;

id1 = [1:numSites, 1:numSites, 1:numSites];
id2 = [ 
    1+1:numSites+1,...
    1+sizeData(1):numSites+sizeData(1),...
    1+sizeData(1)*sizeData(2):numSites+sizeData(1)*sizeData(2)];
value = [Fx, Fy,Fz];
W = sparse(id1,id2,round(value(:)));
W = W(1:numSites,1:numSites);
end
    
