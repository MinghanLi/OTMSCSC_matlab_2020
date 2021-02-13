function [model,param,DeRain1, DeRain2, Rain_Fold, Rains_Fold, Curr_Mask_Fold] = OnlineMSCSC_trans(X_Fold, param, model) 
    if (~isfield(param,'tol'))
        tol = 1.0e-7;
    else
        tol = param.tol;
    end
    if (~isfield(param,'alpha'))
        param.lambda = 10;
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
    B_Fold = model.B_Fold; 
    f_size = model.f_size; 
    Filters = reshape(model.D,[max(f_size),max(f_size),length(f_size)]);
    FInd = FilterInd(f_size); 
    tau=zeros(model.type,1);
    
    % align background for current frame
    if model.preAlignment == 1
        [TB_Fold,tau,MaskOut1] = regMGNC(X_Fold,B_Fold,tau,fix(log(height*wide)/log(2)/2-3),2);
    else
        [TB_Fold,tau,~,MaskOut1] = regImg(X_Fold,B_Fold,tau,ones(height,wide),1);
    end
    TB_Fold(MaskOut1)=X_Fold(MaskOut1);
    TB_Mask=repmat(TB_Fold(:),[1,2]);
    MaskOut_Mask=repmat(MaskOut1,[1,1,2]);

    X_pre = model.lastframe;            
    [TX_pre,~] = preAlign(X_pre,X_Fold,zeros(model.type,1));
    X = X_Fold(:); X_Mask = [TX_pre(:), X];
    T = zeros(size(X)); R = zeros(size(X));
    rho=1; frame = model.frame;
        
    % graph cuts initialization
    % GCO toolbox is called
    Mask = param.Mask; 
    ObjArea = sum(~Mask(:)); minObjArea = numel(X_Fold)/1e4; % minimum number of outliers
    Psigma = param.sigma; 
    Palpha = param.alpha; 
    Pbeta = 0.5*(std(X))^2;                              % Start from a big value
    minbeta = 0.5*(3*std(X)/20)^2;
    converged = false;
    hMRF = GCO_Create(numel(X_Mask),2);
    GCO_SetSmoothCost( hMRF, [0 1;1 0] );                % Smooth Cost(3DTV)
    AdjMatrix = getAdj(reshape(X_Mask,[size(X_Fold),2]));
    GCO_SetNeighbors( hMRF, Palpha*AdjMatrix );
    energy_cut = 0; energy_old = inf;
  
    iter = 1; 
    while ~converged && iter <= model.MaxIter
        E_Mask = X_Mask-TB_Mask;
        if model.isMask == 0
            Mask = ones(size(X_Mask));
        else
            if isempty(Psigma)
                sigma_old = Psigma;
                residue = sort(E_Mask(Mask(:)));
                truncate = 0.005;
                idx1 = round(truncate*length(residue))+1;
                idx2 = round((1-truncate)*length(residue));
                Psigma = std(residue(idx1:idx2));
                if abs(sigma_old-Psigma)/abs(sigma_old) < 0.01
                    Psigma = sigma_old;
                end
            end
            
            % update beta
            if ObjArea < minObjArea
                Pbeta = Pbeta/2;
            else
                Pbeta =min(max([Pbeta/2,3*(param.rho*Psigma)^2,minbeta]),Pbeta);
            end
            % comment these part if there is no moving object
            if Palpha > 0
                % call GCO to run graph cuts
                GCO_SetDataCost( hMRF, int32((10/Pbeta)*[ 0.5*(E_Mask(:)).^2, ~MaskOut_Mask(:)*Pbeta + MaskOut_Mask(:)*0.5*max(E_Mask(:)).^2 ]'));
                % GCO_SetDataCost( hMRF, int32(5/std(X)*[ 5*abs(max(E_Mask(:),0)), ~MaskOut_Mask(:)*std(X) + MaskOut_Mask(:)*max(E_Mask(:))]'));
%                 GCO_SetDataCost( hMRF, int32(5/(std(X))^2*[ 2*(E(:)).^2, ~OmegaOut(:)*(std(X))^2 + OmegaOut(:)*max(E(:)).^2]'));
                GCO_Expansion(hMRF);
                Mask = ( GCO_GetLabeling(hMRF) == 1 )';
                energy_cut = energy_cut + double( GCO_ComputeEnergy(hMRF) );
                ObjArea = sum(Mask(:)==0);
            else
                % direct hard thresholding if no smoothness
                Mask = 0.7*E_Mask.^2 < Pbeta;
                ObjArea = sum(Mask(:)==0);
                energy_cut =  0.5*norm(X_Mask-TB_Mask-E_Mask)^2+Pbeta*ObjArea;
            end
        end
        Mask = reshape(Mask,size(X_Mask));
        Curr_Mask = Mask(:,end);
        Curr_Mask_Fold = reshape(~Curr_Mask, size(X_Fold));
        param.Mask = Mask;
        
        %% initial R
        if iter==1
            % R = max(X_Mask(:,2)-X_Mask(:,1), 0);
            R = Curr_Mask.*reshape(E_Mask(:,2),[height*wide,1]);
        end 
        
        %% Update Foreground F
        Temp = reshape(~Curr_Mask.*(X-R),size(X_Fold));
        F_Fold = FAD(Temp, model.Flamda, 10, 50, [1e-4, 1e-4],0);
        F = F_Fold(:);

        %% Updata dtau, tau
        temp = reshape(Curr_Mask.*(X-R),size(X_Fold));
        dtau = updataTau(TB_Fold,temp,tau);
        tau = tau + dtau;
        
        %% Updata Background
        [TB_Fold, MaskOut1] = warpImg(B_Fold,tau);
        MaskOut=repmat(MaskOut1,[1,1,2]);
        TB_Fold(MaskOut1)=X_Fold(MaskOut1);
        TB = TB_Fold(:);
        
        %% upadata feature map
        if model.use_gpu
            [Map, R_conv, Rains]= CSC_ADMM_GPU(reshape(R-T,[height,wide]), Filters, model.b,100);
        else
            [Map, R_conv]= CSC_ADMM_CPU(reshape(R-T,[height,wide]), Filters, model.b,100);
        end
        % find horizontal filters 
        [max_filter, max_filter,K] = size(Filters);
        n_col = sum(sum(Filters > 0, 1) > 0, 2);    % column
        n_row = sum(sum(Filters > 0, 2) > 0, 1);   % row
        hori_filters = reshape(n_col > n_row, [1, K]);
        for k = 1:K
            if hori_filters(k)== 1
                R_conv = R_conv-Rains(:,:,1,k);
            end
        end
        R_conv = R_conv(:);
        
        %% updata R
        temp =Curr_Mask.*TB+~Curr_Mask.*F-rho*model.Sigma*(R_conv+T);
        R = (X-temp)/(1+rho*model.Sigma);
        rho = rho*1.05;
        
        %% updata filters
        if mod(frame,100)==1
            [Filters,model]=UpFilter_online(reshape(R-T,[height,wide]), model, FInd);
        end

        
        %% updata T
        T = T+R-R_conv;

        %% updata model param
        Error = Curr_Mask.*(X-TB-R_conv);
        model.Sigma = 1/frame* (sum(Error.^2)/sum(Curr_Mask)+(frame-1)*model.Sigma);
        model.b = 1/frame*(mean(abs(Unfold(Map,size(Map),3)'))+(frame-1)*model.b);
        energy = energy_cut + sum(Error(:).^2);
        if ObjArea > minObjArea && abs(energy_old-energy)/energy < tol; break; end
        energy_old = energy; iter = iter+1;
    end
    imshow([X_Fold, TB_Fold, F_Fold, Curr_Mask_Fold])
    % if model.derain == 1
    % DeRain1 = reshape(Curr_Mask.* TB + ~Curr_Mask.*F,[height,wide]);
    DeRain1 =reshape( Curr_Mask.* TB + ~Curr_Mask.*(X) ,[height,wide]);
    % imshow([DeRain1, DeRain_a])
    % elseif model.derain == 2
    R_mask = (R_conv > 0).*Curr_Mask;
    DeRain2 = reshape(~R_mask.*X+R_mask.*TB, [height,wide]);
    % end
    Rain_Fold = reshape(R_conv,[height, wide]);  %% for experiments
    Rains_Fold = reshape(Rains,[height, wide, length(model.f_size)]);  %% for experiments
    model.lastframe = X_Fold;
    model.B_Fold_Omega = TB_Fold;
end

%% function to get the adjacent matirx of the graph
function W = getAdj(img2)
sizeData = size(img2);
numSites = prod(sizeData);
id1 = [1:numSites-1, 1:numSites-sizeData(1), 1:numSites/sizeData(3)];
id2 = [ 1+1:numSites,...
1+sizeData(1):numSites,...
1+sizeData(1)*sizeData(2):numSites];
%value = ones(1,3*numSites);
var_img = var(img2(:)); 
value = ceil(10*exp(-(img2(id1)-img2(id2)).^2/(2*var_img)));
% value = [1*ones(1,2*numSites),weight*ones(1,numSites)];
W = sparse(id1,id2,value,numSites,numSites);
end