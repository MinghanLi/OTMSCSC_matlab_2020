function [model,param,DeRain, Rain_Fold, Rains_Fold, h_omega_Fold] = OnlineMSCSC_t(X_Fold, param, model) 
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
    f_size = model.f_size;  tau=zeros(model.type,1);
    B_Omega = model.U*model.v'; 
    Filters = reshape(model.D,[max(f_size),max(f_size),length(f_size)]);
    FInd = FilterInd(f_size); 
    OmegaOut = param.OmegaOut;

    if model.preAlignment == 1
        [TX_Fold,tau,OmegaOut(:,:,end)] = regMGNC(reshape(B_Omega(:,end),[height,wide]),X_Fold,tau,fix(log(height*wide)/log(2)/2-3),2);
    else
        [TX_Fold,tau,~,OmegaOut(:,:,end)] = regImg(reshape(B_Omega(:,end),[height,wide]),X_Fold,tau,ones(height,wide),1);
    end
    TX_pre = model.lastframe; TX = TX_Fold(:);
    TX_Omega = [TX_pre(:), TX];
    T = zeros(size(TX));rho=1; frame = model.frame;
    R = zeros(size(TX_Omega));
        
    % graph cuts initialization
    % GCO toolbox is called
    Omega = param.Omega; 
    ObjArea = sum(~Omega(:)); minObjArea = numel(X_Fold)/1e4; % minimum number of outliers
    Palpha = param.alpha; 
    Pbeta = 0.5*(std(TX))^2;                              % Start from a big value
    converged = false;
    hMRF = GCO_Create(numel(TX_Omega),2);
    GCO_SetSmoothCost( hMRF, [0 1;1 0] );                % Smooth Cost(3DTV)
    AdjMatrix = getAdj(size(TX_Omega),param.weight);
    amplify = 10 * Palpha;
    GCO_SetNeighbors( hMRF, amplify * AdjMatrix );
    energy_cut = 0; energy_old = inf;
  
    iter = 1; maxloopiter = 2;  
    while ~converged && iter <= model.MaxIter
        loopiter = 1;
        while loopiter <= maxloopiter
            E = TX_Omega-B_Omega-R;
            if model.Mask == 0
                Omega = ones(size(TX_Omega)); 
            else     
              %% update Omega paramcters
                % comment these part if there is no moving object
                if Palpha > 0
                    % call GCO to run graph cuts 
                    GCO_SetDataCost( hMRF, int32(5*[ 5*abs(E(:))/std(TX), ~OmegaOut(:) + OmegaOut(:)*max(E(:))/std(TX)]'));
%                     GCO_SetDataCost( hMRF, int32((amplify/alpha)*[ 0.5*(E(:)).^2, ~OmegaOut(:)*Pbeta + OmegaOut(:)*0.5*max(E(:)).^2]'));  
                    GCO_Expansion(hMRF);
                    Omega = ( GCO_GetLabeling(hMRF) == 1 )';
                    energy_cut = energy_cut + double( GCO_ComputeEnergy(hMRF) ); 
                    ObjArea = sum(Omega(:)==0);
                else
                    % direct hard thresholding if no smoothness
                    Omega = 0.5*E.^2 < Pbeta;
                    ObjArea = sum(Omega(:)==0);   
                    energy_cut =  0.5*norm(TX_Omega-B_Omega-E)^2+Pbeta*ObjArea;
                end  
            end
            Omega = reshape(Omega,size(TX_Omega));
            h_omega = Omega(:,end);
 
            %% initial R
            if iter==1; R = h_omega.*E(:,end); end                      
        
            %% Upadata V
            U=model.U; temp=U'.* repmat(h_omega',model.r,1);
            model.v(end,:)=(temp*U+0.0001*eye(model.r))\temp*(TX-R);
            B_Omega = model.U*model.v'; 
            loopiter = loopiter+1;
        end
        h_omega_Fold = reshape(~h_omega, size(X_Fold));
        param.Omega = Omega;
        B = B_Omega(:,end);    
               
        %% Update Foreground F
        Temp = reshape(~h_omega.*(TX-R),size(X_Fold));
        F_Fold = FAD(Temp, model.Flamda, 10, 50, [1e-4, 1e-4],0);
        F = F_Fold(:);

        %% Updata dtau, tau
        temp = reshape(h_omega.*B+~h_omega.*F+R,size(TX_Fold));
        dtau = updataTau(TX_Fold,temp,tau);
        tau = tau + dtau;
        [TX_Fold, OmegaOut(:,:,end)] = warpImg(X_Fold,tau);
        TX = TX_Fold(:);
        
        %% updata U
        model = update_subspace(TX-R,h_omega,model.v(end,:),model);    % updata subspace
        
        %% upadata feathur map
        if param.use_gpu
            [Map, R_conv, Rains]= CSC_ADMM_GPU(reshape(R-T,[height,wide]), Filters, model.b, 100);
        else
            [Map, R_conv]= CSC_ADMM_CPU(reshape(R-T,[height,wide]), Filters, model.b, 100);
        end
        R_conv = R_conv(:);
        
        %% updata filters
        if mod(frame,20)==1
            [Filters,model]=UpFilter_online(reshape(R-T,[height,wide]), model,FInd);
        end
        
        %% updata R
        temp =h_omega.*B+~h_omega.*F+rho*model.Sigma*(R_conv-T);
        R = (TX-temp)/(1+rho*model.Sigma);
        rho = rho*1.05;
        
        %% updata T
        T = T+R_conv-R;
        
        %% updata model param
        Error = h_omega.*(TX-B-R_conv);
        model.Sigma = 1/frame* (sum(Error.^2)/sum(h_omega)+(frame-1)*model.Sigma);
        model.b = 1/frame*(mean(abs(Unfold(Map,size(Map),3)'))+(frame-1)*model.b);
        energy = energy_cut + sum(Error(:).^2);
        if ObjArea > minObjArea && abs(energy_old-energy)/energy < tol; break; end
        energy_old = energy; iter = iter+1;
        DeRain = reshape(h_omega.* B + ~h_omega.*F,[height,wide]);
        Rain_Fold = reshape(R_conv,[height, wide]);  %% for experiments
        Rains_Fold = reshape(Rains,[height, wide, length(model.f_size)]);  %% for experiments
    end
    model.v(1,:) = model.v(2,:);
    model.lastframe = TX_Fold;
    param.OmegaOut(:,:,1) = OmegaOut(:,:,2); 
end

%% function to get the adjacent matirx of the graph
function W = getAdj(sizeData,weight)
numSites = prod(sizeData);
id1 = [1:numSites, 1:numSites, 1:numSites];
id2 = [ 1+1:numSites+1,...
1+sizeData(1):numSites+sizeData(1),...
1+sizeData(1)*sizeData(2):numSites+sizeData(1)*sizeData(2)];
%value = ones(1,3*numSites);
value = [weight*ones(1,2*numSites),2*ones(1,numSites)];
W = sparse(id1,id2,value);
W = W(1:numSites,1:numSites);
end