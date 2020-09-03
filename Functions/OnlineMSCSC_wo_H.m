 function [model,DeRain, Rain_Fold, Rains_Fold] = OnlineMSCSC_wo_H(X_Fold, param, model)  
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
    B = model.U*model.v(end)';
    Filters = reshape(model.D,[max(f_size),max(f_size),length(f_size)]);

    X_pre = model.lastframe; 
    X = X_Fold(:);
    T = zeros(size(X));rho=1;
    R = zeros(size(X));
        
    % graph cuts initialization
    % GCO toolbox is called
    converged = false;
     
     iter = 1; 
     while ~converged && iter <= model.MaxIter   
        E = X - B;
        sign = E < 0;
        
        %% initial R
        if iter == 1; R = E; end
        %% updata filters
        % if mod(model.frame,500) == 1
        %    [Filters,model] = UpFilter_online(reshape(R-T,[height,wide]), model, FInd);
        % end
        
        %% upadata feathur map 
        if model.use_gpu
            [Map, R_conv, Rains]= CSC_ADMM_GPU(reshape(max(R-T, 0),[height,wide]), Filters, model.b, 100);
        else
            [Map, R_conv]= CSC_ADMM_CPU(reshape(R-T,[height,wide]), Filters, model.b, 100);
        end
        R_conv = R_conv(:);
        
        %% updata R
        Temp = B-rho*model.Sigma*(R_conv+T);
        R = (X-Temp)/(1+rho*model.Sigma);
        rho = rho*1.05;
            
       %% updata T
        T = T+R_conv-R;
    
       %% updata model param
       iter = iter + 1;
        frame = model.frame;
        Error = X-B-R_conv;
        model.Sigma = 1/frame* (sum(Error.^2)/size(X, 1)+(frame-1)*model.Sigma);
        model.b = 1/frame*(mean(abs(Unfold(Map,size(Map),3)'))+(frame-1)*model.b);
        DeRain = reshape(X-R_conv, [height,wide]);
        Rain_Fold = reshape(R_conv,[height, wide]);  %% for experiments
        Rains_Fold = reshape(Rains,[height, wide, length(model.f_size)]);  %% for experiments
     end
    model.lastframe = X_Fold;
 end
