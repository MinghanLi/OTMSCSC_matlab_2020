function [Map,RainS]=UpdataMap(X,Filters,Q,H,model)
 T = zeros(size(X)); rho=1;
 iter=1;MaxIter=2;
     while iter<=MaxIter
         %% upadata feathur map 
         [Map, RainS ] = SynthesisSC(Q+T, Filters, model, rho ); % Z[(h+f_size-1),(w+f_size-1),K]; RainS [h,w,K]
         Rain = reshape(sum(RainS,3), size(X));
         
         %% updata Q
         Q = Rain-T;
         QM = (rho*model.Sigma*Q+X)./(1+rho*model.Sigma);
         Q(H==1) = QM(H==1);
         rho = rho*1.05;
         
         %% updata T
         T = T+Rain-Q;
         iter=iter+1;
     end
end


