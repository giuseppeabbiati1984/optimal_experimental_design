function Q = fun_obj_CPSD(theta,M,Km,Dm,Gmat,w,options)

% Check if 'options' is provided
% ------------------------------------------------------------------------
if ~exist('options', 'var') || isempty(options)
    options = struct(); 
end

if ~isfield(options, 'wcut')
    options.wcut = w(end); 
end

if ~isfield(options, 'idx')
    options.idx = [1,2]; 
end

if ~isfield(options, 'diff')
    options.diff = 'NLSE'; 
end

list = abs(w(1,:))<options.wcut ;
n = length(Gmat) ;

% Some index manipulation
% ------------------------------------------------------------------------
p = options.idx(1) ;
q = options.idx(2) ;

% make sure that we are only considering cross spectra 
if p == q 
    
    idx = setdiff(1:1:n,p) ;

    q = idx(randi(length(idx))) ;

    disp(['One index has been changed to ',num2str(q)])

end

% To be consistent we are only considering one part of matrix
if p > q

    [p,q] = deal(q,p) ;

end

% creating transfer functions 
% ------------------------------------------------------------------------
H_hat = cell(n,1) ;

for i = 1:1:n

    H_hat{i} =...
        (-w.^2*(theta(1) + M) + 1i*w*(theta(2) + Dm(i)) + Km(i)).^-1 ;

end

% create cross spectral estimates
% ------------------------------------------------------------------------
%fac = cell(N,N) ;

GN_hat = cell(n-1,1) ;
Z = cell(n-1,1) ;

for i = 1:1:n-1

        fac = conj(H_hat{i}).*H_hat{n}./(conj(H_hat{p}).*H_hat{q}) ;


        %Zpq{i,j} = ...
        % conj(H_hat{p}).*H_hat{q}./(conj(H_hat{i}).*H_hat{j}) ;
                Z{i} = fac ;
                GN_hat{i} = Gmat{i,n}.*list ;

end

% Compute residuals                                                          
% ------------------------------------------------------------------------

switch options.diff                                                 

    case 'NLSE'
    % --------------------------------------------------------------------
        
        idx = sum(list) ;
        
        %Z{i} = [] ;
        Z_reshape = cell2mat(Z(~cellfun('isempty', Z))) ;
        %Zpq_reshape(:,end) = [] ;

        %Gpq_hat{p,q} = [] ;
        GN_reshape = cell2mat(GN_hat(~cellfun('isempty', GN_hat))) ;
        GN_reshape(:,end) = [] ;
        Gmat{p,q}(end) = [] ;
        
        Q = 0 ;

         for i = 1:1:n-1       
                
            E = sum(abs(Z_reshape(i,1:idx) ...
                - (GN_reshape(i,1:idx)./Gmat{p,q}(1,1:idx))).^2) ; 
            
            Q = Q + E ;
            
         end
 
end

end
