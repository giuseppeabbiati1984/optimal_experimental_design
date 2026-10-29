function Q = fun_obj_VAR(theta,M,Km,Dm,w,X)

Ntest = length(Km) ;
%H_hat = cell(3,1) ;

F2_hat = zeros(1,length(X{1}),Ntest) ;

for i = 1:1:Ntest

    H_hat = (-w.^2*(theta(1) + M) + 1i*w*(theta(2) + Dm(i)) + Km(i)) ;

    F2_hat(1,:,i) = H_hat.*X{i} ;

end

idx = find(F2_hat(1,:,1)==0); 
F2_hat(:,idx,:) = [] ;

Q = sum(sum(var(F2_hat,0,3),2)) ;%*(Ntest-1) ; 