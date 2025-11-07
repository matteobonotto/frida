function xx = fun_DAE_Crank_Nicolson_State_Space(time,xx_t0,uu,E_SS,A_SS,B_SS)


n_dim = size(A_SS,1);
n_time = length(time);

h_step = (time(2:end)-time(1:end-1));
if any(abs(h_step - h_step(1)) > 1e-8) % constant time steps
    
    time_new = linspace(time(1),time(end),n_time);
    
    uu = interp1(time,uu',time_new)';
    time = time_new;
    
end

h_step = h_step(1);

xx = zeros(n_dim,n_time);

xx(:,1) = xx_t0;



%% solution

EE = eye(size(E_SS));

M_1 = (E_SS - .5*h_step*A_SS);
inv_M_1 = M_1\EE;

M_2 = (E_SS + .5*h_step*A_SS);

A_CN = inv_M_1*M_2;
B_CN = .5*h_step*inv_M_1*B_SS;


xx_t = xx_t0;
uu_t = uu(:,1);

for ii=2:n_time
    
    uu_t1 = uu(:,ii);
    
    xx_t1 = A_CN*xx_t + B_CN*(uu_t + uu_t1);

    xx_t = xx_t1;
    uu_t = uu_t1;
    
    xx(:,ii) = xx_t;

end




end











