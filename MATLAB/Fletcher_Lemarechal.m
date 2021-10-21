function [alpha_i] = Fletcher_Lemarechal(f, grad_f, x_k)
    d_k = - grad_f(x_k);
    infini = 1e15;
    alpha_i = 1e-1;
    beta_1 = 1e-3;
    beta_2 = 0.99;
    alpha_r = infini;
    alpha_l = 0;
    gamma = - beta_1 * ( grad_f(x_k) ).' * d_k;
    lambda = 20;
    
    value_find = false;
    i = 0;
    
    %CW1 = @() f(x_k + alpha_i * d_k) <= f(x_k) - alpha_i * gamma;
    %CW2 = @() ( (grad_f(x_k + alpha_i * d_k)).' * d_k ) / ( (grad_f(x_k)).' * d_k ) <= beta_2; 
    
    while(~value_find)
        if( f(x_k + alpha_i * d_k) <= f(x_k) - alpha_i * gamma )
            if( ( (grad_f(x_k + alpha_i * d_k)).' * d_k ) / ( (grad_f(x_k)).' * d_k ) <= beta_2 )
                value_find = true;
            else
                alpha_l = alpha_i;
                if (alpha_r < infini)
                    alpha_i = ( alpha_l + alpha_r )/2;
                else
                    alpha_i = lambda * alpha_i;
                end
                
            end
        else
            alpha_r = alpha_i;
            alpha_i = ( alpha_l + alpha_r )/2;
        end
        i = i + 1;
    end
    
end