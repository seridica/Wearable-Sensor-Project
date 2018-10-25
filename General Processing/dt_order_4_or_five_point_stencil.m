function dxdt = dt_order_4_or_five_point_stencil(t,x)

    if size(x,2) ~= size(t,2)
        error('Vectors must be of the same length')
    elseif (size(x,2) < 5)
        error('Not enough data points for derivative of accuracy order 4')
    else
        
        % Assume constant time step
        h = t(:,2) - t(:,1);
        
        for i=1:size(t,2)
            
            if i==1 % Forward difference, 4th order
                dxdt(:,i) = (-25/12*x(:,i) + 4*x(:,i+1) - 3*x(:,i+2) + 4/3*x(:,i+3) - 1/4*x(:,i+4))/h;
            elseif i==2 % Forward difference, 3rd order
                dxdt(:,i) = (-11/6*x(:,i) + 3*x(:,i+1) - 3/2*x(:,i+2) + 1/3*x(:,i+3))/h;
            elseif i==length(t)-1 % Backward difference
                dxdt(:,i) = (11/6*x(:,i) - 3*x(:,i-1) + 3/2*x(:,i-2) - 1/3*x(:,i-3))/h;
            elseif i==length(t) % Backward difference
                dxdt(:,i) = (25/12*x(:,i) - 4*x(:,i-1) + 3*x(:,i-2) - 4/3*x(:,i-3) + 1/4*x(:,i-4))/h;
            else % Central difference
                dxdt(:,i) = (1/12*x(:,i-2) - 2/3*x(:,i-1) + 2/3*x(:,i+1) - 1/12*x(:,i+2))/h;
            end

        end
    
    end

end