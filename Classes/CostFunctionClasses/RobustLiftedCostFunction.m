classdef RobustLiftedCostFunction < ZhangMatrixCostFunction
%RobustLiftedCostFunction - This computes the heading direction using the lifted kernel formulation
%   This is the main method we compare ERL to
    
    % Side note for this implementation we will
    properties
        weights % Weights of the individual flow vectors
        weights_initalized % Boolean to check if we should use weights
        maxiters         % Max iterationsn before we stop
        scale            % How much we scale lambda by at each step
        lambda_max       % Minimum trust region radius
        lambda_min       % Maximum trust region radius
        min_step_quality % Minimum deviation from our predicted value we allow
        step_tol         % How small our step has to be before we stop
        grad_tol         % How small our gradient has to be before we stop
        tau              % Cost function width parameter
        lambda0          % Initial value
    end
    
    methods
        function c = RobustLiftedCostFunction(flow,varargin)
            c@ZhangMatrixCostFunction(flow)
            
            % Parameters
            definedOrDefault = @(name,value) definedOrDefault_long(name,value,varargin);
            % Set parameters
            
            c.maxiters = definedOrDefault('maxiters',700);
            c.scale = definedOrDefault('scale',4);
            c.lambda_max = definedOrDefault('lambda_max',1e8);
            c.lambda_min = definedOrDefault('lambda_min',1e-12); 
            c.min_step_quality = definedOrDefault('min_step_quality',1e-3);
            c.step_tol = definedOrDefault('step_tol',10^-7);
            c.grad_tol = definedOrDefault('grad_tol',10^-9);
            c.tau = definedOrDefault('tau',0.05);
            c.lambda0 = definedOrDefault('lambda0',0.1);
            c.verbose = definedOrDefault('verbose',true);

            % Use our initial outlier rejection scheme as initialization of
            % the weights
            c.weights_initalized = false;
            c.weights = ones(c.flow.nPoints,1);
            c.weights = c.outlierRejection;
            c.weights = c.weights - min(c.weights);
            c.weights = c.weights/max(c.weights) + 0.01;
            c.weights_initalized = true;
        end
        
        function f = getFlowResiduals(c,V)
            if c.weights_initalized
                f = c.getOmegaAndWeights(V);
            else
                f = c.getFlowResiduals@ZhangMatrixCostFunction(V);
            end
        end
        
        function [Omega] = getOmega(c,V)
            [~, Omega, ~] = c.getOmegaAndWeights(V);
        end
        
        function [r, Omega, w] = getOmegaAndWeights(c,V)
        % Solve for the omega and the weights for the given V
            % Create the matrices
            m = c.flow.nPoints;
            A1 = normc(reshape(c.Aperp*V,2,[]));
            i = 1:(2*m);
            j = repmat(1:m,2,1);
            AperpV = sparse(i',j(:),A1(:));
            A = AperpV'*c.B;
            b = AperpV'*c.u;
            [r, Omega, w] = c.liftedWeightsOptimization(A,b,ones(m,1));
        end
        
        function [r, x, w] = liftedWeightsOptimization(c,A,b,w0)
        % Optimize the weights using levenberg-marquart
            % Parameters
            [m,n] = size(A);
            lambda = c.lambda0;
            % Helper functions
            % f = @(x,w) sum(((A*x - b).^2).*(w.^2) + ((c.tau^2)/2)*(1-w.^2).^2);
            I = eye(n);
            % Initialization
            w = w0;
            W2 = sparse(1:m,1:m,w.^2);
            x = (A'*W2*A + lambda*I)\(A'*W2*b);
            % Begin optimization
            f_prev = sum(((A*x - b).^2).*(w.^2) + ((c.tau^2)/2)*(1-w.^2).^2);
            updated_xw = true;
            for it = 1:c.maxiters
                % Build residual
                r = A*x - b;
                % If we need to update the jacobian
                if updated_xw
                    % Build diagonal Jacobian matrices
                    k = (c.tau/sqrt(2))*(w.^2 - 1);
                    DK = sqrt(2)*c.tau*w;
                    z = r.^2 + DK.^2 + lambda;
                    % Build diagonal weighing matrices
                    d1 = (w.^2).*(1-(r.^2)./z);
                    d2 = w.*((1-(r.^2)./z) + (sqrt(2)*c.tau)*w./z);
                    D1 = sparse(1:m, 1:m, d1);
                    D2 = sparse(1:m, 1:m, d2);
                    % Solve for change in variables and weights using Schur
                    % complements
                end
                dx = (A'*D1*A + lambda*I) \ (-A'*D2*r);
                dw = -w.*(r.*(r + A*dx) + DK.*k)./z;
                x_next = x + dx;
                w_next = w + dw;
                % Check if we need to change trust region
                f_next = sum(((A*x_next - b).^2).*(w_next.^2) + ...
                                        ((c.tau^2)/2)*(1-w_next.^2).^2);
                f_pred = ([w.*(A*dx) + r.*dw + w.*r; DK.*dw + k]);
                f_pred = f_pred'*f_pred;

                if f_next > f_prev || ...
                        (f_next - f_prev)/(f_pred - f_prev) < c.min_step_quality
                    updated_xw = false;
                    lambda = min((c.scale)*lambda, c.lambda_max);
                    if lambda == c.lambda_max
                        break
                    end
                    continue;
                end
                updated_xw = true;
                lambda = max((1/c.scale)*lambda, c.lambda_min);
                f_prev = f_next;
                x = x_next;
                w = w_next;
                % Check for stopping criterion (2 types)
                W2 = sparse(1:m,1:m,w.^2);
                if max(abs([A'*W2*r; w.*(r.^2) + DK.*k])) < c.grad_tol
                    break
                end
                if norm([dx; dw]) < c.step_tol*(sqrt(eps) + norm([x; w]))
                    break
                end
            end
            r = ((A*x - b).^2).*(w.^2) + ((c.tau^2)/2)*(1-w.^2).^2;
        end
        
        function [rho] = getInverseDepths(c,T,Omega)
            if nargin < 3
                Omega = c.getOmega;
            end
            n = c.flow.nPoints;
            A = [ repmat(eye(2),n,1) c.flow.xy(:) ]*T;
            i = 1:(2*n);
            j = repmat(1:n,2,1);
            AT = sparse(i',j(:),A(:));
            rho = -(AT'*AT) \ (AT'*(c.B*Omega - c.u));
        end
        
        % For brevity in batch computations
        function result = getResults(c,nsamples)
            if nargin < 2; nsamples = 25; end;
            [z, x, y] = c.getSurfaceResiduals(nsamples); 
            [~, mi] = min(z(:));
            initT = [x(mi); y(mi); sqrt(max(0,1 - (x(mi)^2 + y(mi)^2)))];
            finalT = c.gradientDescent(initT);
            angError = acosd(sign(finalT'*c.trueT)*min(abs(finalT'*c.trueT),1));
            [r, Omega_opt, w_opt] = c.getOmegaAndWeights(finalT);
            
            result.trueT = c.trueT;
            result.angularError = angError;
            result.minimumT = finalT;
            result.minimumOmega = Omega_opt; 
            result.minimumResidual = r;
            result.extraData = {w_opt};
        end

        
        function [r, x, w, f_hist, x_hist, w_hist, r_hist] = ...
                                    liftedWeightsOptimizationDebug(c,A,b,w0)
            % Parameters
            [m,n] = size(A);
            lambda = c.lambda0;
            % Debug arrays
            f_hist = [];
            x_hist = [];
            w_hist = [];
            r_hist = [];
            % Helper functions
            % f = @(x,w) sum(((A*x - b).^2).*(w.^2) + ((c.tau^2)/2)*(1-w.^2).^2);
            I = speye(n);
            % Initialization
            w = w0;
            W2 = sparse(1:m,1:m,w.^2);
            x = (A'*W2*A + lambda*speye(3))\(A'*W2*b);
            % Begin optimization
            f_prev = sum(((A*x - b).^2).*(w.^2) + ((c.tau^2)/2)*(1-w.^2).^2);
            updated_xw = true;
            for it = 1:c.maxiters
                % Build residual
                r = A*x - b;
                % If we need to update the jacobian
                if updated_xw
                    % Build diagonal Jacobian matrices
                    k = (c.tau/sqrt(2))*(w.^2 - 1);
                    DK = sqrt(2)*c.tau*w;
                    z = r.^2 + DK.^2 + lambda;
                    % Build diagonal weighing matrices
                    d1 = (w.^2).*(1-(r.^2)./z);
                    d2 = w.*((1-(r.^2)./z) + (sqrt(2)*c.tau)*w./z);
                    D1 = sparse(1:m, 1:m, d1);
                    D2 = sparse(1:m, 1:m, d2);
                    % Solve for change in variables and weights using Schur
                    % complements
                end
                dx = (A'*D1*A + lambda*I) \ (-A'*D2*r);
                dw = -w.*(r.*(r + A*dx) + DK.*k)./z;
                x_next = x + dx;
                w_next = w + dw;
                % Get the debug hist
                x_hist = [x_hist x];
                w_hist = [w_hist w];
                f_hist = [f_hist f_prev];
                r_hist = [r_hist r];
                % Check if we need to change trust region
                f_next = sum(((A*x_next - b).^2).*(w_next.^2) + ...
                                        ((c.tau^2)/2)*(1-w_next.^2).^2);
                f_pred = ([w.*(A*dx) + r.*dw + w.*r; DK.*dw + k]);
                f_pred = f_pred'*f_pred;
                fprintf('------------------ %d --------------------\n',it)
                fprintf('f value: %3.2e\n',f_prev)
                fprintf('f_next value: %3.2e\n', f_next)
                fprintf('f_pred value: %3.2e\n', f_pred)
                fprintf('lambda: %3.2e\n', lambda)
                fprintf('Values: %3.2e vs %3.2e\n',f_next, f_prev)
                fprintf('Standard quality: %3.2e vs %3.2e\n',...
                        (f_next - f_prev)/(f_pred - f_prev), c.min_step_quality)

                if f_next > f_prev || ...
                        (f_next - f_prev)/(f_pred - f_prev) < c.min_step_quality
                    updated_xw = false;
                    lambda = min((c.scale)*lambda, c.lambda_max);
                    if lambda == c.lambda_max
                        disp('max lambda stop')
                        break
                    end
                    fprintf('Did not take step\n')
                    continue;
                end
                updated_xw = true;
                lambda = max((1/c.scale)*lambda, c.lambda_min);
                f_prev = f_next;
                x = x_next;
                w = w_next;
                % Check for stopping criterion (2 types)
                W2 = sparse(1:m,1:m,w.^2);
                if norm([A'*W2*r; w.*(r.^2) + DK.*k], Inf) < c.grad_tol
                    disp('gradient stop')
                    break
                end
                if norm([dx; dw]) < c.step_tol*(sqrt(eps) + norm([x; w]))
                    disp('norm stop')
                    break
                end
            end
            r = ((A*x - b).^2).*(w.^2) + ((c.tau^2)/2)*(1-w.^2).^2;
        end
    end
    
    methods(Access = private)
        function [r, J] = mycost(~,input,X,y)
            theta = input(1:size(X,2));
            w = input((size(X,2)+1):end);
            r = [ w.*(y - X*theta);
                     (1-w.^2)      ];
            J = [    diag(w)*X   diag(y - X*theta);
                  zeros(size(X))    2*diag(w)      ];
        end
    end

end

