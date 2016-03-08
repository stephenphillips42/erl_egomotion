classdef RobustERLCostFunction < CostFunction
%RobustERLCostFunction - This implements the Expected Residual Likelihood to estimate heading
%   This is the main result of the paper
    
    % Side note for this implementation we will
    properties
        weights % Weights of the individual flow vectors
        Aperp
        B
        u
    end
    
    methods
        function c = RobustERLCostFunction(flow)
            c@CostFunction(flow)
            % For ease of reading
            n = c.flow.nPoints;
            c.weights = ones(n,1);
            J = [0 1; -1 0];
            Jxy = J*c.flow.xy;
            c.Aperp = [ -repmat(J,n,1) Jxy(:) ];
            col1 = reshape([c.flow.xy(1,:).*c.flow.xy(2,:);
                                     1 + c.flow.xy(2,:).^2 ],[],1);
            col2 = reshape([         -(1+c.flow.xy(1,:).^2);
                            -c.flow.xy(1,:).*c.flow.xy(2,:) ],[],1);
            c.B = [col1 col2 Jxy(:)];
            c.u = c.flow.uv(:);
            % Use our initial outlier rejection scheme as initialization of
            % the weights
            c.weights = c.outlierRejection;
            c.weights = (c.weights - min(c.weights));
            c.weights = (c.weights)/max(c.weights);
        end
        
        function f = getFlowResiduals(c,V)
            n = c.flow.nPoints;
            A1 = normc(reshape(c.Aperp*V,2,[]));
            i = 1:(2*n);
            j = repmat(1:n,2,1);
            AperpVT = sparse(j(:),i',A1(:));
            A = bsxfun(@times, AperpVT*c.B, c.weights);
            b = bsxfun(@times, AperpVT*c.u, c.weights);
            Omega = (A'*A) \ (A'*b);
            f = (A*Omega - b).^2;
        end

        function [Omega] = getOmega(c,V)
            n = c.flow.nPoints;
            A1 = normc(reshape(c.Aperp*V,2,[]));
            i = 1:(2*n);
            j = repmat(1:n,2,1);
            AperpVT = sparse(j(:),i',A1(:));
            A = bsxfun(@times, AperpVT*c.B, c.weights);
            b = bsxfun(@times, AperpVT*c.u, c.weights);
            G = (A'*A);
            H = (A'*b);
            Omega = (G\H);
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
            result = getResults@CostFunction(c,nsamples);
            result.extraData = {c.weights};
        end
    end

end

