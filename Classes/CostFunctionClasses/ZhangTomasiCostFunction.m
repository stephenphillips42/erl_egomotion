classdef ZhangTomasiCostFunction < CostFunction
%ZhangTomasiCostFunction - Uses formulation from Zhang/Tomasi cited in paper
%   No outlier rejection used in this scheme
    
    properties
        B % Matrix coefficients (n x 1 cell array of 3 x 3 matrices)
        Aperp % Flow coefficients (3 x n matrix)
        u % Denominator matrix cofficients (n x 1 cell array of 3 x 3 matrices)
    end
    
    methods
        function c = ZhangTomasiCostFunction(flow)
            c@CostFunction(flow);
            n = c.flow.nPoints;
            J = [0 1; -1 0];
            Jxy = J*c.flow.xy;
            c.Aperp = [ -repmat(J,n,1) Jxy(:) ];
            col1 = reshape([c.flow.xy(1,:).*c.flow.xy(2,:);
                                     1 + c.flow.xy(2,:).^2 ],[],1);
            col2 = reshape([         -(1+c.flow.xy(1,:).^2);
                            -c.flow.xy(1,:).*c.flow.xy(2,:) ],[],1);
            c.B = [col1 col2 Jxy(:)];
            c.u = c.flow.uv(:);
        end
        function [flowResids] = getFlowResiduals(c,V)
            flowResids = allResiduals(c,V);
        end
        function [Omega] = getOmega(c,V)
            n = c.flow.nPoints;
            A1 = normc(reshape(c.Aperp*V,2,[]));
            i = 1:(2*n);
            j = repmat(1:n,2,1);
            AperpVT = sparse(j(:),i',A1(:));
            A = AperpVT*c.B;
            b = AperpVT*c.u;
            Omega = (A.'*A)\(A.'*b);
        end
    end
    
    methods (Access = private)
        function [errs] = allResiduals(c,V)
            % Compute Omega, based on Aperp*(B*Omega - v) least squares
            n = c.flow.nPoints;
            A1 = normc(reshape(c.Aperp*V,2,[]));
            i = 1:(2*n);
            j = repmat(1:n,2,1);
            AperpVT = sparse(j(:),i',A1(:));
            A = AperpVT*c.B;
            b = AperpVT*c.u;
            Omega = (A.'*A)\(A.'*b);
            errs = (A*Omega - b).^2;
        end
    end
    
end

