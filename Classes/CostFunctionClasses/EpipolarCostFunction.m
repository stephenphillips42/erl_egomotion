classdef EpipolarCostFunction < CostFunction
%EpipolarCostFunction - Using Epipolar RANSAC to estimate heading direction
%   Uses Peter Corke's Robotics tool box as implementation of RANSAC
    
    properties
        p % points
    end
    
    methods
        function c = EpipolarCostFunction(flow)
            c@CostFunction(flow);
            c.p = [ c.flow.xy;
                    c.flow.xy + c.flow.uv ];
        end
        
        % TODO: Make this some optimization (just for fun)
        % Residual for each flow vector
        function [flowResids] = getFlowResiduals(c,T)
            flowResids = zeros(c.flow.nPoints,min(size(T)));
        end
        
        function result = getResults(c,~)
            % Run epipolar nonsense
            [F, inliers, r] = ransac(@fmatrix, c.p, 0.00001);
            % Get Candidate Omega and T
            [U,~,V] = svd(F);
            RZ = [ 0 -1  0;
                   1  0  0;
                   0  0  1 ];
            R1 = U*RZ'*V';
            R1 = sign(det(R1))*R1;
            R2 = U*RZ*V';
            R2 = sign(det(R2))*R2;
            T1 = U(:,3);
            
            % Get T
            finalT = T1*sign(T1(3));
            angError = acosd(sign(finalT'*c.flow.trueT)*min(abs(finalT'*c.flow.trueT),1));
            
            % Get Omega
            q1 = [c.flow.xy(:,inliers(1)); 1];
            q2 = [c.flow.xy(:,inliers(1)) + c.flow.uv(:,inliers(1)); 1];
            QQ1 = [ -R1*q1 q2 ];
            QQ2 = [ -R2*q1 q2 ];
            finalOmega = [0;0;0];
            if (all((QQ1'*QQ1)\(QQ1'*T1) > 0) || all((QQ1'*QQ1)\(QQ1'*T1) < 0))
                omega_hat = (R1-R1')/2;
                finalOmega = [ omega_hat(3,2);
                               omega_hat(1,3);
                               omega_hat(2,1) ]; % Extracting skew symmetric 
            elseif (all((QQ2'*QQ2)\(QQ2'*T1) > 0) || all((QQ2'*QQ2)\(QQ2'*T1) < 0))
                omega_hat = (R2-R2')/2;
                finalOmega = [ omega_hat(3,2);
                               omega_hat(1,3);
                               omega_hat(2,1) ]; % Extracting skew symmetric 
            else
                disp('Woah woah woah no omega???')
                disp((QQ2'*QQ2)\(QQ2'*T1))
                disp((QQ1'*QQ1)\(QQ1'*T1))
            end
            
            % Save out results
            result.trueT = c.trueT;
            result.minimumOmega = finalOmega;
            result.angularError = angError;
            result.minimumT = finalT;
            result.minimumResidual = r;
            result.extraData = {inliers};
        end
    end
    
end

