classdef SoattoBrockettCostFunction < CostFunction
%SoattoBrockettCostFunction -  Residual calcuated using the algebra from Soatto and Brockett
%   Algebraic simplifications used to compute the cost function described in the supplement. No
%   outlier rejection 
    
    properties
        X % Spherical positions
        Y % Spherical flow
        
        % Precomputation stuff
        % Omega inversion
        G11 = zeros(3);
        G12 = zeros(3);
        G13 = zeros(3);
        G22 = zeros(3);
        G23 = zeros(3);
        G33 = zeros(3);
        H11 = zeros(3,1);
        H12 = zeros(3,1);
        H13 = zeros(3,1);
        H22 = zeros(3,1);
        H23 = zeros(3,1);
        H33 = zeros(3,1);
        % Residual expantion
        S = zeros(3);
        % Coefficients for the rational functions
        OmegaCoeffs
        DetCoeffs
        HCoeffs

    end
    
    methods
        function c = SoattoBrockettCostFunction(flow)
            c@CostFunction(flow)
            [ c.X, c.Y ] = planarToSphericalFlow( flow.xy, flow.uv );
            
            X1 = [ -(c.X(2,:).^2)-(c.X(3,:).^2);
                        c.X(1,:).*c.X(2,:);
                        c.X(1,:).*c.X(3,:)      ];
            X2 = [       c.X(1,:).*c.X(2,:);
                   -(c.X(1,:).^2)-(c.X(3,:).^2);
                         c.X(2,:).*c.X(3,:)     ];
            X3 = [       c.X(1,:).*c.X(3,:);
                         c.X(2,:).*c.X(3,:);
                   -(c.X(1,:).^2)-(c.X(2,:).^2) ];
            c.G11 = X1*X1';
            c.G12 = X1*X2' + X2*X1';
            c.G13 = X1*X3' + X3*X1';
            c.G22 = X2*X2';
            c.G23 = X2*X3' + X3*X2';
            c.G33 = X3*X3';
            c.H11 = sum(bsxfun(@times,X1,c.Y(1,:)),2);
            c.H12 = sum(bsxfun(@times,X1,c.Y(2,:)),2) + ...
                            sum(bsxfun(@times,X2,c.Y(1,:)),2);
            c.H13 = sum(bsxfun(@times,X1,c.Y(3,:)),2) + ...
                            sum(bsxfun(@times,X3,c.Y(1,:)),2);
            c.H22 = sum(bsxfun(@times,X2,c.Y(2,:)),2);
            c.H23 = sum(bsxfun(@times,X2,c.Y(3,:)),2) + ...
                            sum(bsxfun(@times,X3,c.Y(2,:)),2);
            c.H33 = sum(bsxfun(@times,X3,c.Y(3,:)),2);
            c.S = c.Y*c.Y';

            % Compute coefficients for the Omega term
            % Compute the numerators of omega term
            c.OmegaCoeffs = [ SymOmegaCoefficients(c,2,3) ...
                              SymOmegaCoefficients(c,3,1) ...
                              SymOmegaCoefficients(c,1,2) ];

            % Compute coefficients of the determinant (denominator of Omega term)
            c.DetCoeffs = DeterminantCoefficients(c);
            
            % The 'H' coeffcients
            c.HCoeffs = [ c.H11 ...
                          c.H12 ...
                          c.H13 ...
                          c.H22 ...
                          c.H23 ...
                          c.H33 ]';
            
        end
        function [resid] = getResidual(c,T)
            p1 = [  T(1)^6; 
                    T(1)^5*T(2); 
                    T(1)^5*T(3); 
                    T(1)^4*T(2)^2; 
                    T(1)^4*T(2)*T(3); 
                    T(1)^4*T(3)^2; 
                    T(1)^3*T(2)^3; 
                    T(1)^3*T(2)^2*T(3); 
                    T(1)^3*T(2)*T(3)^2; 
                    T(1)^3*T(3)^3; 
                    T(1)^2*T(2)^4; 
                    T(1)^2*T(2)^3*T(3); 
                    T(1)^2*T(2)^2*T(3)^2; 
                    T(1)^2*T(2)*T(3)^3; 
                    T(1)^2*T(3)^4; 
                    T(1)*T(2)^5; 
                    T(1)*T(2)^4*T(3); 
                    T(1)*T(2)^3*T(3)^2; 
                    T(1)*T(2)^2*T(3)^3; 
                    T(1)*T(2)*T(3)^4; 
                    T(1)*T(3)^5; 
                    T(2)^6; 
                    T(2)^5*T(3); 
                    T(2)^4*T(3)^2; 
                    T(2)^3*T(3)^3; 
                    T(2)^2*T(3)^4; 
                    T(2)*T(3)^5; 
                    T(3)^6         ];
            p2 = [  T(1)^2;
                    T(1)*T(2);
                    T(1)*T(3);
                    T(2)^2;
                    T(2)*T(3);
                    T(3)^2  ];
            Omega = (c.OmegaCoeffs'*p1) /(c.DetCoeffs'*p1);
            H = c.HCoeffs'*p2;

            resid = (T'*c.S*T - H'*Omega); % Implicitly divided by T'*T
        end
        
        function [Omega] = getOmega(c,T)
            p1 = [  T(1)^6; 
                    T(1)^5*T(2); 
                    T(1)^5*T(3); 
                    T(1)^4*T(2)^2; 
                    T(1)^4*T(2)*T(3); 
                    T(1)^4*T(3)^2; 
                    T(1)^3*T(2)^3; 
                    T(1)^3*T(2)^2*T(3); 
                    T(1)^3*T(2)*T(3)^2; 
                    T(1)^3*T(3)^3; 
                    T(1)^2*T(2)^4; 
                    T(1)^2*T(2)^3*T(3); 
                    T(1)^2*T(2)^2*T(3)^2; 
                    T(1)^2*T(2)*T(3)^3; 
                    T(1)^2*T(3)^4; 
                    T(1)*T(2)^5; 
                    T(1)*T(2)^4*T(3); 
                    T(1)*T(2)^3*T(3)^2; 
                    T(1)*T(2)^2*T(3)^3; 
                    T(1)*T(2)*T(3)^4; 
                    T(1)*T(3)^5; 
                    T(2)^6; 
                    T(2)^5*T(3); 
                    T(2)^4*T(3)^2; 
                    T(2)^3*T(3)^3; 
                    T(2)^2*T(3)^4; 
                    T(2)*T(3)^5; 
                    T(3)^6         ];
            Omega = (c.OmegaCoeffs'*p1) /(c.DetCoeffs'*p1);
        end
        
        function [flowResids] = getFlowResiduals(c,T)
            p1 = [  T(1)^6; 
                    T(1)^5*T(2); 
                    T(1)^5*T(3); 
                    T(1)^4*T(2)^2; 
                    T(1)^4*T(2)*T(3); 
                    T(1)^4*T(3)^2; 
                    T(1)^3*T(2)^3; 
                    T(1)^3*T(2)^2*T(3); 
                    T(1)^3*T(2)*T(3)^2; 
                    T(1)^3*T(3)^3; 
                    T(1)^2*T(2)^4; 
                    T(1)^2*T(2)^3*T(3); 
                    T(1)^2*T(2)^2*T(3)^2; 
                    T(1)^2*T(2)*T(3)^3; 
                    T(1)^2*T(3)^4; 
                    T(1)*T(2)^5; 
                    T(1)*T(2)^4*T(3); 
                    T(1)*T(2)^3*T(3)^2; 
                    T(1)*T(2)^2*T(3)^3; 
                    T(1)*T(2)*T(3)^4; 
                    T(1)*T(3)^5; 
                    T(2)^6; 
                    T(2)^5*T(3); 
                    T(2)^4*T(3)^2; 
                    T(2)^3*T(3)^3; 
                    T(2)^2*T(3)^4; 
                    T(2)*T(3)^5; 
                    T(3)^6         ];
            Omega = (c.OmegaCoeffs'*p1) /(c.DetCoeffs'*p1);

            flowResids = zeros(c.flow.nPoints,1);
            for i = 1:c.flow.nPoints
                x = c.X(:,i);
                Xhat = [     0, -x(3),  x(2);
                          x(3),     0, -x(1);
                         -x(2),  x(1),     0 ];
                flowResids(i) = (T'*((Xhat^2)*Omega - c.Y(:,i)))^2;
            end
        end
        
        [f,Df,Hf] = ComputeGradients(residual,T,upto)
        
        function plotSphericalGradients(c,nsamples,scale)
            if nargin < 2; nsamples = 30; end
            if nargin < 3; scale = 10; end
            
            % Get the residuals
            [z,x,y] = c.getSurfaceResiduals(nsamples);

            % 3d plot
            gx = zeros(size(x));
            gy = zeros(size(x));
            gz = zeros(size(x));
            for i = 1:length(x(:))
                T = [ x(i); y(i); sqrt(1 - min(1,(x(i).^2 + y(i).^2))) ];
                [~, g_i] = c.ComputeGradients(T);
                gx(i) = g_i(1);
                gy(i) = g_i(2);
                gz(i) = g_i(3);
            end

            surf(x,y,sqrt(1 - min(1,(x.^2 + y.^2))),z,'LineStyle','none')
            hold on
            scatter3(c.trueT(1),c.trueT(2),c.trueT(3),600,'g.');
            quiver3(...
                x,y,sqrt(1 - min(1,(x.^2 + y.^2))),...
                -scale*gx,-scale*gy,-scale*gz,...
                'Autoscale','off')
            axis equal
            view(2)
            hold off
        end
        % DEPRECATED - using gradients on (x,y) coordinates instead of (x,y,z) coordinates
        function plotPlanarGradients(c,nsamples,scale)
            if nargin < 2; nsamples = 30; end
            if nargin < 3; scale = 0.1; end
            
            % Get the residuals
            [z,x,y] = c.getSurfaceResiduals(nsamples);

            % 2d plot
            gx = zeros(size(x));
            gy = zeros(size(y));
            for i = 1:length(x(:))
                T = [ x(i); y(i); 0 ];
                if norm(T) < 0.9
                    [~, g] = c.ComputeGradients2D(T);
                    gx(i) = g(1);
                    gy(i) = g(2);
                else
                    gx(i) = 0;
                    gy(i) = 0;
                end
            end

            figure
            contourf(x,y,z,30,'LineStyle','none');
            colormap('parula')
            axis equal
            hold on
            % Get our estimate
            [~,t_min_guess] = min(z(:));

            scatter(c.trueT(1),c.trueT(2),600,'g.');
            scatter(x(t_min_guess),y(t_min_guess),600,'k.');

            quiver(x,y,-scale*gx,-scale*gy,'g','Autoscale','off')
            hold off

        end
    end

    methods (Access = protected)
        function [ X, Y ] = planarToSphericalFlow( xy, uv )
        %PLANARTOSPHERICALFLOW Converts planar flow into spherical flow
        % Inputs:
        % params - paramters of the current session
        % xy - size (2 x n), calibrated positions of flow on the plane
        % uv - size (2 x n), corresponding calibrated flow on the plane
        % Outputs:
        % X - size (3 x n), spherical positions of the flow of corresponding planar
        %     flow
        % Y - size (2 x n), spherical flow, tangent to the corresponding points X
        % 
            xy = double(xy);
            uv = double(uv);
            n = size(xy,2);
            X = [ xy; ones(1,n) ];
            norms = repmat(sqrt(sum(X.^2)),3,1);
            X = X ./ norms;
            Y = [               uv(2,:)             ;
                               -uv(1,:)             ;
                  xy(2,:).*uv(1,:) - xy(1,:).*uv(2,:) ] ./ (norms.^2);
        end
        
        OmegaCoeffs = SymOmegaCoefficients(coeffs,i,j)

        DetCoeffs = DeterminantCoefficients(coeffs)
    end
end

