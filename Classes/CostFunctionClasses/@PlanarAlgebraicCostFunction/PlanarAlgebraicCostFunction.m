classdef PlanarAlgebraicCostFunction < CostFunction
%PlanarAlgebraicCostFunction - Fast computation of heading using algebraic simplifications
%   The algebraic simplifications are given in the supplemental material of the paper. No outlier
%   rejection is used here
    
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
        OmegaCoeffs1
        OmegaCoeffs2
        OmegaCoeffs3
        DetCoeffs
        HCoeffs1
        HCoeffs2
        HCoeffs3

    end
    
    methods
        function c = PlanarAlgebraicCostFunction(flow)
            c@CostFunction(flow)
            
            x = flow.xy(1,:)';
            y = flow.xy(2,:)';
            uv = flow.uv;
            % u = uv(1,:)';
            % v = uv(2,:)';

            % First the dumb way: Compute the coefficients every time.
            A = @(x,y) [ -1  0  x;
                          0 -1  y ];
            B = @(x,y) [   x*y  , -(1+x^2),  y;
                         (1+y^2),   -x*y  , -x ];
            J = [ 0 1;
                 -1 0 ];
            BJA = @(x,y)  [ y^2 + 1,  -x*y  ,    -x    ;
                             -x*y  , x^2 + 1,    -y    ;
                              -x   ,    -y  , x^2 + y^2  ];

            c.G11 = zeros(3);
            c.G12 = zeros(3);
            c.G13 = zeros(3);
            c.G22 = zeros(3);
            c.G23 = zeros(3);
            c.G33 = zeros(3);
            c.H11 = zeros(3,1);
            c.H12 = zeros(3,1);
            c.H13 = zeros(3,1);
            c.H22 = zeros(3,1);
            c.H23 = zeros(3,1);
            c.H33 = zeros(3,1);

            for i = 1:flow.nPoints
                BJAx = BJA(x(i),y(i));
                AJv = A(x(i),y(i))'*J'*uv(:,i);
                % Invertion matrix for Omega
                c.G11 = c.G11 + BJAx*[1 0 0;
                                      0 0 0;
                                      0 0 0]*BJAx;
                c.G12 = c.G12 + BJAx*[0 1 0;
                                      1 0 0;
                                      0 0 0]*BJAx;
                c.G13 = c.G13 + BJAx*[0 0 1;
                                      0 0 0;
                                      1 0 0]*BJAx;
                c.G22 = c.G22 + BJAx*[0 0 0;
                                      0 1 0;
                                      0 0 0]*BJAx;
                c.G23 = c.G23 + BJAx*[0 0 0;
                                      0 0 1;
                                      0 1 0]*BJAx;
                c.G33 = c.G33 + BJAx*[0 0 0;
                                      0 0 0;
                                      0 0 1]*BJAx;
                % Vector for Omega
                c.H11 = c.H11 + BJAx*[1 0 0;
                                      0 0 0;
                                      0 0 0]*AJv;
                c.H12 = c.H12 + BJAx*[0 1 0;
                                      1 0 0;
                                      0 0 0]*AJv;
                c.H13 = c.H13 + BJAx*[0 0 1;
                                      0 0 0;
                                      1 0 0]*AJv;
                c.H22 = c.H22 + BJAx*[0 0 0;
                                      0 1 0;
                                      0 0 0]*AJv;
                c.H23 = c.H23 + BJAx*[0 0 0;
                                      0 0 1;
                                      0 1 0]*AJv;
                c.H33 = c.H33 + BJAx*[0 0 0;
                                      0 0 0;
                                      0 0 1]*AJv;
                % Quadratic matrix
                c.S = c.S + AJv*AJv';
            end

            % Finish computing the coefficients
            c.OmegaCoeffs1 = SymOmegaCoefficients(c,2,3);
            c.OmegaCoeffs2 = SymOmegaCoefficients(c,3,1);
            c.OmegaCoeffs3 = SymOmegaCoefficients(c,1,2);

            % Compute coefficients of the determinant
            c.DetCoeffs = DeterminantCoefficients(c);
            HCoeffs = [ c.H11 ...
                        c.H12 ...
                        c.H13 ...
                        c.H22 ...
                        c.H23 ...
                        c.H33 ]';
            size(HCoeffs)
            c.HCoeffs1 = HCoeffs(:,1);
            c.HCoeffs2 = HCoeffs(:,2);
            c.HCoeffs3 = HCoeffs(:,3);

        end
        function [resid] = getResidual(c,T)
        % Takes in a 3 x N vector
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
            Omega = [ (p1'*c.OmegaCoeffs1);
                      (p1'*c.OmegaCoeffs2);
                      (p1'*c.OmegaCoeffs3) ]/(p1'*c.DetCoeffs);
            H = ([ p2'*c.HCoeffs1;
                   p2'*c.HCoeffs2;
                   p2'*c.HCoeffs3 ]);

            resid = (T'*c.S*T - H'*Omega);
        end
        function [Omega] = getOmega(c,T)
        % Takes in a 3 x N vector
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
            Omega = [ (p1'*c.OmegaCoeffs1);
                      (p1'*c.OmegaCoeffs2);
                      (p1'*c.OmegaCoeffs3) ]/(p1'*c.DetCoeffs);
        end
    end
    
    methods (Access = private)
        OmegaCoeffs = SymOmegaCoefficients(coeffs,i,j)

        DetCoeffs = DeterminantCoefficients(coeffs)
    end
end

