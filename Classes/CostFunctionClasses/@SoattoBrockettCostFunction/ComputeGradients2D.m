% For more complex usage of the function
function [f,Df,Hf] = ComputeGradients2D(c,T,upto)
    if nargin < 3; upto = 2; end;
    T(3) = sqrt(1 - (T(1)^2 + T(2)^2));
    if ~isreal(T(3)); T(3) = 0; end;

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
    OmegaDet = (c.DetCoeffs'*p1);
    Omega = (c.OmegaCoeffs'*p1) / OmegaDet;
    H = c.HCoeffs'*p2;

    f = (T'*c.S*T - H'*Omega);
    if upto == 1; return; end;
    % Compute polynomial gradients
    % Here we use that T(3) = f(T(1),T(2)) = sqrt(1 - (T(1)^2 + T(2)^2))
    % And we use the derivatives of that:
    % dT(3)/dT(1) = -T(1)/sqrt(1 - (T(1)^2 + T(2)^2)) = -T(1)/T(3)
    % dT(3)/dT(2) = -T(2)/sqrt(1 - (T(1)^2 + T(2)^2)) = -T(2)/T(3)
    dT3dT1 = -T(1)/T(3);
    dT3dT2 = -T(2)/T(3);
    g1_1 = [  6*T(1)^5; 
              5*T(1)^4*T(2); 
              5*T(1)^4*T(3) + T(1)^5*dT3dT1; 
              4*T(1)^3*T(2)^2; 
              4*T(1)^3*T(2)*T(3) + T(1)^4*T(2)*dT3dT1; 
              4*T(1)^3*T(3)^2; 
              3*T(1)^2*T(2)^3; 
              3*T(1)^2*T(2)^2*T(3) + T(1)^3*T(2)^2*dT3dT1; 
              3*T(1)^2*T(2)*T(3)^2 + 2*T(1)^3*T(2)*T(3)*dT3dT1; 
              3*T(1)^2*T(3)^3 + 3*T(1)^3*T(3)^2*dT3dT1; 
              2*T(1)*T(2)^4; 
              2*T(1)*T(2)^3*T(3) + T(1)^2*T(2)^3*dT3dT1; 
              2*T(1)*T(2)^2*T(3)^2 + 2*T(1)^2*T(2)^2*T(3)*dT3dT1; 
              2*T(1)*T(2)*T(3)^3 + 3*T(1)^2*T(2)*T(3)^2*dT3dT1;
              2*T(1)*T(3)^4 + 4*T(1)^2*T(3)^3*dT3dT1; 
              T(2)^5; 
              T(2)^4*T(3) + T(1)*T(2)^4*dT3dT1; 
              T(2)^3*T(3)^2 + 2*T(1)*T(2)^3*T(3)*dT3dT1;
              T(2)^2*T(3)^3 + 3*T(1)*T(2)^2*T(3)^2*dT3dT1; 
              T(2)*T(3)^4 + 4*T(1)*T(2)*T(3)^3*dT3dT1; 
              T(3)^5 + 5*T(1)*T(3)^4*dT3dT1; 
              0; 
              T(2)^5*dT3dT1; 
              2*T(2)^4*T(3)*dT3dT1; 
              3*T(2)^3*T(3)^2*dT3dT1; 
              4*T(2)^2*T(3)^3*dT3dT1; 
              5*T(2)*T(3)^4*dT3dT1; 
              6*T(3)^5*dT3dT1                                      ];
    g1_2 = [  0; 
              T(1)^5; 
              T(1)^5*dT3dT2; 
              2*T(1)^4*T(2); 
              T(1)^4*T(3) + T(1)^4*T(2)*dT3dT2; 
              2*T(1)^4*T(3)*dT3dT2; 
              3*T(1)^3*T(2)^2; 
              2*T(1)^3*T(2)*T(3) + T(1)^3*T(2)^2*dT3dT2; 
              T(1)^3*T(3)^2 + 2*T(1)^3*T(2)*T(3)*dT3dT2; 
              3*T(1)^3*T(3)^2*dT3dT2; 
              4*T(1)^2*T(2)^3; 
              3*T(1)^2*T(2)^2*T(3) + T(1)^2*T(2)^3*dT3dT2; 
              2*T(1)^2*T(2)*T(3)^2 + 2*T(1)^2*T(2)^2*T(3)*dT3dT2; 
              T(1)^2*T(3)^3 + 3*T(1)^2*T(2)*T(3)^2*dT3dT2; 
              4*T(1)^2*T(3)^3*dT3dT2; 
              5*T(1)*T(2)^4; 
              4*T(1)*T(2)^3*T(3) + T(1)*T(2)^4*dT3dT2; 
              3*T(1)*T(2)^2*T(3)^2 + 2*T(1)*T(2)^3*T(3)*dT3dT2; 
              2*T(1)*T(2)*T(3)^3 + 3*T(1)*T(2)^2*T(3)^2*dT3dT2; 
              T(1)*T(3)^4 + 4*T(1)*T(2)*T(3)^3*dT3dT2; 
              5*T(1)*T(3)^4*dT3dT2; 
              6*T(2)^5; 
              5*T(2)^4*T(3) + T(2)^5*dT3dT2; 
              4*T(2)^3*T(3)^2 + 2*T(2)^4*T(3)*dT3dT2; 
              3*T(2)^2*T(3)^3 + 3*T(2)^3*T(3)^2*dT3dT2; 
              2*T(2)*T(3)^4 + 4*T(2)^2*T(3)^3*dT3dT2; 
              T(3)^5 + 5*T(2)*T(3)^4*dT3dT2; 
              6*T(3)^5*dT3dT2                                        ];

    g2_1 = [  2*T(1);
              T(2);
              T(3) + T(1)*dT3dT1;
              0;
              T(2)*dT3dT1;
              2*T(3)*dT3dT1  ];
    g2_2 = [  0
              T(1);
              T(1)*dT3dT2;
              2*T(2);
              T(3) + T(2)*dT3dT2;
              2*T(3)*dT3dT2  ];
    g1 = [g1_1 g1_2];
    g2 = [g2_1 g2_2];
    % Computing derivative of a quotient, x/y
    DOmega = ((c.OmegaCoeffs'*g1)*(c.DetCoeffs'*p1) - ...
              (c.OmegaCoeffs'*p1)*(c.DetCoeffs'*g1)) ./ ...
              (OmegaDet)^2;
    DH = c.HCoeffs'*g2;
    DT = [    1      0  ;
              0      1  ;
           dT3dT1 dT3dT2 ];
     
    Df = (2*DT'*c.S*T - (DH'*Omega + DOmega'*H));
    if upto == 2; return; end;
    % Now we compute the 2x2 Hessian. We use 
    % T(3) = sqrt(1 - (T(1)^2 - T(2)^2))
    % again, along with these:
    % d^2(T(3))/(dT(1)^2) = (T(2)^2 - 1)/(sqrt(1 - (T(1)^2 + T(2)^2))^3) ...
    %                     = (T(2)^2 - 1)/(T(3)^3)
    % d^2(T(3))/(dT(2)^2) = (T(1)^2 - 1)/(sqrt(1 - (T(1)^2 + T(2)^2))^3) ...
    %                     = (T(1)^2 - 1)/(T(3)^3)
    % d^2(T(3))/(dT(1) dT(2)) = -(T(1)*T(2))/(sqrt(1 - (T(1)^2 + T(2)^2))^3) ...
    %                         = -(T(1)*T(2))/(T(3)^3)
%     d2T3dT11 = (T(2)^2 - 1)/(T(3)^3);
%     d2T3dT22 = (T(1)^2 - 1)/(T(3)^3);
%     d2T3dT12 = -(T(1)*T(2))/(T(3)^3);

    Hf = eye(3); %  when we get around to it.
end