function c = CostFunctionFactory( type, flow, verbose )
%COSTFUNCTIONFACTORY Summary of this function goes here
%   Detailed explanation goes here
% 

if strcmp(type,'SoattoBrocket')
    c = SoattoBrocketCostFunction(flow);
elseif strcmp(type,'Epipolar')
    c = EpipolarCostFunction(flow);
elseif strcmp(type,'ZhangTomasi')
    c = ZhangTomasiCostFunction(flow);
elseif strcmp(type,'PlanarAlgebraicCostFunction')
    c = PlanarAlgebraicCostFunctionCostFunction(flow);
elseif strcmp(type,'RobustERL')
    c = RobustERLCostFunction(flow);
elseif strcmp(type,'RobustLifted')
    c = RobustLiftedCostFunction(flow);
else
    error('Error: Unknown cost function type')
end

if nargin == 3
    c.verbose = verbose;
end

end

