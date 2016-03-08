function c = CostFunctionFactory( type, flow, verbose )
%CostFunctionFactory Used to produce a cost function abstractly. Useful in batch computations
% comparing different types of residual computations.
% Input:
% flow - instance of sub-class of the Flow class from FlowClasses
% verbose - boolean, specifies if the cost function prints all output (optional, default false)        

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

