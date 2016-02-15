function [ estT, estOmega, trueT, trueOmega ] = LightMultiBatchTest( types, varargin )
%BATCHTEST Test over may flow examples of a type of cost function

assert(iscell(types))
definedOrDefault = @(name,default) definedOrDefault_long(name,default,varargin);
% Define parameters
ids = definedOrDefault('ids',2:1199);
flowHandle_default = @(i) ImageFlow('image_index',i,'MaxBidirectionalError',0.01);
flowHandle = definedOrDefault('flowHandle', flowHandle_default);
nsamples = definedOrDefault('nsamples', 25);
verbose = definedOrDefault('verbose', false);

% Initialize results arrays
estT = zeros(length(types),length(ids),3);
estOmega = zeros(length(types),length(ids),3);
trueT = zeros(3,length(ids));
trueOmega = zeros(3,length(ids));

% Iterate over all the flows
t0 = CTimeleft(length(ids));
for i = 1:length(ids)
    t0.timeleft();
    id = ids(i);
    % fprintf('Image %d\n',id);
    flow = flowHandle(id);
    trueT(:,i) = flow.trueT;
    trueOmega(:,i) = flow.trueOmega;
    for k = 1:length(types)
        c = CostFunctionFactory(types{k},flow,verbose);
        result = c.getResults(nsamples);
        estT(k,i,:) = result.minimumT;
        estOmega(k,i,:) = result.minimumOmega;
    end
end


end

