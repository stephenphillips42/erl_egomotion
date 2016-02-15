function [results, flows]= ParameterBatchTest( costFunctions, varargin )
%PARAMETERBATCHTEST Test over may flow examples with varying parametrs

assert(iscell(costFunctions))
assert(~isempty(costFunctions))
assert(isa(costFunctions{1},'function_handle'))
definedOrDefault = @(name,default) definedOrDefault_long(name,default,varargin);
% Define parameters
ids = definedOrDefault('ids',2:1199);
flowHandle_default = @(i) ImageFlow('image_index',i,'MaxBidirectionalError',0.01);
flowHandle = definedOrDefault('flowHandle', flowHandle_default);
nsamples = definedOrDefault('nsamples', 25);
verbose = definedOrDefault('verbose', false);

% Initialize results array
c = costFunctions{1}(flowHandle(ids(1)));
results_proto = c.getResults;
results_proto(length(ids)) = results_proto(1,1);
results = cell(length(costFunctions),1);
for k = 1:length(costFunctions)
    results{k} = results_proto;
end

% Store the flows
flows = cell(length(ids));

% Iterate over all the flows
t0 = CTimeleft(length(ids));
for i = 1:length(ids)
    t0.timeleft();
    id = ids(i);
    flows{i} = flowHandle(id);
    for k = 1:length(costFunctions)
        c = costFunctions{k}(flows{i});
        c.verbose = verbose;
        results{k}(i) = c.getResults(nsamples);
    end
end



end

