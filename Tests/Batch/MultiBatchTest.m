function [ results, flows, x, y ] = MultiBatchTest( types, varargin )
%BATCHTEST Test over may flow examples of a type of cost function

assert(iscell(types))
definedOrDefault = @(name,default) definedOrDefault_long(name,default,varargin);
% Define parameters
ids = definedOrDefault('ids',2:1199);
flowHandle_default = @(i) ImageFlow('image_index',i,'MaxBidirectionalError',0.01);
flowHandle = definedOrDefault('flowHandle', flowHandle_default);
nsamples = definedOrDefault('nsamples', 25);
verbose = definedOrDefault('verbose', false);


% Get positions
c = CostFunctionFactory(types{1},flowHandle(ids(1)));
[~,x,y] = c.getSurfaceResiduals(nsamples);

% Initialize results array
results_proto(1) = c.getResults;
results_proto(length(ids)) = results_proto(1,1);
results = cell(length(types),1);
for k = 1:length(types)
    results{k} = results_proto;
end

% Storing the flows
flows = cell(size(ids));

% Iterate over all the flows
t0 = CTimeleft(length(ids));
for i = 1:length(ids)
    t0.timeleft();
    id = ids(i);
    % fprintf('Image %d\n',id);
    flows{i} = flowHandle(id);
    for k = 1:length(types)
        c = CostFunctionFactory(types{k},flows{i},verbose);
        results{k}(i) = c.getResults(nsamples);
    end
end


end

