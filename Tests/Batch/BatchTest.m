function [ results, flows, x, y ] = BatchTest( type, varargin )
%BATCHTEST Test over may flow examples of a type of cost function

definedOrDefault = @(name,default) definedOrDefault_long(name,default,varargin);
% Define parameters
ids = definedOrDefault('ids',2:1199);
flowHandle_default = @(i) ImageFlow('image_index',i,'MaxBidirectionalError',0.01);
flowHandle = definedOrDefault('flowHandle', flowHandle_default);
nsamples = definedOrDefault('nsamples', 25);
verbose = definedOrDefault('verbose', true);



% Storing the flows
flows = cell(size(ids));

% Get positions
c = CostFunctionFactory(type,flowHandle(ids(1)));
[~,x,y] = c.getSurfaceResiduals(nsamples);
% Build results struct
results_proto = c.getResults;
results(1) = results_proto;
results(length(ids)) = results_proto;

% Iterate over all the flows
for i = 1:length(ids)
    fprintf('Image %d\n',ids(i));
    flows{i} = flowHandle(ids(i));
    c = CostFunctionFactory(type,flows{i});
    c.verbose = verbose;
    % Store results
    a = c.getResults(nsamples);
    results(i) = a;
end


end

