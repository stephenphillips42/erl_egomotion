function [ results, flows, x, y ] = RandomizedBatchTest( type, varargin )
%RANDOMIZEDBATCHTEST Test over may flow examples of a type of cost
%function, with some added randomized trials

definedOrDefault = @(name,default) definedOrDefault_long(name,default,varargin);
% Define parameters
ids = definedOrDefault('ids',2:1199);
flowHandle_default = @(i) ImageFlow('image_index',i,'MaxBidirectionalError',0.5);
flowHandle = definedOrDefault('flowHandle', flowHandle_default);
rngSeed = definedOrDefault('rngSeed', -1); 
nsamples = definedOrDefault('nsamples', 25);
ntrials = definedOrDefault('ntrials', 100);
ntrialsamples = definedOrDefault('ntrialsamples', 800);
verbose = definedOrDefault('verbose', false);

% Set randomness
if rngSeed > 0
    rng(rngSeed)
end

% Get positions
c = CostFunctionFactory(type,flowHandle(ids(1)));
[~,x,y] = c.getSurfaceResiduals(nsamples);

% Initialize results array
results(1,1) = c.getResults;
results(length(ids),ntrials) = results(1,1);

% Store the flows
flows = cell(length(ids),ntrials);

% Iterate over all the flows
for i = 1:length(ids)
    fprintf('Image %d of %d (value %d)\n',i,length(ids),ids(i));
    flow = flowHandle(ids(i));
    t0 = CTimeleft(ntrials);
    for k = 1:ntrials
        t0.timeleft();
        % Randomize the flow
        flow_sample = flow;
        samples = randperm(flow.nPoints,ntrialsamples);
        flow_sample.xy = flow_sample.xy(:,samples);
        flow_sample.uv = flow_sample.uv(:,samples);
        flow_sample.xy_pixel = flow_sample.xy_pixel(:,samples);
        flow_sample.uv_pixel = flow_sample.uv_pixel(:,samples);
        flow_sample.nPoints = ntrialsamples;
        flows{i,k} = flow_sample;
        % Build the cost function
        c = CostFunctionFactory(type,flow_sample,verbose);
        % Store results
        results(i,k) = c.getResults(nsamples);
    end
end


end

