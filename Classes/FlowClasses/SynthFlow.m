classdef SynthFlow < Flow
    %SYNTHFLOW Synthetically generated flow.
    %   The SynthFlow class generates idea flow vectors, plus noise specified by the noiseLevel and
    %   outliers specified by the outlierRate, to generate syntetic flow. A seed can be specifed to
    %   reporoduce results, the parameter rng_seed. Also how many points are used can be specified
    %   by nPoints.
    
    properties (SetAccess = private)
        % Parameters needed to generate the flow
        % Random number parameters
        rngSeed % RNG seed to be able to replicate experiments
        % Noise parameters
        noiseLevel % Noise in standard deviation units
        outlierRate % Fraction outliers
        outliers % Which ones are outliers
        % Scaling parameters
        lambda % Scale to make change more visible (order of the focal length)
        kappa % Scale for rotation
        % Sampling parameters
        depthMin % Minimum depth expected
        depthMax % Maximum depth expected
        % Other flow parameters
        depths % All the detphs of all the points xy
        scaledT % Translation with scale
        uv_noiseless % The 'ideal' flow
    end
    
    methods
        % TODO: Make complicated function where all parameters are
        % optional... this is getting too confusing...
        function synthFlow = SynthFlow(varargin)
            synthFlow@Flow
            definedOrDefault = @(name,default) ...
                                 definedOrDefault_long(name,default,varargin);
            
            % Set up random number generator -- for repeatability
            synthFlow.rngSeed = definedOrDefault('rngSeed',-1);
            if synthFlow.rngSeed > 0
                rng(synthFlow.rngSeed);
            end
            
            % Set up sample points
            synthFlow.nPoints = definedOrDefault('nPoints',300);
            x = randi(synthFlow.imgsize(2),[synthFlow.nPoints,1]);
            y = randi(synthFlow.imgsize(1),[synthFlow.nPoints,1]);
            synthFlow.xy_pixel = [x y]'; % Uncalibrated pixels
            xy = ([ 1 0 0;
                    0 1 0 ] * (synthFlow.K\[ x(:) y(:) ones(size(x))]'));
            synthFlow.xy = xy; % Calibrated pixels

            % Create depths
            synthFlow.depthMin = definedOrDefault('depthMin',5);
            synthFlow.depthMax = definedOrDefault('depthMax',15);
            synthFlow.depths = ...
                (synthFlow.depthMax - synthFlow.depthMin)*rand(synthFlow.nPoints,1) ...
                    + synthFlow.depthMin;

            % Generate motions
            synthFlow.lambda = definedOrDefault('lambda',1);
            synthFlow.kappa = definedOrDefault('kappa',0.2);
            % Translation
            synthFlow.scaledT = definedOrDefault('scaledT',synthFlow.lambda*randn([3 1]));
            synthFlow.trueT = normc(synthFlow.scaledT);
            if synthFlow.trueT(3) < 0; synthFlow.trueT = -synthFlow.trueT; end;
            % Rotation
            synthFlow.trueOmega = ...
                definedOrDefault('trueOmega',synthFlow.kappa*randn([3 1]));

            % Generate optical flow
            A = @(XY) [ -1  0 XY(1);
                         0 -1 XY(2) ];
            B = @(XY) [ XY(1)*XY(2), -(1 + XY(1)^2),  XY(2);
                        1 + XY(2)^2,   -XY(1)*XY(2), -XY(1) ];
            
            uv = zeros(2,synthFlow.nPoints);
            for i = 1:length(uv)
               uv(:,i) = ...
                   (1./synthFlow.depths(i))*A(xy(:,i))*synthFlow.scaledT + ...
                   B(xy(:,i))*synthFlow.trueOmega;
            end
            synthFlow.uv_noiseless = uv;

            % Generate the flow noise
            synthFlow.noiseLevel = definedOrDefault('noiseLevel',0.1);
            sgm = std(uv(:));
            uv = uv + randn(size(uv))*(sgm*synthFlow.noiseLevel);
            
            % Generate the flow outliers
            synthFlow.outlierRate = definedOrDefault('outlierRate',0);
            synthFlow.outliers = false(synthFlow.nPoints,1);
            % We redraw outlier flow from the current magnitude
            % distribution, but in random directions
            if synthFlow.outlierRate > 0
                % Choose the outliers
                nOutliers = floor(synthFlow.outlierRate*(synthFlow.nPoints));
                outlierInds = randperm(synthFlow.nPoints,nOutliers);
                synthFlow.outliers(outlierInds) = true;
                
                % Generate outlier vectors
                magnitudes = sqrt(sum(uv.^2));
                outlier_mags = abs(mean(magnitudes) + std(magnitudes)*randn(1,nOutliers));
                outlier_theta = 2*pi*rand(1,nOutliers);
                uv(:,outlierInds) = [ cos(outlier_theta).*outlier_mags;
                                      sin(outlier_theta).*outlier_mags ];
            end
            
            % Finally set the flow values
            synthFlow.uv = uv;
            synthFlow.uv_pixel = synthFlow.K(1,1)*uv;

        end
        
        function plotFlow(synthFlow,plotOutliers,newFigure)
            if nargin < 2
                plotOutliers = true;
            end
            if nargin < 3
                newFigure = true;
            end
            if plotOutliers
                plotFlow@Flow(synthFlow,synthFlow.outliers,newFigure,true);
            else
                plotFlow@Flow(synthFlow,[],newFigure,false);
            end
            hold on;
            axis equal;
            set(gca,'ydir','reverse')
            ylim([0 synthFlow.imgsize(1)])
            xlim([0 synthFlow.imgsize(2)])
            hold off
        end
        function plotFlowOutliers(synthFlow,outlierWeights,selectedOutliers,newFigure)
            if nargin < 3
                selectedOutliers = ...
                    outlierWeights < mean(outlierWeights) - 2*std(outlierWeights);
            end
            if nargin < 4
                newFigure = true;
            end
            % Begin plot
            if newFigure; figure; end;
            % Plot points (blue)
            hold on
            colormap('summer')
            scatter(synthFlow.xy_pixel(1,:),synthFlow.xy_pixel(2,:),300,outlierWeights,'.')
            plot(synthFlow.xy_pixel(1,selectedOutliers),...
                 synthFlow.xy_pixel(2,selectedOutliers),...
                 'ko','markersize',30)
            plot(synthFlow.xy_pixel(1,synthFlow.outliers),...
                 synthFlow.xy_pixel(2,synthFlow.outliers),...
                 'ro','markersize',25)
             
            % Plot the arrows to make it look pretty
            quiver(...
                synthFlow.xy_pixel(1,:),...
                synthFlow.xy_pixel(2,:),...
                synthFlow.uv_pixel(1,:),...
                synthFlow.uv_pixel(2,:),...
                'AutoScale','off')
            axis equal;
            set(gca,'ydir','reverse')
            ylim([0 synthFlow.imgsize(1)])
            xlim([0 synthFlow.imgsize(2)])
            legend('Flow positions','Selected Outliers','True Outliers')

            hold off
        end
    end
        
end

