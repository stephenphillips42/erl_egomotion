classdef ImageFlow < Flow
    %IMAGEFLOW Flow based on images in KITTI dataset
    %   Flow extracted from the image sequence in the KITTI dataset, the location specified by
    %   baseDir. What flow you get is complete specified by the sequence and the image index
    
    properties % (SetAccess = private)
        img0
        img1
    end
    
    methods
        function imgFlow = ImageFlow(img0,img1,varargin)
            imgFlow@Flow
            definedOrDefault = @(name,default) ...
                                 definedOrDefault_long(name,default,varargin);
            imgFlow.img0 = img0;
            imgFlow.img1 = img1;

            % Extract the KLT flow from images
            minQuality = definedOrDefault('MinQuality',0.01);
            nPyrLevels = definedOrDefault('NumPyramidLevels',4);
            maxBiError = definedOrDefault('MaxBidirectionalError',1);
            points = detectMinEigenFeatures(img0,'MinQuality',minQuality);
            pointTracker = vision.PointTracker('NumPyramidLevels',nPyrLevels,...
                                               'MaxBidirectionalError',maxBiError);
            initialize(pointTracker, points.Location,img0);
            [new_points,valid] = step(pointTracker, img1);
            % Set up the points in calibrated coordinates
            new_points = double(new_points(valid,:));
            points = double(points.Location(valid,:));

            % Calibrate the points for Cperp computation and residual computation
            imgFlow.xy = [1,0,0;
                          0,1,0]*(imgFlow.K\[ points';
                                              ones(1,size(points,1)) ]);
            xy2 = [1,0,0;
                   0,1,0]*(imgFlow.K\[ new_points';
                                       ones(1,size(new_points,1)) ]);
            % Get the optical flow
            imgFlow.uv = xy2 - imgFlow.xy;
            % Get pixel values
            imgFlow.xy_pixel = points';
            imgFlow.uv_pixel = imgFlow.K(1,1)*imgFlow.uv;
            
            % Get the number of flow value we have
            imgFlow.nPoints = length(points);
            
            % Get true pose
            imgFlow.trueT = definedOrDefault('trueT',[0;0;0]);
            imgFlow.trueOmega = definedOrDefault('trueOmega',[0;0;0]);
        end
        function plotFlow(imgFlow,isOutlier,newFigure)
            if nargin < 2
                isOutlier = [];
            end
            if nargin < 3
                newFigure = true;
            end
            if newFigure; figure; end;
            if ~isempty(imgFlow.img0)
                imshow(imgFlow.img0)
            end
            hold on;
            plotFlow@Flow(imgFlow,isOutlier,false,false);
            hold off;
        end
    end
    
end

