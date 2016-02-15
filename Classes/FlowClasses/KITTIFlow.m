classdef KITTIFlow < Flow
    %IMAGEFLOW Flow based on images in KITTI dataset
    %   Flow extracted from the image sequence in the KITTI dataset, the location specified by
    %   baseDir. What flow you get is complete specified by the sequence and the image index
    
    properties % (SetAccess = private)
        sequence_index = 10 % TODO: Make these zero indexed
        image_index = 15 % TODO: Make these zero indexed
        img = []
        baseDir
        imgDir
        poseDir
    end
    
    methods
        function flow = KITTIFlow(baseDir, varargin)
            flow@Flow
            definedOrDefault = @(name,default) ...
                                 definedOrDefault_long(name,default,varargin);
            
            flow.sequence_index = definedOrDefault('sequence_index',10);
            flow.image_index = definedOrDefault('image_index',15);
            flow.baseDir = baseDir;
            flow.poseDir = definedOrDefault('poseDir',...
                                    fullfile(flow.baseDir,'poses'));
            flow.imgDir = definedOrDefault('imgDir',...
                                fullfile(flow.baseDir,'sequences',...
                                         sprintf('%02d',flow.sequence_index),'image_0'));
            imgFile0 = ...
                fullfile(flow.imgDir,...
                            sprintf('%06d.png',flow.image_index-1));
            imgFile1 = ...
                fullfile(flow.imgDir,...
                            sprintf('%06d.png',flow.image_index));
            img0 = imread(imgFile0);
            img1 = imread(imgFile1);

            flow.img = img0;
            % Extract the KLT flow from images
            minQuality = definedOrDefault('MinQuality',0.01);
            nPyrLevels = definedOrDefault('NumPyramidLevels',4);
            maxBiError = definedOrDefault('MaxBidirectionalError',0.01);
            points = detectMinEigenFeatures(img0,'MinQuality',minQuality);
            pointTracker = vision.PointTracker('NumPyramidLevels',nPyrLevels,...
                                               'MaxBidirectionalError',maxBiError);
            initialize(pointTracker, points.Location,img0);
            [new_points,valid] = step(pointTracker, img1);
            % Set up the points in calibrated coordinates
            new_points = double(new_points(valid,:));
            points = double(points.Location(valid,:));

            % Calibrate the points for Cperp computation and residual computation
            flow.xy = [1,0,0;
                          0,1,0]*(flow.K\[ points';
                                              ones(1,size(points,1)) ]);
            xy2 = [1,0,0;
                   0,1,0]*(flow.K\[ new_points';
                                       ones(1,size(new_points,1)) ]);
            % Get the optical flow
            flow.uv = xy2 - flow.xy;
            % Get pixel values
            flow.xy_pixel = points';
            flow.uv_pixel = flow.K(1,1)*flow.uv;
            
            % Get the number of flow value we have
            flow.nPoints = length(points);
            
            % Get true pose
            [ flow.trueT, flow.trueOmega ] = ...
                flow.extractPose();
        end
        function plotFlow(imgFlow,isOutlier,newFigure)
            if nargin < 2
                isOutlier = [];
            end
            if nargin < 3
                newFigure = true;
            end
            if newFigure; figure; end;
            if ~isempty(imgFlow.img)
                imshow(imgFlow.img)
            end
            hold on;
            plotFlow@Flow(imgFlow,isOutlier,false,false);
            hold off;
        end
        
        function [ trueT, trueOmega ] = extractPose( imgFlow )
            if ~exist(fullfile(imgFlow.poseDir,'poses.mat'),'file')
                imgFlow.createPosesFile;
            end
            load(fullfile(imgFlow.poseDir,'poses.mat'))
            
            pose = ...
                reshape(...
                   P{imgFlow.sequence_index+1}(imgFlow.image_index-1,:),...
                   4,3)';
            pose_next = ...
                reshape(...
                   P{imgFlow.sequence_index+1}(imgFlow.image_index,:),...
                   4,3)';
            if imgFlow.image_index == 1
                trueT = zeros(3,1);
                trueOmega = zeros(3,1);
            else
                R_change = pose(:,1:3)'*pose_next(:,1:3);
                T_change = pose(:,1:3)'*(pose_next(:,4) - pose(:,4));
                trueT = T_change/norm(T_change);
                omega = (R_change-R_change')/2;
                trueOmega = [ omega(3,2);
                              omega(1,3);
                              omega(2,1) ]; % Extracting skew symmetric
            end
            
        end
        
        function createPosesFile(imgFlow)
            P = cell(11,1);
            for i = 1:11
                P{i} = imgFlow.readInPoses(i-1);
            end
            save(fullfile(imgFlow.poseDir, 'poses.mat'),'P')
        end
        
        function posesfull = readInPoses(imgFlow,seq_index)
            % Read files in from KITTI Odometry data format
            % Get the file to read in the poses
            pose_file = ...
                fullfile(imgFlow.poseDir,sprintf('%02d.txt',seq_index));

            pid = fopen(pose_file);

            % Get the info from the file
            posesinfo = textscan(pid, ...
                '%f %f %f %f %f %f %f %f %f %f %f %f', ...
                'delimiter', '\n');
            
            % Format pose_info into what we use
            posesfull = zeros(length(posesinfo{1}),12);
            for i = 1:length(posesinfo{1})
                posesfull(i,:) = ...
                    [ posesinfo{1}(i) posesinfo{2}(i) posesinfo{3}(i) ...
                      posesinfo{4}(i) posesinfo{5}(i) posesinfo{6}(i) ...
                      posesinfo{7}(i) posesinfo{8}(i) posesinfo{9}(i) ...
                      posesinfo{10}(i) posesinfo{11}(i) posesinfo{12}(i) ];
            end
            
            % Clean up the file
            fclose(pid);
        end
    end
    
end

