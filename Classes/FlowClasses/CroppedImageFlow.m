classdef CroppedImageFlow < Flow
    %IMAGEFLOW Flow based on images in KITTI dataset
    %   Flow extracted from the image sequence in the KITTI dataset, the location specified by
    %   baseDir. What flow you get is complete specified by the sequence and the image index
    
    properties % (SetAccess = private)
        sequence_index = 10 % TODO: Make these zero indexed
        image_index = 15 % TODO: Make these zero indexed
        img0
        img1
        baseDir
        imgDir
        poseDir
    end
    
    methods
        function imgFlow = CroppedImageFlow(varargin)
            imgFlow@Flow
            definedOrDefault = @(name,default) ...
                                 definedOrDefault_long(name,default,varargin);
            
            imgFlow.sequence_index = definedOrDefault('sequence_index',10);
            imgFlow.image_index = definedOrDefault('image_index',15);
            imgFlow.baseDir = definedOrDefault('baseDir','/scratch/KITTIOdometry/');
            imgFlow.poseDir = definedOrDefault('poseDir',...
                                    fullfile(imgFlow.baseDir,'matlabPoses'));
            imgFlow.imgDir = definedOrDefault('imgDir',...
                                fullfile(imgFlow.baseDir,'dataset','sequences',...
                                         sprintf('%02d',imgFlow.sequence_index),'image_0'));
            imgFile0 = ...
                fullfile(imgFlow.imgDir,...
                            sprintf('%06d.png',imgFlow.image_index-1));
            imgFile1 = ...
                fullfile(imgFlow.imgDir,...
                            sprintf('%06d.png',imgFlow.image_index));
            imgFlow.img0 = imread(imgFile0);
            imgFlow.img1 = imread(imgFile1);

            % Crop images
            % imgFlow.imgsize = [370 1226];
            height = definedOrDefault('height',370);
            width = definedOrDefault('width',370);
            imgcropsize = [(imgFlow.imgsize(2)-width)/2 ...
                           (imgFlow.imgsize(1)-height)/2 ...
                           width height ];
            imgFlow.img0 = imcrop(imgFlow.img0,imgcropsize);
            imgFlow.img1 = imcrop(imgFlow.img1,imgcropsize);
            imgFlow.imgsize = [ height width ];
            
            % Extract the KLT flow from images
            minQuality = definedOrDefault('MinQuality',0.01);
            nPyrLevels = definedOrDefault('NumPyramidLevels',4);
            maxBiError = definedOrDefault('MaxBidirectionalError',0.01);
            points = detectMinEigenFeatures(imgFlow.img0,'MinQuality',minQuality);
            pointTracker = vision.PointTracker('NumPyramidLevels',nPyrLevels,...
                                               'MaxBidirectionalError',maxBiError);
            initialize(pointTracker, points.Location,imgFlow.img0);
            [new_points,valid] = step(pointTracker, imgFlow.img1);
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
            % poses = imgFlow.readInPoses;
            % [ imgFlow.trueT, imgFlow.trueOmega ] = ...
            %     imgFlow.getPose(poses,imgFlow.image_index);
            [ imgFlow.trueT, imgFlow.trueOmega ] = ...
                imgFlow.extractPose();
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
%     end
%     
%     methods (Access = private)
        
        function [ trueT, trueOmega ] = extractPose( c )
            load(fullfile(c.poseDir,'poses.mat'))
            
            pose = reshape(P{c.sequence_index+1}(c.image_index-1,:),4,3)';
            pose_next = reshape(P{c.sequence_index+1}(c.image_index,:),4,3)';
            if c.image_index == 1
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
        
        function [ trueT, trueOmega, pose ] = getPose( ~, poses, i )
            %GETPOSE Loading pose i based on the param struct, and the number i
            pose = poses{i};

            if i == 1
                trueT = zeros(3,1);
                trueOmega = zeros(3,1);
            else
                R_change = poses{i-1}(:,1:3)'*poses{i}(:,1:3);
                T_change = poses{i-1}(:,1:3)'*(poses{i}(:,4) - poses{i-1}(:,4));
                trueT = T_change/norm(T_change);
                omega = (R_change-R_change')/2;
                % TODO: Verify this is right
                trueOmega = [ omega(3,2);
                              omega(1,3);
                              omega(2,1) ]; % Extracting skew symmetric
            end
        end
        
        function posesfull = readInPoses(imgFlow)
            % Get the file to read in the poses
            pose_file = ...
                fullfile(imgFlow.poseDir,sprintf('%02d.txt',imgFlow.sequence_index));

            % WARNING: This does NOT work on Windows
            [~, numLinesStr] = system(['wc -l ' pose_file]);
            nposes = str2num(strtok(numLinesStr));

            % nposes = 1201; % Number of poses in a file
            % sequence_id = 10;
            posesfull = cell(nposes,1);
            pid = fopen(pose_file);

            % Get the first line
            pose_info = textscan(pid, '%f %f %f %f %f %f %f %f %f %f %f %f', 1, 'delimiter', '\n');
            posesfull{1} = reshape(cell2mat(pose_info),4,3)';
            % Read in every line
            for i = 2:nposes %1:2:length(poses_info)

                % Read out this line and the following one to get the R,T
                pose_info = textscan(pid, '%f %f %f %f %f %f %f %f %f %f %f %f', 1);
                posesfull{i} = reshape(cell2mat(pose_info),4,3)';
            end
            
            assert(numel(posesfull) == nposes)

            % Clean up the file
            fclose(pid);
        end
    end
    
end

