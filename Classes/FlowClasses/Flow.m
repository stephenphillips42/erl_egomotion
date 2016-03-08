classdef Flow
%Flow - Superclass that extracts flow
%   The FLOW class extracts the optical flow uv (a 2 x N array) and the flow's positions 
%   (a 2 x N array). It also gives the pixel coordinates of said flow, uv_pixel and
%   xy_pixel, along with the associated calibration matrix K and image size imgsize. The
%   plotFlow function is for convenience to visualize the flow
    
    properties
        xy = zeros(2,1);       % Calibrated flow locations
        uv = zeros(2,1);       % Calibrated flow
        xy_pixel = zeros(2,1); % Uncalibrated flow locations
        uv_pixel = zeros(2,1); % Uncalibrated flow
        nPoints                % Number of flow points

        % For testing purposes
        trueT = [0; 0; 0]; % Totally wrong defaults of zero
        trueOmega = [0; 0; 0];

        % Parameters
        K = [ 707.0912         0  601.8873;
                     0  707.0912  183.1104;
                     0         0    1.0000 ];
        imgsize = [370 1226];
    end
    
    methods
        function plotFlow(flow,isOutlier,newFigure,scaleFlow)
            % Begin plot
            if newFigure; figure; end;
            % Get the flow points
            new_points = flow.xy_pixel + flow.uv_pixel;
            % Plot points (blue) and new points (green)
            hold on
            scatter(flow.xy_pixel(1,:),flow.xy_pixel(2,:),'b.')
            if ~isempty(isOutlier) % Circle the outliers of the points
                plot(flow.xy_pixel(1,isOutlier),flow.xy_pixel(2,isOutlier),'ko','markersize',30)
            end
            scatter(new_points(1,:),new_points(2,:),'g.')
            % Plot the arrows to make it look pretty
            if scaleFlow % For debugging
                quiver(...
                    flow.xy_pixel(1,:),...
                    flow.xy_pixel(2,:),...
                    flow.uv_pixel(1,:),...
                    flow.uv_pixel(2,:))
            else
                quiver(...
                    flow.xy_pixel(1,:),...
                    flow.xy_pixel(2,:),...
                    flow.uv_pixel(1,:),...
                    flow.uv_pixel(2,:),0)
            end

            hold off
        end
        
    end

end

