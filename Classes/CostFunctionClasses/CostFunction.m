classdef CostFunction
    %RESIDUAL Superclass for the residual calculations
    %   The basic functions for getting the residual from a given method. All a subclass really
    %   need to implement is the getResidual function and the constructor 
    
    properties
        flow % Optical flow we are using the residual on
        trueT % Must know to do actual error analysis
        % For outlier rejection
        flow_outliers % Flow before outlier rejection (only used if outliers removed)
        inliers % Which flow vectors are inliers?
        verbose % Print out progress
    end
    
    methods
        % Constructor - all precomputation done here
        function c = CostFunction(flow)
            c.flow = flow;
            c.flow_outliers = flow;
            c.trueT = flow.trueT;
            c.verbose = true;
        end
        % Residual of a translation
        function [resid] = getResidual(c,T)
            resid = sum(c.getFlowResiduals(T));
        end
        % Residual for each flow vector
        function [flowResids] = getFlowResiduals(c,T)
            flowResids = zeros(c.flow.nPoints,1);
        end
        function [Omega] = getOmega(c,T)
            Omega = zeros(size(T));
        end
        function [rho] = getInverseDepths(c,T,Omega)
            rho = zeros(flow.nPoints,1);
        end
        % Get residuals for each translation on the sphere
        function [residuals, translations] = getSphereResiduals(c,sphereDensity)
            if nargin < 2
                sphereDensity = 14;
            end
            [x,y,z] = sphere(sphereDensity);
            xp = x(z>=0);
            yp = y(z>=0);
            zp = z(z>=0);
            translations = unique([xp(:) yp(:) zp(:)],'rows')';

            % Now for each translation compute the Cperp
            residuals = zeros(1,length(translations));
            t0 = CTimeleft(length(translations));
            for t = 1:length(translations)
                if c.verbose; t0.timeleft(); end;
                residuals(t) = c.getResidual(translations(:,t));
            end
        end
        function [z,x,y] = getSurfaceResiduals(c,nsamples)
            if nargin < 2
                nsamples = 25;
            end
            theta = linspace(0,2*pi,nsamples);
            r = linspace(0,1,nsamples);
            [TH,R] = meshgrid(theta,r);
            Z = zeros(size(TH)); % z=(x^2)-(y^2)
            [x,y,z] = pol2cart(TH,R,Z);

            % Compute the residuals
            t0 = CTimeleft(length(x(:)));
            for i = 1:length(x(:))
                if c.verbose; t0.timeleft(); end;
                zval = sqrt(1 - (x(i)^2 + y(i)^2));
                if ~isreal(zval); zval = 0; end;
                T = [x(i); y(i); zval];
                z(i) = c.getResidual(T);
            end
        end
        
        % Go to local minima starting at initPos
        function [finalPos] = gradientDescent(c,initPos)
            cost = @(T) c.getResidual(T);
            options = optimset('Display','off');
            ceq = @(T) T'*T - 1;
            cineq = @(T) -1;
            finalPos = fmincon(cost, initPos, [],[],[],[],[],[], ...
                                @(T) deal(cineq(T),ceq(T)),options);
        end

        % Give score of flow vectors to remove as outliers - using avg.
        % likelihood
        function [O] = outlierRejection(c)
            nsamples = 11; % Parameters to choose
            theta = linspace(0,2*pi,nsamples);
            r = linspace(0,1,nsamples);
            [TH,R] = meshgrid(theta,r);
            [x,y] = pol2cart(TH,R);
            translations = unique([x(:) y(:) sqrt(abs(1 - (x(:).^2 + y(:).^2)))],'rows')';
            % Build error matrix
            E = zeros(c.flow.nPoints,length(translations));
            for j = 1:length(translations)
                E(:,j) = c.getFlowResiduals(translations(:,j));
            end
            % Renormalize
            E = diag(sqrt(sum(c.flow.uv.^2)))\E;
            % Get mean vector
            M = mean(E)';
            O = sum(exp(-E/diag(M))/diag(M),2);
        end
        
        % Give score of flow vectors to remove as outliers - using Joint.
        % log likelihood
        function [O] = outlierRejectionJoint(c)
            nsamples = 11; % Parameters to choose
            theta = linspace(0,2*pi,nsamples);
            r = linspace(0,1,nsamples);
            [TH,R] = meshgrid(theta,r);
            [x,y] = pol2cart(TH,R);
            translations = unique([x(:) y(:) sqrt(abs(1 - (x(:).^2 + y(:).^2)))],'rows')';
            % Build error matrix
            E = zeros(c.flow.nPoints,length(translations));
            for j = 1:length(translations)
                E(:,j) = c.getFlowResiduals(translations(:,j));
            end
            % Renormalize
            E = spdiags(sqrt(sum(c.flow.uv.^2)),0,c.flow.nPoints,c.flow.nPoints)\E;
            % Get mean vector
            M = mean(E)';
            O = zeros(c.flow.nPoints,1);
            for j = 1:length(translations)
                O = O - E(:,j)/M(j) + log(M(j));
            end
        end
        
        % For brevity in batch computations
        function [result] = getResults(c,nsamples)
            if nargin < 2; nsamples = 25; end;
            [z, x, y] = c.getSurfaceResiduals(nsamples); 
            [~, mi] = min(z(:));
            initT = [x(mi); y(mi); sqrt(max(0,1 - (x(mi)^2 + y(mi)^2)))];
            finalT = c.gradientDescent(initT);
            angError = acosd(sign(finalT'*c.trueT)*min(abs(finalT'*c.trueT),1));
            
            result.trueT = c.trueT;
            result.angularError = angError;
            result.minimumT = finalT;
            result.minimumOmega = c.getOmega(finalT);
            result.minimumResidual = c.getResidual(finalT);
            result.extraData = {};
        end
        
        % Different plotting functions
        function plotResiduals(c,useGradientDescent,sphereDensity)
            if nargin < 2; useGradientDescent = false; end
            if nargin < 3; sphereDensity = 14; end
            
            % Get our estimate
            [sphereResiduals, translations] = c.getSphereResiduals(sphereDensity);
            [~,t_min_orig] = min(sphereResiduals);
            gridResidualsMinT = translations(:,t_min_orig);

            % Create the plot
            figure;
            hold on;
            scatter3( ...
                translations(1,:),...
                translations(2,:),...
                translations(3,:),...
                500,sphereResiduals,'.');

            colormap('parula')
            axis equal;

            % Plot the true translation
            scatter3(c.trueT(1),c.trueT(2),c.trueT(3),500,...
                c.getResidual(c.trueT),'.');
            scatter3(c.trueT(1),c.trueT(2),c.trueT(3),400,'go')

            % Plot minimum value for computed residuals
            scatter3(gridResidualsMinT(1),gridResidualsMinT(2),gridResidualsMinT(3),500,'ro')
            % Use gradient descent
            if useGradientDescent
                finalT = c.gradientDescent(gridResidualsMinT);
                scatter3(finalT(1),finalT(2),finalT(3), 500, ...
                            c.getResidual(finalT), '.');
                scatter3(finalT(1),finalT(2),finalT(3),500,'co');
            end
            
            colorbar
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            view(2)
        end
        function [min_value] = plotResidualsSurface(c,newFigure,useGradientDescent,nsamples)
            if nargin < 2; newFigure = true; end
            if nargin < 3; useGradientDescent = false; end
            if nargin < 4; nsamples = 25; end
            
            [z,x,y] = c.getSurfaceResiduals(nsamples);
            % Plot the residuals
            if newFigure; figure('units','pixels','position',[0 0 1001 1001]); end;
            hold on;
            surf(x,y,z)
            colormap('jet')

            alpha(0.4)
            % Get our estimate
            [~,t_min_guess] = min(z(:));

            scatter3(c.trueT(1),c.trueT(2),...
                c.getResidual(c.trueT),600,'g.');
            guessedT = [ x(t_min_guess);
                         y(t_min_guess);
                         sqrt(1 - (x(t_min_guess)^2 + y(t_min_guess)^2))];
            scatter3(x(t_min_guess),y(t_min_guess),...
                c.getResidual(guessedT),600,'k.');
            % Use gradient descent
            if useGradientDescent
                finalT = c.gradientDescent(guessedT);
                scatter3(finalT(1),finalT(2),...
                    c.getResidual(finalT),500,'c.');
            end
            xlabel('X')
            ylabel('Y')
            zlabel('Residual')
            hold off
            if nargout > 0
                min_value = guessedT;
            end
        end
        function plotResidualsHeatmap(c,newFigure,useGradientDescent,nsamples,nlevels)
            if nargin < 2
                newFigure = true;
            end
            if nargin < 3
                useGradientDescent = false;
            end
            if nargin < 4
                nsamples = 25;
            end
            if nargin < 5
                nlevels = 30;
            end
            
            [z,x,y] = c.getSurfaceResiduals(nsamples);
            % Plot the residuals
            if newFigure; figure; end;
            hold on;
            contourf(x,y,z,nlevels)
            colormap('jet')

            alpha(0.4)
            % Get our estimate
            [~,t_min_guess] = min(z(:));

            scatter(c.trueT(1),c.trueT(2),600,'g.');
            guessedT = [x(t_min_guess);
                             y(t_min_guess);
                             sqrt(1 - (x(t_min_guess)^2 + y(t_min_guess)^2))];
            scatter(x(t_min_guess),y(t_min_guess),600,'k.');
            % Use gradient descent
            if useGradientDescent
                finalT = c.gradientDescent(guessedT);
                scatter3(finalT(1),finalT(2),500,'c.');
            end
            xlabel('X')
            ylabel('Y')
            zlabel('Residual')
            axis equal
            hold off
        end
    end
end

