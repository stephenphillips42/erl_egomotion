function RunBatchComputationTmp(types,seqids,myids)

% types = {'Algebraic','ZhangMatrix','Epipolar','RobustERRL','RobustLifted'};
% types = {'Algebraic'};
% ids = [ 144 258 301 390 511 653 743 888 913 1163 1171 1200 ];
load('/scratch/KITTIOdometry/matlabPoses/poses.mat'); % Load P
for i = 1:length(seqids)
        seqid = seqids(i);
        flowHandle = @(i) ...
            ImageFlow('image_index',i,'sequence_index',seqid,...
                             'width',370,'height',370,'MaxBidirectionalError',1);
        % fprintf('Sequence Number %d\n',seqid);
        if nargin == 3
                ids = myids{i};
        else
                ids = 2:(length(P{seqid+1})-1);
        end
        [ estT, estOmega, trueT, trueOmega ] = ...
            LightMultiBatchTest( types,'flowHandle',flowHandle,'ids',ids );
        save(sprintf('~/Outputs/results_epi_seq%02d_rng%d_%d.mat', seqid, min(ids), max(ids)),...
             'estT', 'estOmega', 'trueT', 'trueOmega')
end

end

