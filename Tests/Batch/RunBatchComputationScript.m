
%%
% types = {'Algebraic','ZhangMatrix','Epipolar','RobustERRL','RobustLifted'};
types = {'RobustLifted'};
% ids = [ 144 258 301 390 511 653 743 888 913 1163 1171 1200 ];
load('/scratch/KITTIOdometry/matlabPoses/poses.mat'); % Load P
% for seqid = 0:10
seqid = 0;
profile on
    fprintf('Sequence Number %d\n',seqid);
%     ids = 2:(length(P{seqid+1})-1);
    ids = 2:46;
    flowHandle = @(i) ImageFlow('image_index',i,'sequence_index',seqid,'MaxBidirectionalError',1);
    [ estT, estOmega, trueT, trueOmega ] = ...
        LightMultiBatchTest( types,'flowHandle',flowHandle,'ids',ids );
profile of
profile view
%     save(sprintf('Outputs/2015-09-11/results_opt_seq%02d.mat',seqid), 'estT', 'estOmega', 'trueT', 'trueOmega')
% end


%%
for k = 1:size(angErrs,1)
    ecdf(angErrs(k,:)); 
    hold on
end

hold off
legend(types)
title('Translational Velocity Error')
xlabel('Error (deg)')
ylabel('Probability of Error')


%%
load(sprintf('Outputs/2015-09-11/results_seq%02d.mat',0), 'estT', 'estOmega', 'trueT', 'trueOmega')
angErrs = zeros(length(types),length(trueT));
for k = 1:length(types)
    angErrs(k,:) = acosd(trueT(1,:).*squeeze(estT(k,:,1)) + ...
                            trueT(2,:).*squeeze(estT(k,:,2)) + ...
                            trueT(3,:).*squeeze(estT(k,:,3)));
    ecdf(angErrs(k,:)); 
    hold on
end
hold off
legend(types)
title('Translational Velocity Error')
xlabel('Error (deg)')
ylabel('Probability of Error')

%%

omegaErrs = zeros(length(types),length(trueOmega));
for k = 1:length(types)
    omegaErrs(k,:) = sqrt((trueOmega(1,:) - squeeze(estOmega(k,:,1))).^2 + ...
                          (trueOmega(2,:) - squeeze(estOmega(k,:,2))).^2 + ...
                          (trueOmega(3,:) - squeeze(estOmega(k,:,3))).^2) ./ ...
                     sqrt(trueOmega(1,:).^2 + trueOmega(2,:).^2 + trueOmega(3,:).^2);
    ecdf(100.*omegaErrs(k,:));
    hold on
end
hold off
legend(types)
title('Translational Velocity Error')
xlabel('Percent Error')
ylabel('Probability of Error')



%% 

for k = 1:length(types)
    plot(velNorms(inds), angErrs(k,inds)); 
    hold on
end
legend(types)
hold off

