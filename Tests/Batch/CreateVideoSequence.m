
seqno = 10;
inds_all = 350:460;
inds_multi = [2 3 4 5 6 7 8];
inds_single = 2;
params = createParams();

% profile on
for i = inds_all
    trueT = getPose(params,i);
    flow = ImageFlow(seqno,i);
    residual = AlgebraicResidual(flow,trueT);
    figure
    subplot(2,1,1)
    imshow(flow.img)
    subplot(2,1,2)
    flow.plotFlow([],false)
    hold on
    foe = flow.K*(trueT/trueT(3));
    scatter(foe(1),foe(2),700,'r.')
    hold off
    print(sprintf('/scratch/KITTIOdometry/2015-06-24/FlowAndImage_%04d.png',i),'-dpng')
    close all
    residual.plotResidualsSurface;
    print(sprintf('/scratch/KITTIOdometry/2015-06-24/ResidualType1_%04d.png',i),'-dpng')
    close all
end
profile off
profile viewer


