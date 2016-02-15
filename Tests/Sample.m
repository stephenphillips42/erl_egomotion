% Example of how the classes typically work
Im1 = imread('Data/Sample/000014.png');
Im2 = imread('Data/Sample/000015.png');
ld = load('Data/Sample/samplepose.mat');

% Create Flow Objects
flow = ImageFlow(Im1,Im2,'K',ld.K,'trueT',ld.trueT,'trueOmega',ld.trueOmega);

% Create Cost Function Object
c = CostFunctionFactory('RobustERL',flow);

% Plot everything
flow.plotFlow
c.plotResidualsSurface;


