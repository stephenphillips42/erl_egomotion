% Air Ground Dataset
basdir = '/home/stephen/Documents/Research/data/air-ground-data';
gt = load(fullfile(basdir,'ground_truth.mat'));

% for i = 1:404
%     Im1 = imread(sprintf('left-%04d00',i));
%     Im2 = imread(sprintf('left-%04d00',i+1));
%     flow = ImageFlow(Im1,Im2,)
    
% end

% Kagaru Dataset
K1 = [ 1641.99751,          0, 642.15139;
                0, 1642.30964, 470.34929;
                0,          0,         1 ];
K2 = [ 1646.07299,          0, 620.74483;
                0, 1645.39302, 477.47527;
                0,          0,         1 ];

inds = 2700:2800;
results = cell(100,1);
basedir = '/home/stephen/Documents/Research/data/kagaru_airborne/100831_155323_MultiCamera0_subset_db';
for i = 1:100
    Im1 = double(mean(...
        imread(fullfile(basedir, sprintf('cam0_image%05d.png',inds(i)))),3))/255;
    Im2 = double(mean(...
        imread(fullfile(basedir,sprintf('cam0_image%05d.png',inds(i+1)))),3))/255;
    flow = ImageFlow(Im1,Im2,'K',K1);
    c = CostFunctionFactory('ZhangTomasi',flow);
    results{i} = c.getResults;
end

save('results.mat','results')

