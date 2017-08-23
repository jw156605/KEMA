addpath('general_routine')
load('simulated_data1')
domain1_labeled = struct('X',X1,'Y',Y1);
domain2_labeled = struct('X',X2,'Y',Y2);
labeled = {domain1_labeled,domain2_labeled};

domain1_un = struct('X',U1,'Y',Y1U);        
domain2_un = struct('X',U2,'Y',Y2U);      
unlabeled = {domain1_un,domain2_un};

options=struct('numDomains',2,'kernelt','rbf','printing',1,'lambda',0,'projections',500);
[evects,evals,K,opts]=KMA(labeled,unlabeled,options);

domain1_test = struct('X',XT1,'Y',YT1);
domain2_test = struct('X',XT2,'Y',YT2);
test = {domain1_test,domain2_test};
[proj] = KMAproject(labeled,unlabeled,test,evects,opts);
P1 = proj{1,1}.train;
P2 = proj{1,2}.train;
%P1 = proj{1,1}.test;
%P2 = proj{1,2}.test;

close all;
scatter(P1(1,:)',P1(2,:)');
hold on
scatter(P2(1,:)',P2(2,:)');
print('plot1','-dpng');

close all                      
scatter(P1(1,:)',P1(2,:)',[],YT1);
hold on                           
scatter(P2(1,:)',P2(2,:)',[],YT2);
colormap(jet);                    
print('plot2','-dpng');