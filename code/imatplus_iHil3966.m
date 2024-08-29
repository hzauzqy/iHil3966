addpath ./codeMERGE-master(iMAT++)/bins/
addpath ./code/MERGE-master(iMAT++)/input/
addpath ./code/MERGE-master(iMAT++)/1_iMAT++/scripts
initCobraToolbox(false);
load('iHil3966.mat');
model.csense = repmat('E', size(model.S, 1), 1);
model.dsense = repmat('E', length(model.ctrs), 1);

parsedGPR = GPRparser_xl(model);% Extracting GPR data from the model
model.parsedGPR = parsedGPR;
model = buildRxnGeneMat(model); 
model = creategrRulesField(model);
%[epsilon_f, epsilon_r] = makeEpsilonSeq2(model, model.rxns, 0.01, 0.5);
%save('.Input/epsilon_HL.mat','epsilon_f', 'epsilon_r');
load('epsilon_HL.mat');

%Fatty acid reaction
c12={'MAR00188','MAR09777','MAR10033','MAR00187','MAR10033_Human-GEM','MAR02172'}';
c12_index=ismember(model.rxns,c12);
c14={'MAR02177','MAR00196','MAR10033_Human-GEM','MAR10033','MAR00195','MAR09783'}';
c14_index=ismember(model.rxns,c14);
c16={'MAR02182','MAR09678','MAR10033_Human-GEM','MAR08939','MAR00216','MAR00217','MAR10033'}';
c16_index=ismember(model.rxns,c16);
c180={'MAR00248','MAR09792','MAR10033_Human-GEM','MAR08941','MAR02329','MAR00249','MAR10033'}';
c180_index=ismember(model.rxns,c180);
c181={'MAR09786','MAR00262','MAR10033_Human-GEM','MAR08940','MAR10033'}';
c181_index=ismember(model.rxns,c181);
c182={'MAR00396','MAR09779','MAR10033_Human-GEM','MAR00397','MAR10033'}';
c182_index=ismember(model.rxns,c182);

solution = optimizeCbModel(model, 'max');
f1=solution.f*0.9;

%¸÷Ä£ÐÍÔ¼Êø
model1=model;
model2=model;
model3=model;
model4=model;
model5=model;
model6=model;
model7=model;
model8=model;
%model1
model1.lb(ismember(model.rxns,{'Biomass'})) = f1;%biomass
model1.lb(ismember(model.rxns,{'MAR10033'})) = f1*0.05;%Fatty acids 279.816
model1.lb(ismember(model.rxns,{'MAR00187'})) = f1*0.05*0.72*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model1.lb(ismember(model.rxns,{'MAR00195'})) = f1*0.05*0.06*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model1.lb(ismember(model.rxns,{'MAR00216'})) = f1*0.05*0.02*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model1.lb(ismember(model.rxns,{'MAR00248'})) = f1*0.05*0.01*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model1.lb(ismember(model.rxns,{'MAR00262'})) = f1*0.05*0.05*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model1.lb(ismember(model.rxns,{'MAR00396'})) = f1*0.05*0.03*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model1.lb(ismember(model.rxns,{'MAR00397'})) = f1*0.05*0.03*(279.816/279.428);
model1.ub(ismember(model.rxns,{'MAR10033'})) = f1*0.05;%Fatty acids 279.816
model1.ub(ismember(model.rxns,{'MAR00187'})) = f1*0.05*0.72*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model1.ub(ismember(model.rxns,{'MAR00195'})) = f1*0.05*0.06*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model1.ub(ismember(model.rxns,{'MAR00216'})) = f1*0.05*0.02*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model1.ub(ismember(model.rxns,{'MAR00248'})) = f1*0.05*0.01*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model1.ub(ismember(model.rxns,{'MAR00262'})) = f1*0.05*0.05*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model1.ub(ismember(model.rxns,{'MAR00396'})) = f1*0.05*0.03*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model1.ub(ismember(model.rxns,{'MAR00397'})) = f1*0.05*0.03*(279.816/279.428);

%model2
model2.lb(ismember(model.rxns,{'Biomass'})) = f1;%biomass
model2.lb(ismember(model.rxns,{'MAR10033'})) = f1*0.08;%Fatty acids
model2.lb(ismember(model.rxns,{'MAR00187'})) = f1*0.08*0.08*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model2.lb(ismember(model.rxns,{'MAR00195'})) = f1*0.08*0.04*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model2.lb(ismember(model.rxns,{'MAR00216'})) = f1*0.08*0.21*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model2.lb(ismember(model.rxns,{'MAR00248'})) = f1*0.08*0.04*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model2.lb(ismember(model.rxns,{'MAR00262'})) = f1*0.08*0.34*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model2.lb(ismember(model.rxns,{'MAR00396'})) = f1*0.08*0.19*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model2.lb(ismember(model.rxns,{'MAR00397'})) = f1*0.08*0.19*(279.816/279.428);
model2.ub(ismember(model.rxns,{'MAR10033'})) = f1*0.08;%Fatty acids
model2.ub(ismember(model.rxns,{'MAR00187'})) = f1*0.08*0.08*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model2.ub(ismember(model.rxns,{'MAR00195'})) = f1*0.08*0.04*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model2.ub(ismember(model.rxns,{'MAR00216'})) = f1*0.08*0.21*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model2.ub(ismember(model.rxns,{'MAR00248'})) = f1*0.08*0.04*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model2.ub(ismember(model.rxns,{'MAR00262'})) = f1*0.08*0.34*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model2.ub(ismember(model.rxns,{'MAR00396'})) = f1*0.08*0.19*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model2.ub(ismember(model.rxns,{'MAR00397'})) = f1*0.08*0.19*(279.816/279.428);

%model3
model3.lb(ismember(model.rxns,{'Biomass'})) = f1;%biomass
model3.lb(ismember(model.rxns,{'MAR10033'})) = f1*0.18;%Fatty acids
model3.lb(ismember(model.rxns,{'MAR00187'})) = f1*0.18*0.38*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model3.lb(ismember(model.rxns,{'MAR00195'})) = f1*0.18*0.08*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model3.lb(ismember(model.rxns,{'MAR00216'})) = f1*0.18*0.1*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model3.lb(ismember(model.rxns,{'MAR00248'})) = f1*0.18*0.02*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model3.lb(ismember(model.rxns,{'MAR00262'})) = f1*0.18*0.15*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model3.lb(ismember(model.rxns,{'MAR00396'})) = f1*0.18*0.14*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model3.lb(ismember(model.rxns,{'MAR00397'})) = f1*0.18*0.14*(279.816/279.428);
model3.ub(ismember(model.rxns,{'MAR10033'})) = f1*0.18;%Fatty acids
model3.ub(ismember(model.rxns,{'MAR00187'})) = f1*0.18*0.38*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model3.ub(ismember(model.rxns,{'MAR00195'})) = f1*0.18*0.08*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model3.ub(ismember(model.rxns,{'MAR00216'})) = f1*0.18*0.1*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model3.ub(ismember(model.rxns,{'MAR00248'})) = f1*0.18*0.02*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model3.ub(ismember(model.rxns,{'MAR00262'})) = f1*0.18*0.15*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model3.ub(ismember(model.rxns,{'MAR00396'})) = f1*0.18*0.14*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model3.ub(ismember(model.rxns,{'MAR00397'})) = f1*0.18*0.14*(279.816/279.428);

%model4
model4.lb(ismember(model.rxns,{'Biomass'})) = f1;%biomass
model4.lb(ismember(model.rxns,{'MAR10033'})) = f1*0.24;%Fatty acids
model4.lb(ismember(model.rxns,{'MAR00187'})) = f1*0.24*0.52*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model4.lb(ismember(model.rxns,{'MAR00195'})) = f1*0.24*0.09*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model4.lb(ismember(model.rxns,{'MAR00216'})) = f1*0.24*0.08*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model4.lb(ismember(model.rxns,{'MAR00248'})) = f1*0.24*0.01*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model4.lb(ismember(model.rxns,{'MAR00262'})) = f1*0.24*0.1*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model4.lb(ismember(model.rxns,{'MAR00396'})) = f1*0.24*0.09*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model4.lb(ismember(model.rxns,{'MAR00397'})) = f1*0.24*0.09*(279.816/279.428);
model4.ub(ismember(model.rxns,{'MAR10033'})) = f1*0.24;%Fatty acids
model4.ub(ismember(model.rxns,{'MAR00187'})) = f1*0.24*0.52*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model4.ub(ismember(model.rxns,{'MAR00195'})) = f1*0.24*0.09*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model4.ub(ismember(model.rxns,{'MAR00216'})) = f1*0.24*0.08*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model4.ub(ismember(model.rxns,{'MAR00248'})) = f1*0.24*0.01*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model4.ub(ismember(model.rxns,{'MAR00262'})) = f1*0.24*0.1*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model4.ub(ismember(model.rxns,{'MAR00396'})) = f1*0.24*0.09*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model4.ub(ismember(model.rxns,{'MAR00397'})) = f1*0.24*0.09*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182

%model5
model5.lb(ismember(model.rxns,{'Biomass'})) = f1;%biomass
model5.lb(ismember(model.rxns,{'MAR10033'})) = f1*0.25;%Fatty acids
model5.lb(ismember(model.rxns,{'MAR00187'})) = f1*0.25*0.59*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model5.lb(ismember(model.rxns,{'MAR00195'})) = f1*0.25*0.09*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model5.lb(ismember(model.rxns,{'MAR00216'})) = f1*0.25*0.06*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model5.lb(ismember(model.rxns,{'MAR00248'})) = f1*0.25*0.01*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model5.lb(ismember(model.rxns,{'MAR00262'})) = f1*0.25*0.06*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model5.lb(ismember(model.rxns,{'MAR00396'})) = f1*0.25*0.08*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model5.lb(ismember(model.rxns,{'MAR00397'})) = f1*0.25*0.08*(279.816/279.428);
model5.ub(ismember(model.rxns,{'MAR10033'})) = f1*0.25;%Fatty acids
model5.ub(ismember(model.rxns,{'MAR00187'})) = f1*0.25*0.59*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model5.ub(ismember(model.rxns,{'MAR00195'})) = f1*0.25*0.09*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model5.ub(ismember(model.rxns,{'MAR00216'})) = f1*0.25*0.06*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model5.ub(ismember(model.rxns,{'MAR00248'})) = f1*0.25*0.01*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model5.ub(ismember(model.rxns,{'MAR00262'})) = f1*0.25*0.06*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model5.ub(ismember(model.rxns,{'MAR00396'})) = f1*0.25*0.08*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model5.ub(ismember(model.rxns,{'MAR00397'})) = f1*0.25*0.08*(279.816/279.428);

%model6
model6.lb(ismember(model.rxns,{'Biomass'})) = f1;%biomass
model6.lb(ismember(model.rxns,{'MAR10033'})) = f1*0.26;%Fatty acids
model6.lb(ismember(model.rxns,{'MAR00187'})) = f1*0.26*0.75*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model6.lb(ismember(model.rxns,{'MAR00195'})) = f1*0.26*0.09*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model6.lb(ismember(model.rxns,{'MAR00216'})) = f1*0.26*0.04*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model6.lb(ismember(model.rxns,{'MAR00248'})) = f1*0.26*0.01*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model6.lb(ismember(model.rxns,{'MAR00262'})) = f1*0.26*0.02*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model6.lb(ismember(model.rxns,{'MAR00396'})) = f1*0.26*0.05*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model6.lb(ismember(model.rxns,{'MAR00397'})) = f1*0.26*0.05*(279.816/279.428);
model6.ub(ismember(model.rxns,{'MAR10033'})) = f1*0.26;%Fatty acids
model6.ub(ismember(model.rxns,{'MAR00187'})) = f1*0.26*0.75*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model6.ub(ismember(model.rxns,{'MAR00195'})) = f1*0.26*0.09*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model6.ub(ismember(model.rxns,{'MAR00216'})) = f1*0.26*0.04*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model6.ub(ismember(model.rxns,{'MAR00248'})) = f1*0.26*0.01*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model6.ub(ismember(model.rxns,{'MAR00262'})) = f1*0.26*0.02*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model6.ub(ismember(model.rxns,{'MAR00396'})) = f1*0.26*0.05*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model6.ub(ismember(model.rxns,{'MAR00397'})) = f1*0.26*0.05*(279.816/279.428);

%model7
model7.lb(ismember(model.rxns,{'Biomass'})) = f1;%biomass
model7.lb(ismember(model.rxns,{'MAR10033'})) = f1*0.09;%Fatty acids
model7.lb(ismember(model.rxns,{'MAR00187'})) = f1*0.09*0.7*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model7.lb(ismember(model.rxns,{'MAR00195'})) = f1*0.09*0.08*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model7.lb(ismember(model.rxns,{'MAR00216'})) = f1*0.09*0.05*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model7.lb(ismember(model.rxns,{'MAR00248'})) = f1*0.09*0.01*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model7.lb(ismember(model.rxns,{'MAR00262'})) = f1*0.09*0.03*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model7.lb(ismember(model.rxns,{'MAR00396'})) = f1*0.09*0.05*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model7.lb(ismember(model.rxns,{'MAR00397'})) = f1*0.09*0.05*(279.816/279.428);
model7.ub(ismember(model.rxns,{'MAR10033'})) = f1*0.09;%Fatty acids
model7.ub(ismember(model.rxns,{'MAR00187'})) = f1*0.09*0.7*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model7.ub(ismember(model.rxns,{'MAR00195'})) = f1*0.09*0.08*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model7.ub(ismember(model.rxns,{'MAR00216'})) = f1*0.09*0.05*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model7.ub(ismember(model.rxns,{'MAR00248'})) = f1*0.09*0.01*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model7.ub(ismember(model.rxns,{'MAR00262'})) = f1*0.09*0.03*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model7.ub(ismember(model.rxns,{'MAR00396'})) = f1*0.09*0.05*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model7.ub(ismember(model.rxns,{'MAR00397'})) = f1*0.09*0.05*(279.816/279.428);

%model8
model8.lb(ismember(model.rxns,{'Biomass'})) = f1;%biomass
model8.lb(ismember(model.rxns,{'MAR10033'})) = f1*0.08;%Fatty acids
model8.lb(ismember(model.rxns,{'MAR00187'})) = f1*0.08*0.75*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model8.lb(ismember(model.rxns,{'MAR00195'})) = f1*0.08*0.08*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model8.lb(ismember(model.rxns,{'MAR00216'})) = f1*0.08*0.05*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model8.lb(ismember(model.rxns,{'MAR00248'})) = f1*0.08*0.01*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model8.lb(ismember(model.rxns,{'MAR00262'})) = f1*0.08*0.01*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model8.lb(ismember(model.rxns,{'MAR00396'})) = f1*0.08*0.02*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model8.lb(ismember(model.rxns,{'MAR00397'})) = f1*0.08*0.02*(279.816/279.428);
model8.ub(ismember(model.rxns,{'MAR10033'})) = f1*0.08;%Fatty acids
model8.ub(ismember(model.rxns,{'MAR00187'})) = f1*0.08*0.75*(279.816/199.304);%70 pp2;%80 HP2;%70;%50;%40; %8 4d;% C12
model8.ub(ismember(model.rxns,{'MAR00195'})) = f1*0.08*0.08*(279.816/227.356);%13;%10;%10;%10;%6;%2;% C14
model8.ub(ismember(model.rxns,{'MAR00216'})) = f1*0.08*0.05*(279.816/255.408);%4;%4;%4;%7;%10;%20;% C16
model8.ub(ismember(model.rxns,{'MAR00248'})) = f1*0.08*0.01*(279.816/283.46);%0;%0;%0;%1;%2;%4;% C180
model8.ub(ismember(model.rxns,{'MAR00262'})) = f1*0.08*0.01*(279.816/281.444);%2;%2;%2;%10;%15;%35;% C181
%model8.ub(ismember(model.rxns,{'MAR00396'})) = f1*0.08*0.02*(279.816/279.428);%3;%3;%3;%10;%15;%15;% C182
model8.ub(ismember(model.rxns,{'MAR00397'})) = f1*0.08*0.02*(279.816/279.428);
%%
folder_path = './Input/';
folder_path2 = './Result/';
files = {'categ_HL1d.mat','categ_HL4d.mat','categ_HL8d.mat','categ_HL12d.mat','categ_HPp1.mat','categ_HPp2.mat','categ_HP1.mat','categ_HP2.mat'};
%model_file={model1,model2,model3,model4,model5,model6,model7,model8};
myCSM1 = struct(); 
myCSM2 = struct(); 
myCSM3 = struct(); 
myCSM4 = struct(); 
myCSM5 = struct(); 
myCSM6 = struct(); 
myCSM7 = struct(); 
myCSM8 = struct(); 
%for i = 1:length(files)
    file_name = [folder_path,files{1}];
    file_name2 = [folder_path2,'myCSMHL1d.mat'];
    modelType = 3;
    %currentVar = eval(['model', num2str(1)]);
    %myCSM1 = eval(['myCSM', num2str(1)]);
    %currentVar2 = struct();
    load(file_name);
    [myCSM1.OFD,...
    myCSM1.PFD,...
    myCSM1.N_highFit,...
    myCSM1.N_zeroFit,...
    myCSM1.minLow,...
    myCSM1.minTotal,...
    myCSM1.minTotal_OFD,...
    myCSM1.MILP,...
    myCSM1.MILP_PFD,...
    myCSM1.HGenes,...
    myCSM1.RLNames,...
    myCSM1.OpenGene,...
    myCSM1.latentRxn,...
    myCSM1.Nfit_latent,...
    myCSM1.wasteDW]...
    = IMATplusplus(model1,epsilon_f,epsilon_r, ExpCateg, modelType);

save(file_name2,'myCSM1');
%end
file_name = [folder_path,files{2}];
    file_name2 = [folder_path2,'myCSMHL4d.mat'];
    modelType = 3;
    %currentVar = eval(['model', num2str(1)]);
    %myCSM1 = eval(['myCSM', num2str(1)]);
    %currentVar2 = struct();
    load(file_name);
    [myCSM2.OFD,...
    myCSM2.PFD,...
    myCSM2.N_highFit,...
    myCSM2.N_zeroFit,...
    myCSM2.minLow,...
    myCSM2.minTotal,...
    myCSM2.minTotal_OFD,...
    myCSM2.MILP,...
    myCSM2.MILP_PFD,...
    myCSM2.HGenes,...
    myCSM2.RLNames,...
    myCSM2.OpenGene,...
    myCSM2.latentRxn,...
    myCSM2.Nfit_latent,...
    myCSM2.wasteDW]...
    = IMATplusplus(model2,epsilon_f,epsilon_r, ExpCateg, modelType);

save(file_name2,'myCSM2');
%
file_name = [folder_path,files{3}];
    file_name2 = [folder_path2,'myCSMHL8d.mat'];
    modelType = 3;
    %currentVar = eval(['model', num2str(1)]);
    %myCSM1 = eval(['myCSM', num2str(1)]);
    %currentVar2 = struct();
    load(file_name);
    [myCSM3.OFD,...
    myCSM3.PFD,...
    myCSM3.N_highFit,...
    myCSM3.N_zeroFit,...
    myCSM3.minLow,...
    myCSM3.minTotal,...
    myCSM3.minTotal_OFD,...
    myCSM3.MILP,...
    myCSM3.MILP_PFD,...
    myCSM3.HGenes,...
    myCSM3.RLNames,...
    myCSM3.OpenGene,...
    myCSM3.latentRxn,...
    myCSM3.Nfit_latent,...
    myCSM3.wasteDW]...
    = IMATplusplus(model3,epsilon_f,epsilon_r, ExpCateg, modelType);

save(file_name2,'myCSM3');
%
file_name = [folder_path,files{4}];
    file_name2 = [folder_path2,'myCSMHL12d.mat'];
    modelType = 3;
    %currentVar = eval(['model', num2str(1)]);
    %myCSM1 = eval(['myCSM', num2str(1)]);
    %currentVar2 = struct();
    load(file_name);
    [myCSM4.OFD,...
    myCSM4.PFD,...
    myCSM4.N_highFit,...
    myCSM4.N_zeroFit,...
    myCSM4.minLow,...
    myCSM4.minTotal,...
    myCSM4.minTotal_OFD,...
    myCSM4.MILP,...
    myCSM4.MILP_PFD,...
    myCSM4.HGenes,...
    myCSM4.RLNames,...
    myCSM4.OpenGene,...
    myCSM4.latentRxn,...
    myCSM4.Nfit_latent,...
    myCSM4.wasteDW]...
    = IMATplusplus(model4,epsilon_f,epsilon_r, ExpCateg, modelType);

save(file_name2,'myCSM4');
%
file_name = [folder_path,files{5}];
    file_name2 = [folder_path2,'myCSMHPp1.mat'];
    modelType = 3;
    %currentVar = eval(['model', num2str(1)]);
    %myCSM1 = eval(['myCSM', num2str(1)]);
    %currentVar2 = struct();
    load(file_name);
    [myCSM5.OFD,...
    myCSM5.PFD,...
    myCSM5.N_highFit,...
    myCSM5.N_zeroFit,...
    myCSM5.minLow,...
    myCSM5.minTotal,...
    myCSM5.minTotal_OFD,...
    myCSM5.MILP,...
    myCSM5.MILP_PFD,...
    myCSM5.HGenes,...
    myCSM5.RLNames,...
    myCSM5.OpenGene,...
    myCSM5.latentRxn,...
    myCSM5.Nfit_latent,...
    myCSM5.wasteDW]...
    = IMATplusplus(model5,epsilon_f,epsilon_r, ExpCateg, modelType);

save(file_name2,'myCSM5');
%
file_name = [folder_path,files{6}];
    file_name2 = [folder_path2,'myCSMHPp2.mat'];
    modelType = 3;
    %currentVar = eval(['model', num2str(1)]);
    %myCSM1 = eval(['myCSM', num2str(1)]);
    %currentVar2 = struct();
    load(file_name);
    [myCSM6.OFD,...
    myCSM6.PFD,...
    myCSM6.N_highFit,...
    myCSM6.N_zeroFit,...
    myCSM6.minLow,...
    myCSM6.minTotal,...
    myCSM6.minTotal_OFD,...
    myCSM6.MILP,...
    myCSM6.MILP_PFD,...
    myCSM6.HGenes,...
    myCSM6.RLNames,...
    myCSM6.OpenGene,...
    myCSM6.latentRxn,...
    myCSM6.Nfit_latent,...
    myCSM6.wasteDW]...
    = IMATplusplus(model6,epsilon_f,epsilon_r, ExpCateg, modelType);

save(file_name2,'myCSM6');
%
file_name = [folder_path,files{7}];
    file_name2 = [folder_path2,'myCSMHP1.mat'];
    modelType = 3;
    %currentVar = eval(['model', num2str(1)]);
    %myCSM1 = eval(['myCSM', num2str(1)]);
    %currentVar2 = struct();
    load(file_name);
    [myCSM7.OFD,...
    myCSM7.PFD,...
    myCSM7.N_highFit,...
    myCSM7.N_zeroFit,...
    myCSM7.minLow,...
    myCSM7.minTotal,...
    myCSM7.minTotal_OFD,...
    myCSM7.MILP,...
    myCSM7.MILP_PFD,...
    myCSM7.HGenes,...
    myCSM7.RLNames,...
    myCSM7.OpenGene,...
    myCSM7.latentRxn,...
    myCSM7.Nfit_latent,...
    myCSM7.wasteDW]...
    = IMATplusplus(model7,epsilon_f,epsilon_r, ExpCateg, modelType);

save(file_name2,'myCSM7');
%
file_name = [folder_path,files{8}];
    file_name2 = [folder_path2,'myCSMHP2.mat'];
    modelType = 3;
    %currentVar = eval(['model', num2str(1)]);
    %myCSM1 = eval(['myCSM', num2str(1)]);
    %currentVar2 = struct();
    load(file_name);
    [myCSM8.OFD,...
    myCSM8.PFD,...
    myCSM8.N_highFit,...
    myCSM8.N_zeroFit,...
    myCSM8.minLow,...
    myCSM8.minTotal,...
    myCSM8.minTotal_OFD,...
    myCSM8.MILP,...
    myCSM8.MILP_PFD,...
    myCSM8.HGenes,...
    myCSM8.RLNames,...
    myCSM8.OpenGene,...
    myCSM8.latentRxn,...
    myCSM8.Nfit_latent,...
    myCSM8.wasteDW]...
    = IMATplusplus(model8,epsilon_f,epsilon_r, ExpCateg, modelType);

save(file_name2,'myCSM8');

%%Export Results
newTable = table(model.rxns,myCSM1.OFD, myCSM2.OFD, myCSM3.OFD, myCSM4.OFD, myCSM5.OFD, myCSM6.OFD, myCSM7.OFD, myCSM8.OFD);
newTable.Properties.VariableNames = {'rxn_id', 'HL1d', 'HL4d', 'HL8d', 'HL12d', 'HPp1', 'HPp2', 'HP1', 'HP2'};
writetable(newTable, './Result/geneCons.xlsx');



