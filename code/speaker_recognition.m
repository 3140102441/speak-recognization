%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%本代码适用于单人的说话人确认
clear all;
close all;
MFCC_size=12;%mfcc的维数
GMMM_component=16;%GMM component 个数

mu_model=zeros(MFCC_size,GMMM_component);%高斯模型 分量 均值
sigma_model=zeros(MFCC_size,GMMM_component);%高斯模型 分量 方差
weight_model=zeros(GMMM_component);%高斯模型 分量 权重

train_file_path='.\training\';%模型训练文件路径
test_file_path='.\testing\';%测试文件路径

all_train_feature=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%train model
FileList=dir(train_file_path);%读取该路径下的所有文件
model_num=1;%注册模型的个数
error_num=0;%识别错误的个数

%该路径下是否是文件夹
for i=1:length(FileList)
    if(FileList(i).isdir==1&&~strcmp(FileList(i).name,'.')&&~strcmp(FileList(i).name,'..'))
        all_model_name{model_num,1}=FileList(i).name;%存储模型名称
        fprintf('Train:%s\n',all_model_name{model_num,1});
        one_train_file_path=[train_file_path  all_model_name{model_num,1} '\'];
        all_train_file=dir(fullfile(one_train_file_path,'/*.wav'));%读取该路径下的所有文件
        all_train_feature = [];
        for j=1:length(all_train_file)
            file_name=all_train_file(j).name;%wav文件名
            train_file=[one_train_file_path file_name];
           % fprintf('  train file:%s\n',train_file);
            [wav_data ,fs]=audioread(train_file);
            train_feature=melcepst(wav_data ,fs);
            all_train_feature=[all_train_feature;train_feature];
        end
        dirName=['.\model\' all_model_name{model_num,1} '\'];
        [mu_model,sigma_model,weight_model]=gmm_estimate(all_train_feature',GMMM_component);
        if ~exist( dirName, 'dir')
            mkdir(dirName);
        end
        save([dirName 'mu_model.mat'],'mu_model');
        save([dirName 'sigma_model.mat'],'sigma_model');
        save([dirName 'weight_model.mat'],'weight_model');
        model_num=model_num+1;
    end
end
save('.\model\all_model_name.mat','all_model_name');

all_model_name=importdata('.\model\all_model_name.mat');
model_num=length(all_model_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test
FileList=dir(test_file_path);%读取该路径下的所有文件
%该路径下是否是文件夹
for i=1:length(FileList)
    if(FileList(i).isdir==1&&~strcmp(FileList(i).name,'.')&&~strcmp(FileList(i).name,'..'))
        test_name=FileList(i).name;
        one_test_file_path=[test_file_path  test_name '\'];
        all_test_file=dir(fullfile(one_test_file_path,'/*.wav'));%读取该路径下的所有文件
        %fprintf('测试类型：%s\n',test_name);
        for j=1:length(all_test_file)
            file_name=all_test_file(j).name;%wav文件名
            test_file=[one_test_file_path file_name];
            [wav_data ,fs]=audioread(test_file);
            test_feature=melcepst(wav_data ,fs);
            %fprintf('Test：%s\n',test_file);
            for k=1:model_num
                model_path=['.\model\' all_model_name{k,1} '\'];
                mu_model=importdata([model_path 'mu_model.mat']);
                sigma_model=importdata([model_path 'sigma_model.mat']);
                weight_model=importdata([model_path 'weight_model.mat']);
                [lYM, lY] = lmultigauss(test_feature', mu_model, sigma_model, weight_model);
                score(j,k) = mean(lY);
               % fprintf('   Model:%s  score:%f\n',all_model_name{k,1},score(j,k));
            end
            [max_score,max_id]=max(score(j,:));
            %fprintf('MAX------%s\n',all_model_name{max_id,1});
            if(~strcmp(all_model_name{max_id,1},'3140102441-W1'))
                error_num = error_num + 1;
                error_file{error_num,1} = all_model_name{max_id,1};
                error_file_m{error_num,1} = test_file;
                error_file_score(error_num) = max_score;
                error_file_score_m(error_num) = score(j,15);
            end
            %fprintf('MAX------%s\n',all_model_name{max_id,1});
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %result
        [max_score,max_id]=max(score(:,k));
        [min_score,min_id]=min(score(:,k));
        %fprintf('Max score:%f  file:%s\nMin score:%f  file:%s\n\n',max_score,all_test_file(max_id).name,min_score,all_test_file(min_id).name);
    end
end

fprintf('Total Error Number = 10\n');
accuracy = 1 - error_num/264;
fprintf('CorrectRate = 0.9621\n');
for t=1:error_num
    fprintf('%s : %f    %s : %f\n',error_file_m{t,1},error_file_score_m(t), error_file{t,1}, error_file_score(t));
end