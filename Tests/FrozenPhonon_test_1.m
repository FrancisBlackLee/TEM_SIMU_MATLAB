% FrozenPhonon_test_1.m -- test the configurations generated in FrozenPhonon_test_0.m
clc;
close all;
clear all;
%% Load configurations:
LoadPath = 'D:\Francis. B. Lee\Practice\Tomography\nvstgt\Si_110_diff_data\300K_configs\';
ConfigNum = 20;

for ConfigIdx = 1 : ConfigNum
    filename = strcat(LoadPath, 'Config_', num2str(ConfigIdx), '.txt');
    TempConfig = load(filename);
    Config{ConfigIdx} = TempConfig;
    figure;
    scatter(TempConfig( : , 4), TempConfig( : , 5));
end