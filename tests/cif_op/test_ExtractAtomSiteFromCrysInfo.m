%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2019 - 2020  Francis Black Lee and Li Xian

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.

%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

%   Email: warner323@outlook.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_ExtractAtomSiteFromCrysInfo.m
clc;
clear;
close all;
%% main:
rawData = [13      130     0.000000    0.000000    0.999885        0.00000     0.00000   -19.87112;
13      130     0.666667    0.333333    0.333218       -2.12940     0.97794    -6.62218;
13      130     0.333333    0.666667    0.666552       -0.21778     2.33308   -13.24665;
13      130     1.000000    1.000000    0.999885       -2.34718     3.31102   -19.87112;
13      130     1.000000    0.000000    0.999885       -4.04102    -0.37720   -19.87112;
13      130     0.000000    1.000000    0.999885        1.69384     3.68823   -19.87112;
31      311     1.000000    0.000000    0.500067       -4.04102    -0.37720    -9.93803;
31      311     0.000000    1.000000    0.500067        1.69384     3.68823    -9.93803;
31      311     0.666667    0.333333    0.833400       -2.12940     0.97794   -16.56250;
31      311     0.333333    0.666667    0.166734       -0.21778     2.33308    -3.31357;
31      311     0.000000    0.000000    0.500067        0.00000     0.00000    -9.93803;
31      311     1.000000    1.000000    0.500067       -2.34718     3.31102    -9.93803;
33      332     0.000000    1.000000    0.124797        1.69384     3.68823    -2.48014;
33      332     0.333333    0.666667    0.791464       -0.21778     2.33308   -15.72908;
33      332     1.000000    0.000000    0.124797       -4.04102    -0.37720    -2.48014;
33      332     0.666667    0.333333    0.458130       -2.12940     0.97794    -9.10461;
33      332     0.000000    0.000000    0.124797        0.00000     0.00000    -2.48014;
33      332     1.000000    1.000000    0.124797       -2.34718     3.31102    -2.48014;
33      333     0.333333    0.666667    0.291919       -0.21778     2.33308    -5.80142;
33      333     1.000000    1.000000    0.625252       -2.34718     3.31102   -12.42588;
33      333     0.666667    0.333333    0.958585       -2.12940     0.97794   -19.05035;
33      333     0.000000    1.000000    0.625252        1.69384     3.68823   -12.42588;
33      333     0.000000    0.000000    0.625252        0.00000     0.00000   -12.42588;
33      333     1.000000    0.000000    0.625252       -4.04102    -0.37720   -12.42588;];

xyz = rawData(:, 6 : 8)';
uvw = rawData(:, 3 : 5)';

convMat = [4.05858319000000,-2.02929159500000,0;0,3.51483614591249,0;0,0,19.8734024600000];

convedUvw = convMat * uvw;

figure;
subplot(1, 2, 1);
scatter(xyz(1, :), xyz(2, :), 'filled');
title('xyz');

subplot(1, 2, 2);
scatter(convedUvw(1, :), convedUvw(2, :), 'filled');
title('conv * uvw');