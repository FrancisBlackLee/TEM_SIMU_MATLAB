% create_space_group_table.m
clc;
clear;
close all;
%% space group Hermann-Mauguin notation and number in the table
symmetry_Int_Tables_number = uint8((1 : 230)');
space_group_IT_number = symmetry_Int_Tables_number;
space_group_H_M = readlines('space_group_H-M.txt');

%% space group crystal system
space_group_crystal_system = cell(230, 1);
% triclinic
for i = 1 : 2
    space_group_crystal_system{i} = 'triclinic';
end

% monoclinic
for i = 3 : 15
    space_group_crystal_system{i} = 'monoclinic';
end

% orthorhombic
for i = 16 : 74
    space_group_crystal_system{i} = 'orthorhombic';
end

% tetragonal
for i = 75 : 142
    space_group_crystal_system{i} = 'tetragonal';
end

% trigonal
for i = 143 : 167
    space_group_crystal_system{i} = 'trigonal';
end

% hexagonal
for i = 168 : 194
    space_group_crystal_system{i} = 'hexagonal';
end

% cubic
for i = 195 : 230
    space_group_crystal_system{i} = 'cubic';
end

space_group_crystal_system = string(space_group_crystal_system);

%% symmetry cell setting
symmetry_cell_setting = space_group_crystal_system;
symmetry_cell_setting(146) = 'rhombohedral';
symmetry_cell_setting(148) = 'rhombohedral';
symmetry_cell_setting(155) = 'rhombohedral';
symmetry_cell_setting(160) = 'rhombohedral';
symmetry_cell_setting(161) = 'rhombohedral';
symmetry_cell_setting(166) = 'rhombohedral';
symmetry_cell_setting(167) = 'rhombohedral';

%% make table
space_group_table = table(symmetry_Int_Tables_number, space_group_IT_number,...
    space_group_H_M, space_group_crystal_system, symmetry_cell_setting);

% writetable(space_group_table, 'cif_op\space_group_table.txt');

