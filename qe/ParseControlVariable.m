function [namelist] = ParseControlVariable(namelist, str)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

charVarList = ["calculation", "title", "verbosity", "restart_mode",...
    "outdir", "wfcdir", "prefix", "disk_io", "pseudo_dir"];
logicVarList = ["tstress", "tprnfor", "tefield", "dipfield", "lelfield",...
    "lorbm", "lberry", "gate", "twochem", "lfcp", "trism"];
intVarList = ["nstep", "iprint", "nberrycyc", "gdir", "nppstr"];
realVarList = ["dt", "max_seconds", "etot_conv_thr", "forc_conv_thr"];

% make sure comment is not parsed
strs = split(str, '=');
if numel(strs) >= 2
    varStr = strtrim(strs{1});
    valueStr = strtrim(strs{2});
else
    varStr = "invalid";
    valueStr = "invalid";
end

if any(contains(charVarList, varStr))
    namelist.(varStr) = valueStr;
elseif any(contains(logicVarList, varStr))
    if strcmpi(valueStr, '.true.')
        namelist.(varStr) = true;
    elseif strcmpi(valueStr, '.false.')
        namelist.(varStr) = false;
    else
        error('Invalid variable');
    end
elseif any(contains(intVarList, varStr))
    namelist.(varStr) = str2num(valueStr);
elseif any(contains(realVarList, varStr))
    namelist.(varStr) = str2double(valueStr);
end

end