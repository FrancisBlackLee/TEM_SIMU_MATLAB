function [namelist] = ParseElectronsVariable(namelist, str)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

intVarList = ["electron_maxstep", "exx_maxstep", "mixing_ndim",...
    "mixing_fixed_ns", "diago_cg_maxiter", "diago_ppcg_maxiter",...
    "diago_david_ndim", "diago_rmm_ndim", "diago_gs_nblock"];

logicVarList = ["scf_must_converge", "adaptive_thr", "diago_rmm_conv",...
    "diago_full_acc", "tqr", "real_space"];

realVarList = ["conv_thr", "conv_thr_init", "conv_thr_multi", "mixing_beta",...
    "diago_thr_init", "efield", "efield_cart(i)"];

charVarList = ["mixing_mode", "diagonalization", "efield_phase",...
    "startingpot", "startingwfc"];

% make sure comment is not parsed
strs = split(str, '=');
if numel(strs) >= 2
    varStr = strtrim(strs{1});
    valueStr = strtrim(strs{2});
else
    varStr = "invalid";
    valueStr = "invalid";
end

varPat = replaceBetween(varStr, "(", ")", "");
varPat = erase(varPat, "(");
varPat = erase(varPat, ")");

varLabel = erase(varStr, " ");
varLabel = strrep(varLabel, "(", "_");
varLabel = strrep(varLabel, ",", "_");
varLabel = strrep(varLabel, ")", "");

valuePat = erase(valueStr, "(");
valuePat = erase(valuePat, ")");

if any(contains(intVarList, varPat))
    namelist.(varLabel) = str2num(valuePat);
elseif any(contains(realVarList, varPat))
    namelist.(varLabel) = str2double(valuePat);
elseif any(contains(logicVarList, varPat))
    if strcmpi(valuePat, ".true.")
        namelist.(varLabel) = true;
    elseif strcmpi(valuePat, ".false.")
        namelist.(varLabel) = false;
    else
        error('Invalid variable');
    end
elseif any(contains(charVarList, varPat))
    namelist.(varLabel) = valuePat;
end

end