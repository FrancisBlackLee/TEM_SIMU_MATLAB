function [namelist] = ParseSystemVariable(namelist, str)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

intVarList = ["ibrav", "nat", "ntyp", "nbnd", "nbnd_cond", "nr1", "nr2",...
    "nr3", "nr1s", "nr2s", "nr3s", "nspin", "nqx1", "nqx2", "nqx3", "edir",...
    "report", "esm_nfit", "dftd3_version", "space_group", "origin_choice",...
    "nextffield"];

realVarList = ["celldm(i)", "A", "B", "C", "cosAB", "cosAC", "cosBC",...
    "tot_charge", "starting_charge(i)", "tot_magnetization",...
    "starting_magnetization(i)", "ecutwfc", "ecutrho", "ecutfock",...
    "degauss_cond", "nelec_cond", "degauss", "sic_gamma", "sci_vb",...
    "sci_cb", "ecfixed", "qcutz", "q2sigma", "exx_fraction",...
    "screening_parameter", "ecutvcut", "localization_thr",...
    "Hubbard_occ(ityp,i)", "Hubbard_alpha(i)", "Hubbard_beta(i)",...
    "starting_ns_eigenvalue(m,ispin,ityp)", "emaxpos", "eopreg", "eamp",...
    "angle1(i)", "angle2(i)", "fixed_magnetization(i)", "lambda", "esm_w",...
    "esm_efield", "gcscf_mu", "gcscf_conv_thr", "gcscf_beta", "london_s6",...
    "london_c6(i)", "london_rvdw(i)", "london_rcut", "ts_vdw_econv_thr",...
    "xdm_a1", "xdm_a2", "zgate", "block_1", "block_2", "block_height"];

logicVarList = ["nosym", "nosym_evc", "noinv", "no_t_rev", "force_symmorphic",...
    "use_all_frac", "one_atom_occupations", "starting_spin_angle",...
    "sic_energy", "noncolin", "ace", "x_gamma_extrapolation", "dmft",...
    "ensemble_energies", "lforcet", "lspinorb", "lgcscf", "london",...
    "dftd3_threebody", "ts_vdw_isolated", "xdm", "uniqueb", "rhombohedral",...
    "relaxz", "block"];

charVarList = ["occupations", "smearing", "pol_type", "input_dft",...
    "exxdiv_treatment", "dmft_prefix", "constrained_magnetization",...
    "assume_isolated", "esm_bc", "vdw_corr"];

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