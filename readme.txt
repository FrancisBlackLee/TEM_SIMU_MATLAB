A simple version of multislice simulation of BF-CTEM and ADF-STEM, for more details, please refer to
Advanced computing in electron microscopy, E. J. Kirkland.

Please add the folder to your MATLAB path, note that Scattering_Factors.txt must be included as well.

Function and script list:

propTF_1 -- modified Fresnel propagation function;
InteractionCoefficient -- used to calculate the interaction coefficient;
ProjectedPotential -- used to generate the 2D projected potential distribution of the input series of atoms;
AberrationFunction -- used to generate the transfer function of the objective lens, in parallel beam mode or convergent
		beam mode;
ProbeCreate -- used to generate a convergent electron beam;
BandwidthLimit -- used to symmetrically limit the bandwidth of a transfer function, i.e. the transmission function;
BesselAberration -- an adapted function from AberrationFunction.m used to generate Bessel electron beam together with
		BesselBeam.m, just like generate a convergent electron beam with AberrationFunction.m and ProbeCreate/m;
multislice -- perform the multislice process of one electron beam;
CTEM_sample_0 and STEM_sample_0 -- samples for demonstration.