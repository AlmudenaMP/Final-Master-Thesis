%AUTHORS: ALMUDENA MÉNDEZ PÉREZ & JESÚS DAZA GARCÍA

%Scaling: 
%Time is rescaled in units of the mRNA lifetime, the number of repressors necessary to half-maximally repress a promoter 
%Protein concentrations are written in units of K_M
%mRNA concentrations are rescaled by their translation efficiency, the average number of proteins produced per mRNA molecule.

function v1 = v1(t,y)

%ATc concentration
%ATc = 0; % OFF state
ATc = 100000; % ON state  

%Parameters declaration 

%mRNA degradation constant (s^-1: mRNAs per second)
%Adapted from J. A. Bernstein et al., Proc. Natl Acad. Sci. USA 99:9697, 2002.
mtetR = 0.0030303;
mcI = 0.0030303; 
mProt = 0.0030303;
mLacI = 0.0030303;
mEspB = 0.0030303;

%mRNA synthesis without leakiness (transcripts per second). 
%Adapted from:
% -> Bremer, H., Dennis, P. P. (1996) Modulation of chemical composition and other parameters of the cell by growth rate. Neidhardt, et al. eds. Escherichia coli and Salmonella typhimurium: Cellular and Molecular Biology, 2nd ed. chapter 97 Table 3
% -> Vogel U, Jensen KF. The RNA chain elongation rate in Escherichia coli depends on the growth rate. J Bacteriol. 1994 May176(10):2807-13 p.2811 table 2PubMed ID7514589
% -> Proshkin S, Rahmouni AR, Mironov A, Nudler E. Cooperation between translating ribosomes and RNA polymerase in transcription elongation. Science. 2010 Apr 23 328(5977):504-8. p.505 table 1PubMed ID20413502
alphatetR = 0.0937499; 
alphacI = 0.0822784; 
alphaProt = 0.0228395; 
alphaLacI = 0.0531818; 
alphaEspB = 0.0607476; 

%mRNA synthesis by leaky promoter (transcripts per second)
%adjusted
ltetR = 0; 
lcI = 0; 
lProt = 0;
lLacI = 5*10^(-10); 
lEspB = 5*10^(-12); 

%Protein degradation constant (s^-1: degradated proteins per second). 
%adjusted
betatetR = 0.001155; 
betacI = 0.001155; 
betaProt = 0.001155;
betaLacI = 0.001155; 
betaEspB = 0.001155;

%Binding constant (M^-1).
KbATc = 1.26*10^(12); % -> Two mutations in the tetracycline repressor change the inducer anhydrotetracycline to a corepressor
KbtetR = 5.6*10^(9); % -> Two mutations in the tetracycline repressor change the inducer anhydrotetracycline to a corepressor
KbcI = 1.6*10^(9); % -> The Lysis-Lysogeny Decision of Bacteriophage 933W: a 933W Repressor-Mediated Long-Distance Loop Has No Role in Regulating 933W PRM Activity
KbLacI = 10^(11); % -> A  genetic  biosensor  for  identification of transcriptional repressors of target promoters

%mRNA translation efficiency (s^-1) 
%Adapted from (Michael B. Elowitz Nature, 2000)
gammatetR = 5; 
gammacI = 5; 
gammaProt = 5;
gammaLacI = 5; 
gammaEspB = 5; 

%Hill coefficient.
ntetR = 2; 
ncI = 2; 
nLacI = 2;  

% Protease-specific LacI degradation constant (s^-1).
B = 0.001155;

%EQUATIONS. 

if y(2) < 0.00
    f = 0;
else
    f = 1;   
end 

if y(4) < 0.00
    g = 0;
else
    g = 1;   
end 

if y(6) < 0.00
    h = 0;
else
    h = 1;   
end 

if y(10) < 0.00
    i = 0;
else
    i = 1;   
end 

%tetR
v1(1) = -mtetR*y(1) + alphatetR + ltetR; %mRNA
v1(2) = (-betatetR*y(2) + gammatetR*y(1) - KbATc*ATc)*f; %functional tetR protein

%cI
v1(3) = -mcI*y(3) + alphacI/(1+KbtetR*(y(2)^ntetR)) + lcI; %mRNA
v1(4) = (-betacI*y(4) + gammacI*y(3))*g; %protein

%LacI
v1(5) = -mLacI*y(5) + alphaLacI/(1+KbcI*(y(4)^ncI)) + lLacI; %mRNA
v1(6) = (-betaLacI*y(6) - B*y(6)*y(10) + gammaLacI*y(5))*h; %protein

%TSS3 reporter (assembled, functional T3SS) 
v1(7) = -mEspB*y(7) + alphaEspB/(1+KbLacI*(y(6)^nLacI)) + lEspB; %mRNA
v1(8) = -betaEspB*y(8) + gammaEspB*y(7); %protein

%PROTEASE
v1(9) = -mProt*y(9) + alphaProt/(1+KbtetR*(y(2)^ntetR)) + lProt; %mRNA
v1(10) = (-betaProt*y(10) + gammaProt*y(9))*i; %protein

v1=[v1(1); v1(2); v1(3); v1(4); v1(5); v1(6); v1(7); v1(8); v1(9); v1(10);];
end
