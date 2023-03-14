clear;clc;

Materials_names = {'Air','LungInhale','Adiposetissue','BreatMammary','Water','Muscle','Liver','dryribwater','Blood','Redmarrow'};

[linear_atten_mat, density_vec] = PerBinMaterialsAttenuations(Materials_names);

