clear;clc;

x_vec = [0.00 5.00e-03 1.00e-02 1.50e-02 2.00e-02 2.50e-02 3.00e-02 4.00e-02 5.00e-02 7.00e-02 9.00e-02 1.00e-01 1.25e-01 1.50e-01 1.75e-01 2.00e-01 2.50e-01 3.00e-01 4.00e-01 5.00e-01 6.00e-01 7.00e-01 8.00e-01 9.00e-01 1.00e+00 1.25e+00 1.50e+00 2.00e+00 2.50e+00 3.00e+00 3.50e+00 4.00e+00 5.00e+00 6.00e+00 7.00e+00 8.00e+00 1.00e+01 1.50e+01 2.00e+01 5.00e+01 8.00e+01 1.00e+02 1.00e+03 1.00e+06 1.00e+09];
F_x_C = [6.0000e+00 5.9974e+00 5.9898e+00 5.9711e+00 5.9594e+00 5.9368e+00 5.9093e+00 5.8406e+00 5.7544e+00 5.5369e+00 5.2702e+00 5.1225e+00 4.7407e+00 4.3310e+00 3.9371e+00 3.5775e+00 2.9614e+00 2.5015e+00 1.9512e+00 1.6856e+00 1.5353e+00 1.4245e+00 1.3206e+00 1.2165e+00 1.1121e+00 8.6482e-01 6.5662e-01 3.7202e-01 2.1465e-01 1.2832e-01 8.0452e-02 5.2230e-02 2.4330e-02 1.2650e-02 7.1471e-03 4.3194e-03 1.8365e-03 3.7767e-04 1.2157e-04 3.2386e-06 5.0734e-07 2.1123e-07 3.5964e-11 1.6609e-20 1.6810e-29];

E_kev = 60;
angstrem = 12398.52/(E_kev*1e3);

sin_theta_half = angstrem*x_vec;

%%
legend_vec = {};
for E_kev = [1,5,10,20,40,100]
    angstrem = 12398.52/(E_kev*1e3);
    r_e_squared = 7.9524e-26;
    theta_vec = linspace(0,pi,401);

    x_for_calc = sin(theta_vec/2)/angstrem;

    F_theta = interp1(x_vec,F_x_C,x_for_calc,'spline');
    rayleigh_cross = r_e_squared*0.5*(1+cos(theta_vec).^2).*(F_theta.^2);

    polarplot(theta_vec,rayleigh_cross); hold on
    legend_vec = cat(1,legend_vec,sprintf('Incident Energy = %d [kEv]',E_kev));
    thetalim([0,180])
end
legend(legend_vec,'Location','North')
set(gca,'FontSize',15);