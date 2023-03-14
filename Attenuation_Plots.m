clear; clc;
gParams;

water_density_ranges = [ 0.95865, 1];
bone_density_ranges = [ 0.8, 1.1];

%%
energy_low = 50;
energy_high = 80;
water_mass_atten_low = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_low);
water_mass_atten_high = interp1((1e3)*water_atten(:,1),water_atten(:,3),energy_high);
bone_mass_atten_low = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_low);
bone_mass_atten_high = interp1((1e3)*bone_atten(:,1),bone_atten(:,3),energy_high);

water_lin_atten_low = water_density_ranges*water_mass_atten_low;
water_lin_atten_high = water_density_ranges*water_mass_atten_high;
bone_lin_atten_low = bone_density_ranges*bone_mass_atten_low;
bone_lin_atten_high = bone_density_ranges*bone_mass_atten_high;

% area(lin_atten_low,lin_atten_high)
water_pos=[water_lin_atten_low(1) water_lin_atten_high(1) diff(water_lin_atten_low) diff(water_lin_atten_high)];
rectangle('Position',water_pos,'FaceColor','r'); hold on;
text(water_lin_atten_low(2),water_lin_atten_high(2),'Water','Color','r')
bone_pos=[bone_lin_atten_low(1) bone_lin_atten_high(1) diff(bone_lin_atten_low) diff(bone_lin_atten_high)];
rectangle('Position',bone_pos,'FaceColor','b'); hold off;
text(bone_lin_atten_low(2),bone_lin_atten_high(2),'Bone','Color','b')

