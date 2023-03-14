clear; clc;
load('220222\GT_I_for_cross_sections.mat');

sinograms = zeros(400,180);
for ii=1:180
    sinograms(:,ii) = I_GT(1,1+mod((ii-1)*10+(-200:199),1800),40);
end

imshow(sinograms',[]);
%%
clear; clc;

first = load('220301/hu_images.mat');
second = load('220301/hu_images_2.mat');

high_100_kev_hu_images = cat(1,first.high_100_kev_hu_images,second.high_100_kev_hu_images);
high_150_kev_hu_images = cat(1,first.high_150_kev_hu_images,second.high_150_kev_hu_images);

high_100_kev_hu_images(high_100_kev_hu_images<-1000)=-1000;
high_150_kev_hu_images(high_150_kev_hu_images<-1000)=-1000;

[mu_water, ~] = PerEnergyMaterialsAttenuations({'Water'}, [100,150]);
[mu_air, ~] = PerEnergyMaterialsAttenuations({'Air'}, [100,150]);

high_100_kev_mu_images = mu_water(1)+((mu_water(1)-mu_air(1))*high_100_kev_hu_images/1000);
high_150_kev_mu_images = mu_water(2)+((mu_water(2)-mu_air(2))*high_150_kev_hu_images/1000);


%%
tmp_var = [permute(high_100_kev_mu_images,[2,3,1]),permute(high_150_kev_mu_images,[2,3,1])];
FindSolids(tmp_var)
