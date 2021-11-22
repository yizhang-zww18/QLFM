air_pos = zeros(3,5,45);
for j = 1: 3
for i = 1: 5
    filename = strcat('U:\YL\zyforYL\phase\2021_0412_air_bead',int2str(j),'_pos20_phasePic\zernike_para_layer_1_ite',int2str(i),'_rot180.mat');
    temp = load(filename);
    air_pos(j,i,:) = temp.a1;
end
end

water_pos = zeros(3,5,45);
for j = 1: 3
for i = 1: 5
    filename = strcat('U:\YL\zyforYL\phase\2021_0412_water_bead',int2str(j),'_pos20_phasePic\zernike_para_layer_1_ite',int2str(i),'_rot180.mat');
    temp = load(filename);
    water_pos(j,i,:) = temp.a1;
end
end

air_neg = zeros(3,5,45);
for j = 1: 3
for i = 1: 5
    filename = strcat('U:\YL\zyforYL\phase\2021_0412_air_bead',int2str(j),'_neg20_phasePic\zernike_para_layer_1_ite',int2str(i),'_rot180.mat');
    temp = load(filename);
    air_neg(j,i,:) = temp.a1;
end
end

water_neg = zeros(3,5,45);
for j = 1: 3
for i = 1: 5
    filename = strcat('U:\YL\zyforYL\phase\2021_0412_water_bead',int2str(j),'_neg20_phasePic\zernike_para_layer_1_ite',int2str(i),'_rot180.mat');
    temp = load(filename);
    water_neg(j,i,:) = temp.a1;
end
end

air_pos_sum = sum(air_pos,2);
water_pos_sum = sum(water_pos,2);
air_neg_sum = sum(air_neg,2);
water_neg_sum = sum(water_neg,2);


for i = 1:3
    temp = reshape(air_pos_sum(i,1,:),[1,45]);
    subplot(2,2,1),title('air pos');
    plot(temp);
    if i ~= 3
        hold on
    end
end

for i = 1:3
    temp = reshape(water_pos_sum(i,1,:),[1,45]);
    subplot(2,2,2),title('water pos');
    plot(temp);
    if i ~= 3
        hold on
    end
end

for i = 1:3
    temp = reshape(air_neg_sum(i,1,:),[1,45]);
    subplot(2,2,3),title('air neg');
    plot(temp);
    if i ~= 3
        hold on
    end
end

for i = 1:3
    temp = reshape(water_neg_sum(i,1,:),[1,45]);
    subplot(2,2,4),title('water neg');
    plot(temp);
    if i ~= 3
        hold on
    end
end