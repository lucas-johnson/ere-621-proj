
%{ 
    ERE 621 - Fall 2021
    Final Project
    Lucas Johnson
%}

clear all; close all; clc

%%%% FORMAT DATA %%%%

% Only work with 2019
target_year = 2019; 

% Training plots
plots = readmatrix("/Volumes/big_bag/data/CUI/ERE_621/lucas_stuff/training.csv");
plots = plots(plots(:, 2) == target_year, :);
plot_res = plots(:, [5, 6, 3, 4, 2]);
plot_res(:, 3) = plot_res(:, 4) - plot_res(:, 3);

% Test plots
test_plots = readmatrix("/Volumes/big_bag/data/CUI/ERE_621/lucas_stuff/testing.csv");
test_plots = test_plots(test_plots(:, 2) == target_year, :);
test_plot_res = test_plots(:, [5, 6, 3, 4, 2]);
test_plot_res(:, 3) = test_plot_res(:, 4) - test_plot_res(:, 3);


% Average distance between plots
nys = shaperead("../../data/nys_shape/state.shp");
nys_area = polyarea(nys.X(~isnan(nys.X)), nys.Y(~isnan(nys.Y))) / (1000^2);
plot_res_dist = dist_mat(plot_res);
min_dist = min(plot_res_dist(plot_res_dist > 0), [], 'all');
max_dist = max(plot_res_dist(plot_res_dist > 0), [], 'all');
intensity = size(plot_res, 1) / nys_area;
writematrix([1, 2, 3; min_dist, max_dist, intensity], '../../data/train_dist.csv');

test_dist = dist_mat(test_plot_res);
min_dist = min(test_dist(test_dist > 0), [], 'all');
max_dist = max(test_dist(test_dist > 0), [], 'all');
intensity = size(test_plot_res, 1) / nys_area;
writematrix([1, 2, 3; min_dist, max_dist, intensity], '../../data/test_dist.csv');

%%%% FIRST ORDER EFFECTS %%%%

% Check the distribution
figure('Name', 'Training Residuals Histogram');
histogram(plot_res(:, 3), 25);
ylabel("Count");
xlabel("Residual (AGB Mg ha^{-1})");

res_min = min(plot_res(:, 3));
res_max = max(plot_res(:, 3));

% rule of thumb for number of breaks
num_bins = round(1 + 3.3 * log(size(plot_res, 1)));

% Custom color bar to match num_bins
rwb = [
    1, 0, 0; 
    1, .1, .1;
    1, .2, .2; 
    1, .3, .3; 
    1, .4, .4; 
    1, .5, .5; 
    1, .6, .6; 
    1, .7, .7; 
    1, .8, .8;
    1, .9, .9;
    1, 1, 1;  
    .9, .9, 1;
    .8, .8, 1;
    .7, .7, 1;
    .6, .6, 1;
     .5, .5, 1;
    .4, .4, 1;
    .3, .3, 1;
    .2, .2, 1; 
    .1, .1, 1;
    0, 0, 1;
];

% Continuous Point Kernel and Map
bw = 25000;
[Xcent,Ycent,Kval] = johnson_continuous_point_kernel(plot_res, bw, 5000, 5000, 1);
figure('Name',['Image: Plot Residuals - Kernel Bandwidth=' num2str(bw)])
hold on
k_plot = imagesc(Xcent,Ycent,Kval);
set(k_plot, 'AlphaData', ~isnan(Kval))
set(gca,'YDir','normal');
set(gca,'color',[0.5, 0.5, 0.5]);
colormap(rwb)
caxis([-res_max, res_max]);
ylabel("Y (m)");
xlabel("X (m)");
z = colorbar;
set(get(z,'label'),'string',"Residual (AGB Mg ha^{-1})");
set(get(z,'label'),'rotation',270);
set(get(z,'label'),'position',[4, 0, 0]);
set(get(z,'label'),'FontSize',12);

%%%% FIT VARIOGRAM %%%%
[h, vgram] = johnson_variogram(plot_res, 500, 0, 1); 
ylabel('Variance');
xlabel('Distance (m)');

% Use cross validation to find best variogram parameters
rng(123);
nugget_test = 500:100:2000;
range_test = 15000:15000:75000;
sill_test = 500:250:2500;
model_test= [0, 1];
[nt, rt, st, mt] = ndgrid(nugget_test, range_test, sill_test, model_test);
combs = [nt(:), rt(:), st(:), mt(:)];

%results = grid_search_variogram(plot_res, combs);
%writematrix([1, 2, 3, 4, 5; results], '../../data/variogram_params.csv');
results = readmatrix("../../data/variogram_params.csv");
results = results(2:size(combs, 1), :);
results = sortrows(results, 1, 'ascend');
nugget = results(1, 2);
range = results(1, 3);
sill = results(1, 4);
model_type = results(1, 5);

% Best model parameters
test_cv = spherical_model(h, sill, range, nugget);

% Manual best fit
test_manual = exponential_model(h, 1850, 50000, 900);
manual_sill = 1850; 
manual_range = 50000;
manual_nugget = 900;
manual_model_type = 1;

figure('name', 'Variogram Models');
plot(h, vgram, '.');
ylabel('Variance');
xlabel('Distance (m)');
hold on;
plot(h, test_cv, '-r');
plot(h, test_manual, '-b');
hold off;
legend('Estimated Variogram','CV Fit - Spherical Model', 'Manual Fit - Exponential Model');
legend('Location','northwest');


%%%% PLOT KRIGING SURFACES %%%%
Mat = plot_res;
X = Mat(:,1);
Y = Mat(:,2);
Z = Mat(:,3); 

XYZ = [Mat(:, [1 2 3])];

%Create a 500*500 grid
n = 500; 
ti_x = [min(X):(max(X)-min(X))/n:max(X)];
ti_y = [min(Y):(max(Y)-min(Y))/n:max(Y)];

[temp  nx ] = size(ti_x);
[temp  ny ] = size(ti_y);

%Get all the possible coordinates in this grid.
X1 = repmat(ti_x,ny, 1);
X2  =reshape (X1,nx*ny,1);
Y1 = repmat(ti_y, 1, nx);
Y2  =reshape (Y1,nx*ny,1);
Matp =[X2  Y2];

[zp_man, ze_man] = johnson_ordinary_kriging(plot_res, Matp, manual_sill, manual_sill, manual_range, manual_nugget, manual_model_type);

zpman_out= reshape(zp_man,ny,nx);
zeman_out  =reshape (ze_man,ny,nx);

[XI,YI] = meshgrid(ti_x,ti_y);


figure ('name','Kriging Prediction');
man_krig_plot = imagesc(X2, Y2, zpman_out); hold on
set(gca,'YDir','normal');
set(gca,'color',[0.5, 0.5, 0.5]);
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(rwb)
caxis([-max(abs(min(zp_man)), abs(max(zp_man))), max(abs(min(zp_man)), abs(max(zp_man)))]);
ylabel("Y (m)");
xlabel("X (m)");
z = colorbar;
set(get(z,'label'),'string',"Predicted Residual (AGB Mg ha^{-1})");
set(get(z,'label'),'rotation',270);
set(get(z,'label'),'position',[4, 0, 0]);
set(get(z,'label'),'FontSize',12);
plot(test_plot_res(:, 1), test_plot_res(:, 2), '.', 'Color', "#000000")
hold off;

figure ('name','Kriging Variance');
man_var_plot = imagesc(X2, Y2, zeman_out); hold on
set(gca,'YDir','normal');
set(gca,'color',[0.5, 0.5, 0.5]);
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(hsv)
ylabel("Y (m)");
xlabel("X (m)");
z = colorbar;
set(get(z,'label'),'string',"Prediction Variance");
set(get(z,'label'),'rotation',270);
label_pos = get(z,'label').Position;
label_pos(1) = label_pos(1) + 1.5;
set(get(z,'label'),'position', label_pos);
set(get(z,'label'),'FontSize',12);
plot(test_plot_res(:, 1), test_plot_res(:, 2), '.', 'Color', "#000000")
hold off;

%%%% EVALUATE AT TEST DATA LOCATIONS %%%%

%%% Metrics & plots for original performance at test locations %%%
old_rmse = sqrt(mean(test_plot_res(:, 3) .^ 2));
old_mbe = mean(test_plot_res(:, 3));
old_r2 = 1 - (sum((test_plots(:, 4) - test_plots(:, 3)).^2)/sum((test_plots(:, 3) - mean(test_plots(:, 3))).^2));
writematrix([1, 2, 3; old_rmse, old_mbe, old_r2], '../../data/test_perf.csv');

%%% Metrics for CV Fit Improvements %%%
[cvpred_res, cvpred_err] = johnson_ordinary_kriging(plot_res, test_plot_res(:, [1, 2]), sill, sill, range, nugget, model_type);
cv_predictions = max(test_plots(:, 4) - cvpred_res, 0);
cv_res = cv_predictions - test_plots(:, 3); 
cv_rmse = sqrt(mean(cv_res .^ 2));
cv_mbe = mean(cv_res);
cv_r2 = 1 - (sum((cv_predictions - test_plots(:, 3)).^2)/sum((test_plots(:, 3) - mean(test_plots(:, 3))).^2));
writematrix([1, 2, 3; cv_rmse, cv_mbe, cv_r2], '../../data/cv_perf.csv');

%%% Metrics for Manual Fit Improvements %%%
[manpred_res, manpred_err] = johnson_ordinary_kriging(plot_res, test_plot_res(:, [1, 2]), manual_sill, manual_sill, manual_range, manual_nugget, manual_model_type);
man_predictions = max(test_plots(:, 4) - manpred_res, 0);
man_res = man_predictions - test_plots(:, 3); 
man_rmse = sqrt(mean(man_res .^ 2));
man_mbe = mean(man_res);
man_r2 = 1 - (sum((man_predictions - test_plots(:, 3)).^2)/sum((test_plots(:, 3) - mean(test_plots(:, 3))).^2));
writematrix([1, 2, 3; man_rmse, man_mbe, man_r2], '../../data/man_perf.csv');

%%% Scatter plot comparisons - Predictions vs Reference %%%

figure('Name', 'Scatter Comparisons');
hold on;
plot(test_plots(:, 3), test_plots(:, 4), '.', 'Color', [0.75, 0, 0]);
plot(test_plots(:, 3), cv_predictions, 'x', 'Color', [0, 0.75, 0]);
plot(test_plots(:, 3), man_predictions, 's', 'Color', [0, 0, 0.75]);
ref = refline(1, 0);
ref.Color = '#000000';
hold off;
ylabel('Predicted AGB (Mg ha^{-1})');
xlabel('Reference AGB (Mg ha^{-1})');
legend('Original Predictions', 'CV Fit Predictions', 'Manual Fit Predictions', '1:1 Line');
legend('Location','northwest');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% KRIGING FUNCTIONS %%%%

function [Zpr, Zpr_error] = johnson_ordinary_kriging(xy, pts, sigma_sq, sill, range, nugget, cov_function_shape)
    
    d_mat = dist_mat(xy);
    
    function cov = covariogram(h)
        if cov_function_shape == 0
            cov = sill - spherical_model(h, sill, range, nugget);
        else
            cov = sill - exponential_model(h, sill, range, nugget);
        end
    end

    for i = 1:size(d_mat, 1)
        for j = 1:size(d_mat, 1)
            COV1(i, j) = covariogram(d_mat(i, j));
        end
    end
    
    % Ordinary kriging transformation
    COV1(size(d_mat, 1) + 1, :) = ones(1, size(d_mat, 1));
    COV1(:, size(d_mat, 1) + 1) = ones(size(d_mat, 1)+1, 1);
    COV1(size(COV1, 1), size(COV1, 1)) = 0;
    COV1_inv = inv(COV1);
    
    Zpr = zeros(size(pts, 1), 1);
    Zpr_error = zeros(size(pts, 1), 1);
    for i = 1:size(pts, 1)
        c_pred = zeros(size(xy, 1), 1);
        for j = 1:size(xy)
            d = my_dist(pts(i, :), xy(j, :));
            c_pred(j, 1) = covariogram(d);
        end
        c_pred(size(c_pred, 1)+1, 1) = 1;
        w = COV1_inv * c_pred;
        w = w(1:size(w, 1)-1, 1);
        
        Zpr(i, 1) = sum(w.' * xy(:, 3));
        Zpr_error(i, 1) = sill - c_pred.' * COV1_inv * c_pred;
    end
    
end

function sphere_val = spherical_model(h, s, r, nug)
    sphere_val = nug + (s - nug) * (((3.*h)/(2*r)) - ((h.^3)/(2*r^3)));
    sphere_val(h == 0) = 0;
    sphere_val(h > r) = s;
end

function exp_val = exponential_model(h, s, r, nug)
    exp_val = nug + ((s - nug) * (1 - exp(-3.*h/r)));
    exp_val(h == 0) = 0;
end

function [bins, vgram] = johnson_variogram(xy, n_bins, do_plot, cloud_number)
    d_mat = dist_mat(xy);
    a_mat = attr_mat(xy);
    
    bins = build_bins(n_bins, d_mat);
    bins = bins(bins ~= 0).';
    t = abs(bins(1, 1) - bins(2, 1))/2;
    vgram = zeros(size(bins, 1), 1);
    for i = 1:size(bins, 1)
        if i == 1
            lower = 0;
            upper = bins(i, 1) + t;
        elseif i == size(bins, 2)
            lower = bins(i, 1) - t;
            upper = bins(i, 1);
        else 
            lower = bins(i, 1) - t;
            upper = bins(i, 1) + t;
        end
        
        mask = (d_mat >= lower & d_mat <= upper);
        vgram(i, 1) = mean(a_mat(mask == 1), 'all') / 2;
    end
    
    if do_plot
        figure('Name',['Variogram - num_bins =', num2str(n_bins)],'NumberTitle','on')
        plot(bins, vgram, '.');
        grid on
    end
    
    if cloud_number
        cloud = zeros(size(xy, 1), 2);
        upper_tri = triu(ones(size(d_mat, 1), size(d_mat, 2)));
        
        cloud_i = 1;
        
        for i = 1:size(xy, 1)
            for j = 1:size(xy, 1)
                if j == i
                    continue
                elseif upper_tri(i, j) == 0
                    continue
                else 
                    cloud(cloud_i, 1) = d_mat(i, j);
                    cloud(cloud_i, 2) = a_mat(i, j);
                    cloud_i = cloud_i + 1;                        

                end
            end
        end
        figure('Name', 'Variogram Cloud','NumberTitle','on')
        plot(cloud(:, 1), cloud(:, 2), '.');
        axis tight 
        grid on
    end
    
end

function a_mat = attr_mat(xy)
    n = size(xy, 1);
    a = ones(n);
    for i = 1:n
        for j = 1:n
            if j == i
                a(i, j) = inf;
            else 
                a(i, j) = (xy(i, 3) - xy(j, 3))^2;
            end
        end
    end
    a_mat = a;
end


%%%% KERNEL DENSITY FUNCTIONS %%%%

function [x_cent, y_cent, k_vals] = johnson_continuous_point_kernel(xy, bw, x_dist, y_dist, weight_kernel)
    [ksx, ksy] = kernel_coords(xy, x_dist, y_dist);
    
    
     for r = 1:size(ksy, 2)
        for c = 1:size(ksx, 2)
            
            % Loop on the points
            k = zeros(size(xy, 1), 1);
            yi = zeros(size(xy, 1), 1);
            for z = 1:size(xy, 1)
                
                cur_dist = my_dist([ksx(c), ksy(r)], [xy(z, 1), xy(z, 2)]);
                
                % Check the distance and compute the kernel
                if  cur_dist > bw
                    k(z) = 0;
                    yi(z) = 0;
                else
                    if weight_kernel
                        kern = my_kernel(bw, cur_dist);
                    else
                        kern = 1;
                    end
                    
                    yi(z) =  kern * xy(z, 3);
                    k(z) = kern;
                end
            end

            kv(r, c) = sum(yi)/sum(k);
        end
    end
    
    % Produce a columnar set of x,y coords
    centers = xy_lists_to_coords(ksx, ksy);
    x_cent = centers(:, 1);
    y_cent = centers(:, 2);
    k_vals = kv;

end

function [ksx, ksy] = kernel_coords(xy, x_dist, y_dist)
    % Build kernel center locations
    xy_max = max(xy);
    xy_min = min(xy);
    
    ksx = xy_min(1):x_dist:xy_max(1);
    ksy =  xy_max(2):-y_dist:xy_min(2);

end

function [coords] = xy_lists_to_coords(x, y)
    coords = ones(size(y, 2) * size(x, 2), 2);
    c = 1;
    for i = 1:size(x, 2)
        for j = 1:size(y, 2)
            coords(c, 1) = x(i);
            coords(c, 2) = y(j);
            c = c+1;
        end
    end

end

function k = my_kernel(bw, h)
    k = (3 / (pi * bw^2)) * (1 - (h^2/bw^2))^2;
end


%%%% VARIOGRAM CROSS VALIDATION FUNCTIONS %%%%

function [results] = grid_search_variogram(data, params)
    
    for i = 1:size(params, 1)
        param_rmse(i, 1) = k_fold_cv(data, 5, params(i, 1), params(i, 2), params(i, 3), params(i, 4));
    end
    results = [param_rmse, params];
end

function rmse = calc_rmse(x, y)
    rmse = sqrt(mean((x-y).^2));
end

function rmse = k_fold_cv(data, k, nug, range, sill, model_form)
    per_fold = floor(size(data, 1) / k);
    fold_order = randsample(1:size(data, 1), per_fold * k, false);
    fold_rows = reshape(fold_order,per_fold,[]);
    for i = 1:size(fold_rows, 2)
         fold_train = data(setdiff(1:length(fold_order), fold_rows(:, i)), :);
         fold_test = data(fold_rows(:, i), :);
         [zpr, zerr] = johnson_ordinary_kriging(fold_train, fold_test, sill, sill, range, nug, model_form);
         fold_rmse(i, 1) = calc_rmse(zpr, fold_test(:, 3));
    end
    rmse = mean(fold_rmse, 'all');
end


%%%% UTILITY FUNCTIONS %%%%

function bins = build_bins(n_bins, d)
    max_dist = max(d .* ~isinf(d), [], 'all');
    bins = 0:(max_dist/n_bins):max_dist;
end

function d = my_dist(pt1, pt2)
    x_dist = (pt1(1, 1) - pt2(1, 1));
    y_dist = (pt1(1, 2) - pt2(1, 2));
    d = sqrt(x_dist ^ 2 + y_dist ^ 2);
end

function d_mat = dist_mat(xy)
    % Make a distance matrix
    n = size(xy, 1);
    d = ones(n);
    for i = 1:n
        for j = 1:n
            if j == i
                d(i, j) = 0;
            else 
                d(i, j) = my_dist([xy(i, 1), xy(i, 2)], [xy(j, 1), xy(j, 2)]);
            end
        end
    end
    d_mat = d;
end