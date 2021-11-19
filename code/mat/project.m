
%{ 
    ERE 621 - Fall 2021
    Final Project
    Lucas Johnson
%}

clear all; close all; clc

year = 2019; % Only work with one year for now

%%%% FIT VARIOGRAM %%%%
plots = readmatrix("/Volumes/big_bag/data/CUI/ERE_621/training_plots.csv");
plot_res = plots(:, [5, 6, 3, 4, 2]);
plot_res = plot_res(plot_res(:, 5) == 2019, :);
plot_res(:, 3) = plot_res(:, 4) - plot_res(:, 3);
[h, vgram] = johnson_variogram(plot_res, 500, 0, 0); 

sill = std(plot_res(:, 3))^2;
test_sphere = spherical_model(h, sill, 50000, 1000);
test_exp = exponential_model(h, sill, 50000, 1000);
figure('name', 'Variogram Models');
plot(h, vgram, '.');

hold on;
plot(h, test_sphere, '-r');
plot(h, test_exp, '-b');
hold off;

range = 50000;
nugget = 1000;
model_type = 1;

%%%% ORDINARY KRIGING %%%%
%Load dataset of Sudan
Mat = plot_res;
X = Mat(:,1);
Y = Mat(:,2);
Z = Mat(:,3); 

XYZ = [Mat(:, [1 2 3])];

%Create a 1000*1000 grid
n = 1000; % You can change the size of the grid if you want
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


[zpr, zerr] = johnson_ordinary_kriging(plot_res, Matp, sill, sill, range, nugget, model_type);

zpr_out=reshape (zpr,ny,nx);
zerr_out  =reshape (zerr,ny,nx);

[XI,YI] = meshgrid(ti_x,ti_y);

figure ('name','Ordinary Kriging');

subplot(2,2,1)
mesh(XI,YI,zpr_out) ; hold on
colormap hsv
stem3(X,Y,Z,'fill'); hold off
title('Original Z with kriging surface')

subplot(2,2,2)
mesh(XI,YI,zerr_out) ; hold on
colormap hsv
stem3(X,Y,Z,'fill'); hold off
title('Kriging Square of Standard deviation')

subplot(2,2,3)
[C,h] = contour(XI,YI,zpr_out); text_handle = clabel(C,h); set(text_handle,'BackgroundColor',[1 1 .6],'Edgecolor',[.7 .7 .7]); hold on; colormap hsv
plot3(X,Y,Z,'o'), hold off
title('Contour of Original Z with kriging surface')

subplot(2,2,4)
[C,h] = contour(XI,YI,zerr_out); text_handle = clabel(C,h); set(text_handle,'BackgroundColor',[1 1 .6],'Edgecolor',[.7 .7 .7]); hold on; colormap hsv
plot3(X,Y,Z,'o'), hold off
title('Contour of Kriging Square of Standard deviation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
