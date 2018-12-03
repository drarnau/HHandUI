% Creates Latex tables that compare Data vs. Model
clc;
clear all;

% Define path where to store table
    my_path = '/home/arnau/Dropbox/Choi_Valladares_2015/QEresubmission/';

% Read data file
    d_stocks = NaN(2,4); % E-to-pop, U-rate; SF, SM, Married F, MM
    d_trans = NaN(3,3,4);
    
    fid = fopen('data.txt'); 
    for ty = 1:4
        for j = 1:3 % Rows
        for i = 1:3 % Columns
            tline = fgets(fid, 7);
            d_trans(j,i,ty) = str2num(tline);
            fgets(fid);
        end
        end
%       % Read E/Pop
        tline = fgets(fid, 7);
        d_stocks(1,ty) = str2num(tline);
        fgets(fid);
        % Read Urate
        tline = fgets(fid, 7);
        d_stocks(2,ty) = str2num(tline);
        fgets(fid);
    end
    fclose(fid);
    
% Read model file
    stocks = NaN(2,4); % E-to-pop, U-rate; SF, SM, Married F, MM
    trans = NaN(3,3,4);
    
    fid = fopen('model.txt'); 
    for ty = 1:4
        for j = 1:3 % Rows
        for i = 1:3 % Columns
            tline = fgets(fid);
            trans(j,i,ty) = str2num(tline);
        end
        end
%       % Read E/Pop
        tline = fgets(fid);
        stocks(1,ty) = str2num(tline);
        % Read Urate
        tline = fgets(fid);
        stocks(2,ty) = str2num(tline);
    end
    fclose(fid);

% Create LaTeX tables comparing model and data
    
    for t_ag = 1:4
        latex_m_vs_d(d_stocks,d_trans,stocks,trans,my_path,t_ag);
    end



    
    



