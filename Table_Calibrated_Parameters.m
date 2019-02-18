% Creates a Latex Table with calibrated parameters
clc;
clear all;

% Define path where to store table
    my_path = '/home/arnau/Dropbox/Choi_Valladares_2015/QEresubmission/';
    name = 'table_parameters.tex';
% Read column vector with calibrated parameters
    fileID = fopen('calibrated.txt');
    aux = textscan(fileID, '%f %*s %*[^\n]');
    fclose(fileID);
    p = aux{1,1};
    clear aux

% Print to latex table

    name = [my_path name];
    fid = fopen(name, 'w');


    %fprintf(fid, '\\begin{table}[h!] \n');
    fprintf(fid, '\\begin{centering} \n');
    fprintf(fid, '\\begin{tabular}{clc} \n');
    fprintf(fid, 'Symbol  & Description  & Value\\tabularnewline \n');
    fprintf(fid, '\\hline \n');
    fprintf(fid, '\\hline \n');

    fprintf(fid, '$\\beta$  & Discount factor  & %5.4f \\tabularnewline \n', p(1));

    fprintf(fid, '\\multicolumn{3}{c}{} \\tabularnewline \n');
    fprintf(fid, '\\multicolumn{3}{l}{Single Male Households} \\tabularnewline \n');
    fprintf(fid, '\\hline \n');
    fprintf(fid, '$\\alpha_{\\mathcal{S},m}$  & Disutility of work  & %5.4f \\tabularnewline \n', p(2));
    fprintf(fid, '$\\epsilon^{\\gamma}_{\\mathcal{S},m}$  & S.d. search cost  & %5.4f \\tabularnewline \n', p(4));
    fprintf(fid, '$\\lambda^{u}_{\\mathcal{S},m}$  & Arrival probability unemployed & %5.4f \\tabularnewline \n', p(5));
    fprintf(fid, '$\\lambda^{n}_{\\mathcal{S},m}$  & Arrival probability OLF & %5.4f \\tabularnewline \n', p(6));
    fprintf(fid, '$\\sigma_{\\mathcal{S},m}$  & Losing probability & %5.4f \\tabularnewline \n', p(7));

    fprintf(fid, '\\multicolumn{3}{c}{} \\tabularnewline \n');
    fprintf(fid, '\\multicolumn{3}{l}{Single Female Households} \\tabularnewline \n');
    fprintf(fid, '\\hline \n');
    fprintf(fid, '$\\alpha_{\\mathcal{S},f}$  & Disutility of work  & %5.4f \\tabularnewline \n', p(8));
    fprintf(fid, '$\\epsilon^{\\gamma}_{\\mathcal{S},f}$  & S.d. search cost  & %5.4f \\tabularnewline \n', p(10));
    fprintf(fid, '$\\lambda^{u}_{\\mathcal{S},f}$  & Arrival probability unemployed & %5.4f \\tabularnewline \n', p(11));
    fprintf(fid, '$\\lambda^{n}_{\\mathcal{S},f}$  & Arrival probability OLF & %5.4f \\tabularnewline \n', p(12));
    fprintf(fid, '$\\sigma_{\\mathcal{S},f}$  & Losing probability & %5.4f \\tabularnewline \n', p(13));

    fprintf(fid, '\\multicolumn{3}{c}{} \\tabularnewline \n');
    fprintf(fid, '\\multicolumn{3}{l}{Married Households} \\tabularnewline \n');
    fprintf(fid, '\\hline \n');
    fprintf(fid, '$\\alpha_{\\mathcal{M},f}$  & Disutility of work male  & %5.4f \\tabularnewline \n', p(14));
    fprintf(fid, '$\\alpha_{\\mathcal{M},m}$  & Disutility of work female  & %5.4f \\tabularnewline \n', p(15));
    fprintf(fid, '$\\alpha_{\\mathcal{M}}$  & Disutility of joint work  & %5.4f \\tabularnewline \n', p(16));
    fprintf(fid, '$\\underbar{c}$  & Consumption floor  & %5.4f \\tabularnewline \n', p(17));
    fprintf(fid, '$\\epsilon^{\\gamma}_{\\mathcal{M},m}$  & S.d. search cost male  & %5.4f \\tabularnewline \n', p(19));
    fprintf(fid, '$\\epsilon^{\\gamma}_{\\mathcal{M},f}$  & S.d. search cost female  & %5.4f \\tabularnewline \n', p(20));
    fprintf(fid, '$\\lambda^{u}_{\\mathcal{M},m}$  & Arrival probability unemployed male & %5.4f \\tabularnewline \n', p(21));
    fprintf(fid, '$\\lambda^{u}_{\\mathcal{M},f}$  & Arrival probability unemployed female & %5.4f \\tabularnewline \n', p(22));
    fprintf(fid, '$\\lambda^{n}_{\\mathcal{M},m}$  & Arrival probability OLF male & %5.4f \\tabularnewline \n', p(23));
    fprintf(fid, '$\\lambda^{n}_{\\mathcal{M},f}$  & Arrival probability OLF female & %5.4f \\tabularnewline \n', p(24));
    fprintf(fid, '$\\sigma_{\\mathcal{M},m}$  & Losing probability male & %5.4f \\tabularnewline \n', p(25));
    fprintf(fid, '$\\sigma_{\\mathcal{M},f}$  & Losing probability female & %5.4f \\tabularnewline \n', p(26));

    fprintf(fid, '\\hline \n');
    fprintf(fid, '\\multicolumn{3}{c}{} \\tabularnewline \n');
    fprintf(fid, '\\end{tabular} \n');
    fprintf(fid, '\\par\\end{centering} \n');

    fclose(fid);
