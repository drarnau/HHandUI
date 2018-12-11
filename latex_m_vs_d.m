function latex_m_vs_d(d_stocks,d_trans,stocks,trans,my_path,t_ag)
%latex_m Prints a tex file with a table comparing data and model outcomes
%in Choi and Valladares-Esteban

if t_ag == 1
    name = 'table_sf.tex';
    title = 'Single Females';
elseif t_ag == 2
    name = 'table_sm.tex';
    title = 'Single Males';
elseif t_ag == 3
    name = 'table_mf.tex';
    title = 'Married Females';
elseif t_ag == 4
    name = 'table_mm.tex';
    title = 'Married Males';
end

name = [my_path name];
fid = fopen(name, 'w');

% Percentage 
d_trans = d_trans*100;
d_stocks = d_stocks*100;
trans = round((trans*10000))/10000; %Keep only 2 decimal points
trans = trans*100;
stocks = round((stocks*10000))/10000;
stocks = stocks*100;

%fprintf(fid, '\\begin{table}[h!] \n');
fprintf(fid, '\\begin{centering} \n');
fprintf(fid, '\\begin{tabular}{cr@{\\extracolsep{0pt}.}lr@{\\extracolsep{0pt}.}lr@{\\extracolsep{0pt}.}l|cr@{\\extracolsep{0pt}.}lr@{\\extracolsep{0pt}.}lr@{\\extracolsep{0pt}.}l} \n');
fprintf(fid, '\\multicolumn{7}{c|}{\\textbf{Data}} & \\multicolumn{7}{c}{\\textbf{Model}}\\tabularnewline \n');
fprintf(fid, '\\hline \n'); 
fprintf(fid, '\\hline \n'); 
fprintf(fid, '{\\scriptsize{}From / To}  & \\multicolumn{2}{c}{$E$ } & \\multicolumn{2}{c}{$U$ } & \\multicolumn{2}{c|}{$N$ } & {\\scriptsize{}From / To} & \\multicolumn{2}{c}{$E$ } & \\multicolumn{2}{c}{$U$ } & \\multicolumn{2}{c}{$N$} \\tabularnewline \n');
fprintf(fid, '$E$  & %s & %s & %s & $E$  & %s & %s & %s \\tabularnewline \n', Lnum(d_trans(1,1,t_ag)), Lnum(d_trans(1,2,t_ag)), Lnum(d_trans(1,3,t_ag)), Lnum(trans(1,1,t_ag)), Lnum(trans(1,2,t_ag)), Lnum(trans(1,3,t_ag)));   
fprintf(fid, '$U$  & %s & %s & %s & $U$  & %s & %s & %s \\tabularnewline \n', Lnum(d_trans(2,1,t_ag)), Lnum(d_trans(2,2,t_ag)), Lnum(d_trans(2,3,t_ag)), Lnum(trans(2,1,t_ag)), Lnum(trans(2,2,t_ag)), Lnum(trans(2,3,t_ag)));
fprintf(fid, '$N$  & %s & %s & %s & $N$  & %s & %s & %s \\tabularnewline \n', Lnum(d_trans(3,1,t_ag)), Lnum(d_trans(3,2,t_ag)), Lnum(d_trans(3,3,t_ag)), Lnum(trans(3,1,t_ag)), Lnum(trans(3,2,t_ag)), Lnum(trans(3,3,t_ag)));
fprintf(fid, '\\multicolumn{7}{c|}{} & \\multicolumn{7}{c}{} \\tabularnewline \n');
% fprintf(fid, '\\multicolumn{5}{l}{Employment Rate} & %s & \\multicolumn{5}{l}{Employment Rate} & %s \\tabularnewline \n', Lnum(d_stocks(1,t_ag)), Lnum(stocks(1,t_ag)));
fprintf(fid, '\\multicolumn{5}{l}{Unemployment Rate} & %s & \\multicolumn{5}{l}{Unemployment Rate} & %s \\tabularnewline \n', Lnum(d_stocks(2,t_ag)), Lnum(stocks(2,t_ag)));
fprintf(fid, '\\hline \n');
fprintf(fid, '\\multicolumn{14}{c}{} \\tabularnewline \n');
fprintf(fid, '\\end{tabular} \n');
fprintf(fid, '\\par\\end{centering} \n');
%fprintf(fid, '\\protect\\caption{{\\footnotesize{}\\label{tab:%s}Data versus Model (\\%%). %s. CPS 2000-2005.}} \n', title, title);
%fprintf(fid, '\\end{table} \n');

fclose(fid);

end

