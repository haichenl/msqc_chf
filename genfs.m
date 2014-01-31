clear classes;
filebasename1 = 'etch';
fullfilename1 = strcat('C:\working\matlab\msqc_chf\fitting\', filebasename1, '.mat');
% filebasename2 = 'etch';
% fullfilename2 = strcat('C:\working\matlab\msqc_conshf\fitting\', filebasename2, '.mat');

fs = {};
xref = zeros(1,20);
x = zeros(1,20);
for i=1:20
    fs{i} = FitSingle(fullfilename1,i);
    fs{i}.hl.solvetranshf();
    xref(i) = fs{i}.hl.Ehf;
    x(i) = fs{i}.hl.newEhf;
end
plot(xref);
hold;
plot(x, 'r');
% fs2 = FitSingle(fullfilename2,7);
% fs1.ll.solvehft();
% fs1.hl.solvehft();
% fs1.solvehf();
% fs1.ll.solvetranshf();
% fs1.hl.solvetranshf();
% fs2.ll.solvehft();
% fs2.hl.solvehft();
% fs2.solvehf();
% fs2.ll.solvetranshf();
% fs2.hl.solvetranshf();
