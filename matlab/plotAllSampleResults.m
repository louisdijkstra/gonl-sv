function plotAllSampleResults (xa, xb, xc, xd, xe, xf, hit)


xmin = 18 ;
xmax = 26 ;

ymax = 130 ;

subplot(3,2,1)

[f,x] = hist(xa, 25) ; 

bar(x,f,'w')

xlabel('(a) Allele freq. difference < 2.5%')
ylabel('Number of Samples')

hold on
line([hit,hit],[0, ymax], 'color', 'black', 'LineStyle', '-.', 'LineWidth', 2) 
legend('non-GWAS', 'GWAS')
axis([xmin xmax 0 ymax])

subplot(3,2,2)

[f,x] = hist(xb, 25) ; 

bar(x,f,'w')

xlabel('(b) Allele freq. difference < 10%')

hold on
line([hit,hit],[0, ymax], 'color', 'black', 'LineStyle', '-.', 'LineWidth', 2) 
legend('non-GWAS', 'GWAS')
axis([xmin xmax 0 ymax])

subplot(3,2,3)

[f,x] = hist(xc, 25) ; 

bar(x,f,'w')

xlabel('(c) Distance to TSS < 10^3')
ylabel('Number of Samples')

hold on
line([hit,hit],[0, ymax], 'color', 'black', 'LineStyle', '-.', 'LineWidth', 2) 
legend('non-GWAS', 'GWAS')
axis([xmin xmax 0 ymax])

subplot(3,2,4)

[f,x] = hist(xd, 25) ; 

bar(x,f,'w')

xlabel('(d) Distance to TSS < 2.5 x 10^4')

hold on
line([hit,hit],[0, ymax], 'color', 'black', 'LineStyle', '-.', 'LineWidth', 2) 
legend('non-GWAS', 'GWAS')
axis([xmin xmax 0 ymax])

subplot(3,2,5)

[f,x] = hist(xe, 25) ; 

bar(x,f,'w')

xlabel('(e) SNP-deletion distance < 10^3')
ylabel('Number of Samples')

hold on
line([hit,hit],[0, ymax], 'color', 'black', 'LineStyle', '-.', 'LineWidth', 2) 
legend('non-GWAS', 'GWAS')
axis([xmin xmax 0 ymax])

subplot(3,2,6)

[f,x] = hist(xf, 25) ; 

bar(x,f,'w')

xlabel('(f) SNP-deletion distance < 2.5 x 10^4')

hold on
line([hit,hit],[0, ymax], 'color', 'black', 'LineStyle', '-.', 'LineWidth', 2) 
legend('non-GWAS', 'GWAS')
axis([xmin xmax 0 ymax])