function plotSampleResults (samples, hit)

[f,x] = hist(samples, 25) ; 

bar(x, f, 'w') ; 
xlabel('Percentage of Significant SNP-Deletion Pairs')
ylabel('Number of Samples')
title('(non-)GWAS SNP-Deletion Pairs')
hold on 
line([hit,hit],[0, max(f) * 1.02], 'color', 'black', 'LineStyle', '-.', 'LineWidth', 2) 
legend('non-GWAS', 'GWAS')