function [m_hit, m_nonhit] = plotHist (hit, metric, xl)

m_hit = [] ;
m_nonhit = [] ;

for i = 1:numel(hit)
   if hit(i) == 1
      m_hit = [m_hit; metric(i)] ;
   else
      m_nonhit = [m_nonhit; metric(i)] ;
   end
end

subplot(2,2,1)
hist(m_hit, 100, 'FaceColor', 'r')
title('GWAS SNPs')
xlabel(xl)
ylabel('# SNP-Deletion Pairs')
subplot(2,2,2)
hist(m_nonhit, 100)
title('Non-GWAS SNPs')
xlabel(xl)
ylabel('# SNP-Deletion Pairs')

subplot(2,2,3)
cdfplot(m_nonhit)
hold on; h=cdfplot(m_hit); hold off
set(h,'color','r')
title('Cumulative Distribution')
xlabel(xl)
ylabel('')
legend('Non-GWAS SNPs','GWAS SNPs','Location','SE')