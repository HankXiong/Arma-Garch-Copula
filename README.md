# Arma-Garch-Copula

Built ARMA-GARCH-Copula model to model the dependence between SP500 and TSX log return from 2006 - 2018

1. By comparing the different copula's average distance with empirical copula, I find that t copula fits the dependence best which means
that the extremes are more likely to happen and contrary to most results, in relatively short time, the dependece is not asymmetric.

2. Simulate returns based on the model, estimate VaR and compare their performance within different models and with traditional
VaR method through back-test. The comparison confirms the superiority of the Arma-Garch-Clayton copula in estimating VaR.

3. Please see more details about copula by reading attached my math thesis 
'Copula and Its Application in Estimating Portfolio Value-at-risk'
