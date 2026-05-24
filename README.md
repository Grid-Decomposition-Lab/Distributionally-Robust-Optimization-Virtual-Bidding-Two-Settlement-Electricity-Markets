# Distributionally-Robust-Optimization-Virtual-Bidding-Two-Settlement-Electricity-Markets
Independent reproductions of a data-driven virtual bidding strategy designed for two-settlement wholesale electricity markets.   Virtual bidding (or convergence bidding) allows market participants to arbitrage price differences between the Day-Ahead Market (DAM) and Real-Time Market (RTM). 
## Overview
This repository provides a data-driven simulation framework for Virtual Bidding in two-settlement electricity markets (e.g., NYISO). Distributionally Robust Optimization (DRO) strategy is reproduced, which enables market participants to arbitrage price spreads between Day-Ahead and Real-Time markets while managing severe volatility and distributional uncertainty.
The solver is built upon a tractable conic reformulation of the Wasserstein-metric-based DRO model, integrated with Conditional Value at Risk (CVaR) to handle tail risks and ensure robust profitability under extreme price uncertainties.
## Technical Highlights
- **Wasserstein-based Ambiguity Set**: Instead of assuming a fixed probability distribution, the model constructs a Wasserstein ball around historical empirical data to remain robust against price spread forecasting errors.
- **Integrated Risk Management (CVaR)**: The objective function incorporates a risk-aversion factor ($\rho$) and the CVaR measure to explicitly penalize extreme losses, protecting the portfolio during high-volatility periods.
- **Tailored Gurobi Tuning**: The solver is configured with interior-point methods (Method 2) and disabled crossover to ensure deterministic performance and rapid convergence for large-scale SOCP reformulations.
## Numerical Observations and Discussion
- **Sensitivity to Parameters ($\rho, \epsilon$)**: The strategy's performance, measured by the Calmar and Sharpe ratios, exhibits a non-monotonic relationship with the risk-aversion factor $\rho$ and the Wasserstein radius $\epsilon$. An excessively large $\epsilon$ (robustness) leads to overly conservative bids, while an extremely low $\rho$ may expose the portfolio to unacceptable drawdowns. The optimal balance depends on the underlying market's tail-risk characteristics.
- **Bidding Volume Distribution**. The total bidding volume constraint $||q_t||_1<=L$ works like Lasso algorithm, while a larger Wasserstein ball radius $\epsilon$ impose a more distributed bidding acrros different zones.
## References
[1] X. Audet, K. Qako, and A. Lesage-Landry, "A distributionally robust optimization strategy for virtual bidding in two-settlement electricity markets," Sustainable Energy, Grids and Networks, vol. 43, p. 101904, 2025. 
