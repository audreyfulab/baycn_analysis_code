Code repository for Martin et al. "Approximate Bayesian inference of directed acyclic graphs in biology 
with flexible priors on edge states".

The baycn package is available at https://github.com/evanamartin/baycn.

- simulation_analysis/ (true edges as input)

  - Goal: Evaluate inference accuracy when the true are used as input.

  - What it does:

    - Generates synthetic data under Gaussian distribution with different sample sizes or signal strengths.

    - Runs baycn and other methods.

    - Computes metrics (e.g., precision, recall/power, MSE between true and posterior probability adjacency matrices).

    - Produces comparison plots and summary tables.

- simulation_analysis_with_fully_connected_input/ (fully connected start)

  - Goal: Stress-test inference when starting from no prior structure.

  - What it does: Same as above, but the initial graph is fully connected.
    This highlights how baycn behaves with minimal knowledge of the graph topology.

- thinning/ 

  - Goal: Thin the MCMC samples from other methods such that these methods use the same MCMC settings as baycn.
    
- geuvadis_analysis/

  - Goal: Apply baycn to the GEUVADIS data of gene expression and genotypes in humans.

  - What it does: Data preprocessing, model fitting, and result visualization (e.g., correlation heatmaps, 
posterior edge probabilities and performance summaries).

- drosophila_analysis/

  - Goal: Apply baycn to Drosophila cis-regulatory modules (CRMs) and transcription factor (TF) binding data.

  - What it does: Prepares binary indicators for CRM expression patterns across tissues and continuous TF 
binding measurements across time points, runs baycn and other methods, and generates the figures/tables.
