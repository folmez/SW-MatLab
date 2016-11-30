# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) 
and this project adheres to [Semantic Versioning](http://semver.org/).

## 0.0.2 - 2016-11-30
### Changed
- Minor changes in "generate_random_graph.m"

### Added
- "do_SW_drift_diffusion_analysis.m" is added. All related parts of "SW_network.m" are moved there.
- "fit_PL_and_EXP_to_SW_bouts.m" is added. All parts of "SW_network.m" related to bounded power-law and exponential distribution fitting to SW bout durations are moved there. 
- "interval_to_bout.m" was in the RBM folder. It is now copied into SW-MatLab and renamed to "compute_bout_durations.m".
- "scale_free_graph_generation.m", "determine_alpha.m" and "chung_lu.m" are moved into a new folder titled "generate_network_configuration/SF_Chung_Lu".
- "generate_random_graph.m" and "generate_inhibitory_links.m" are moved into a new folder titled "generate_network_configuration".
- "compute_wake_sleep_domains_NEW.m", "compute_wake_sleep_domains.m", "compute_duration_matrix.m", "calculate_ints_tpb.m" are moved into a new folder titled "compute_activity_domains"
- "compute_activity_domains.m" is added. All related parts of "SW_network.m are moved there.
- "record_sim.m" is added. All parts of "SW_network.m" that are related to recording simulation output in a single structure are moved there.
- "generate_network_configuration.m" is added. All graph generation related parts of "SW_network.m" are now moved there.

## 0.0.1 - 2016-11-30
### Added
- This CHANGELOG file will be updated as my MatLab code for the Sleep-Wake project is updated.

