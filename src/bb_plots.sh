#!/bin/sh

echo "Giving plots for title (a.k.a. \$1) = $1"

set -ev

./create_bb_plots $1

mkdir -p plots

./make_plot.py --figtitle="BB_Correlation_function_$1" --xlabel="Tau (s)" --ylabel="Correlation" data/bb_nba_pe_$1_correlation_fn.txt data/bb_nma_pe_$1_correlation_fn.txt data/bb_ta_pe_$1_correlation_fn.txt data/bb_fma_pe_$1_correlation_fn.txt data/bb_fba_pe_$1_correlation_fn.txt data/bb_total_pe_$1_correlation_fn.txt data/bb_config_$1.txt

./make_plot.py --figtitle="BB_Locally averaged PE_vs_time_$1" --ymax=10  --xlabel="Runtime (s)" --ylabel="PE / 0.5*kb*T" --hline=1.0 data/bb_nba_pe_$1.txt data/bb_nma_pe_$1.txt data/bb_ta_pe_$1.txt data/bb_fma_pe_$1.txt data/bb_fba_pe_$1.txt data/bb_total_pe_$1.txt data/bb_config_$1.txt

./make_plot.py --figtitle="BB_PE_average_vs_time_$1" --xlabel="Runtime (s)" --ymax=10 --ylabel="PE / 0.5*kb*T" --hline=1.0 data/bb_nba_pe_$1_eq_ave.txt data/bb_nma_pe_$1_eq_ave.txt data/bb_ta_pe_$1_eq_ave.txt data/bb_fma_pe_$1_eq_ave.txt data/bb_fba_pe_$1_eq_ave.txt data/bb_total_pe_$1_eq_ave.txt data/bb_config_$1.txt

./make_plot.py --logx --logy --figtitle="BB_Log_error_vs_log_time_$1" --xlabel="log(iterations)" --ylabel="log(| PE / ET - 1|)" --hline=1.0 data/bb_nba_pe_$1_log_error.txt data/bb_nma_pe_$1_log_error.txt data/bb_ta_pe_$1_log_error.txt data/bb_fma_pe_$1_log_error.txt data/bb_fba_pe_$1_log_error.txt data/bb_total_pe_$1_log_error.txt data/bb_config_$1.txt

./make_plot.py --figtitle="BB_Locally averaged angle_vs_time_$1" --xlabel="Runtime (s)" --ylabel="Angle" data/bb_nba_angle_$1.txt data/bb_nma_angle_$1.txt data/bb_ta_angle_$1.txt data/bb_fma_angle_$1.txt data/bb_fba_angle_$1.txt data/bb_config_$1.txt

./make_plot.py --figtitle="BB_Angle_n_PE_$1" --xlabel="Runtime (s)" --ylabel="Angle/PE" data/bb_fma_angle_$1.txt data/bb_fma_pe_$1.txt data/bb_ta_angle_$1.txt data/bb_ta_pe_$1.txt data/bb_config_$1.txt
