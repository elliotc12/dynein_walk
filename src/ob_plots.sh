#!/bin/sh

echo "Giving plots for title (a.k.a. \$1) = $1"

set -ev

./create_ob_plots $1

mkdir -p plots

./make_plot.py --figtitle="OB_Correlation_function_$1" --xlabel="Tau (s)" --ylabel="Correlation" data/ob_bba_pe_$1_correlation_fn.txt data/ob_bma_pe_$1_correlation_fn.txt data/ob_ta_pe_$1_correlation_fn.txt data/ob_uma_pe_$1_correlation_fn.txt data/ob_total_pe_$1_correlation_fn.txt data/ob_config_$1.txt

./make_plot.py --figtitle="OB_Locally averaged PE_vs_time_$1" --ymax=10  --xlabel="Runtime (s)" --ylabel="PE / 0.5*kb*T" --hline=1.0 data/ob_bba_pe_$1.txt data/ob_bma_pe_$1.txt data/ob_ta_pe_$1.txt data/ob_uma_pe_$1.txt data/ob_total_pe_$1.txt data/ob_config_$1.txt

./make_plot.py --figtitle="OB_PE_average_vs_time_$1" --xlabel="Runtime (s)" --ymax=10 --ylabel="PE / 0.5*kb*T" --hline=1.0 data/ob_bba_pe_$1_eq_ave.txt data/ob_bma_pe_$1_eq_ave.txt data/ob_ta_pe_$1_eq_ave.txt data/ob_uma_pe_$1_eq_ave.txt data/ob_total_pe_$1_eq_ave.txt data/ob_config_$1.txt

./make_plot.py --logx --logy --figtitle="OB_Log_error_vs_log_time_$1" --xlabel="log(iterations)" --ylabel="log(| PE / ET - 1|)" --hline=1.0 data/ob_bba_pe_$1_log_error.txt data/ob_bma_pe_$1_log_error.txt data/ob_ta_pe_$1_log_error.txt data/ob_uma_pe_$1_log_error.txt data/ob_total_pe_$1_log_error.txt data/ob_config_$1.txt

./make_plot.py --figtitle="OB_Locally averaged angle_vs_time_$1" --xlabel="Runtime (s)" --ylabel="Angle" --ymax=3.0 data/ob_bba_angle_$1.txt data/ob_bma_angle_$1.txt data/ob_ta_angle_$1.txt data/ob_uma_angle_$1.txt data/ob_config_$1.txt

./make_plot.py --figtitle="OB_Angle_n_PE_$1" --xlabel="Runtime (s)" --ylabel="Angle/PE" data/ob_uma_angle_$1.txt data/ob_uma_pe_$1.txt data/ob_ta_angle_$1.txt data/ob_ta_pe_$1.txt data/ob_config_$1.txt

./make_plot.py --figtitle="OB_Force_x_$1" --xlabel="Runtime (s)" --ylabel="Internal force" data/ob_bba_force_$1_x.txt data/ob_bma_force_$1_x.txt data/ob_ta_force_$1_x.txt data/ob_uma_force_$1_x.txt data/ob_uba_force_$1_x.txt data/ob_total_force_$1_x.txt data/ob_config_$1.txt

./make_plot.py --figtitle="OB_Force_y_$1" --xlabel="Runtime (s)" --ylabel="Internal force" data/ob_bba_force_$1_y.txt data/ob_bma_force_$1_y.txt data/ob_ta_force_$1_y.txt data/ob_uma_force_$1_y.txt data/ob_uba_force_$1_y.txt data/ob_total_force_$1_y.txt data/ob_config_$1.txt

