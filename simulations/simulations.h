void get_onebound_PE_correlation_function(double* tau_data, double* corr_data, int iterations, int d_iter, int max_tau_iter);
void get_onebound_equipartition_ratio_per_runtime(double* runtime_data, double* eq_data, int d_runtime_iter, int min_runtime_iter, int max_runtime_iter);

void print_data_to_file(double* data1, double* data2, int iterations, const char* legend, const char* fname);
