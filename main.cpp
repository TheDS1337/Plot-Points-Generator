#include "photon_statistics.h"
#include "mandel_factor.h"
#include "second_order_correlation_parameter.h"
#include "squeezing_factors.h"
#include "wigner.h"
#include "visibility.h"
#include "fidelity.h"

int main()
{
/*		FOR BENCHMARK PURPOSES
*
*	auto start = std::time(NULL);
*	run some code here
*		auto end = std::time(NULL);
*
*	std::cout << "Started at " << time_stamp(localtime_xp(start), "%T") << std::endl;
*	std::cout << "Finished at " << time_stamp(localtime_xp(end), "%T") << std::endl;
*
*	auto duration = difftime(end, start) / 60.0;
*	std::cout << "Duration: " << duration << " minutes" << std::endl;
*/

//	plot_photon_statistics_n(&f_jackson, 0.95, 1.0, 1.0, M_PI);
//	plot_photon_statistics_t(&f_jackson, 0.5, 1.0, 0, M_PI);
//	plot_quantum_visibility_n(&f_jackson, 0.5, 1.0, 1.0); 
//	plot_quantum_visibility_t(&f_jackson, 0.5, 1.0, 0);
//	plot_fidelity_t("fidelity_t.dat", &f_jackson, 0.5, 1.2, M_PI);
//	plot_fidelity_alpha("fidelity_alpha.dat", &f_jackson, 1 - LDBL_EPSILON, 1.0, M_PI);

//	plot_i1(1.0 + LDBL_EPSILON, M_PI / 2, M_PI / 2 );
//	plot_i1(1.5 + LDBL_EPSILON);

//	wigner("wigner_canonical.dat", &f_tsallis, 1.0 - LDBL_EPSILON, 1.0);
//	wigner("wigner_harmonious.dat", &f_tsallis, 2.0 - LDBL_EPSILON, 0.99);
//	wigner_superposition("wigner_superposition_canonical.dat", &f_tsallis, 1.0 - LDBL_EPSILON, 1.0, M_PI / 4);
//	wigner_superposition("wigner_superposition_harmonious.dat", &f_tsallis, 2.0 - LDBL_EPSILON, 0.99, M_PI / 4);

//	plot_photon_statistics_n("photon_statistics_n_2_even.dat", &f_tsallis, 1.5, 1.0, LDBL_EPSILON, 0.0);
//  plot_quantum_visibility_t("visibility_t_2_0.dat", &f_tsallis, 1.5, 1.0, 0);
//	fidelity_superposition_t("fidelity_t_2_even2.dat", &f_tsallis, 1.5, 1.0, 0.0);
//	plot_quantum_visibility_n("visibility_n_2_0.dat", &f_tsallis, 1.5, 1.0, 0.0);






//	plot_mandel_parameter_t("mandel_parameter_t.dat", &f_tsallis, 1.0 + LDBL_EPSILON, 1.0, M_PI);
//	plot_second_order_correlation_parameter_alpha("second_order_correlation_parameter_alpha.dat", &f_tsallis, 1.0 + LDBL_EPSILON, LDBL_EPSILON, false);

//	plot_mandel_parameter_alpha("mandel_parameter_alpha.dat", &f_tsallis, 1.0 + LDBL_EPSILON, LDBL_EPSILON);

//	plot_squeezing_factors_t("squeezing_factors_t.dat", &f_tsallis, 1.0 + LDBL_EPSILON, LDBL_EPSILON, M_PI / 4);
//	squeezing_factors_superposition_alpha("squeezing_factors_alpha.dat", &f_tsallis, 2.0 - LDBL_EPSILON, LDBL_EPSILON);

	system("pause");
	return 0;
}