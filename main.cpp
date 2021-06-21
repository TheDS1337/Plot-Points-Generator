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

//	plot_photon_statistics_n("photon_statistics_n_2_even.dat", &f_tsallis, 1.5, 1.0, LDBL_EPSILON, 0.0);
//  plot_quantum_visibility_t("visibility_t_2_0.dat", &f_tsallis, 1.5, 1.0, 0);
//	fidelity_superposition_t("fidelity_t_2_even2.dat", &f_tsallis, 1.5, 1.0, 0.0);
//	plot_quantum_visibility_n("visibility_n_2_0.dat", &f_tsallis, 1.5, 1.0, 0.0);



//	mandel_parameter_superposition2_t("mandel_parameter2_canonical.dat", &f_tsallis, 1.0 + LDBL_EPSILON, 1.0, M_PI);
//	mandel_parameter_superposition2_t("mandel_parameter2_harmonious.dat", &f_tsallis, 2.0 - LDBL_EPSILON, 0.9, M_PI);

//	second_order_correlation_parameter_superposition2_t("second_order_correlation_parameter2_canonical.dat", &f_tsallis, 1.0 + LDBL_EPSILON, 1.0, M_PI);
//	second_order_correlation_parameter_superposition2_t("second_order_correlation_parameter2_harmonious.dat", &f_tsallis, 2.0 - LDBL_EPSILON, 0.9, M_PI);

//	squeezing_factors_superposition2_t("squeezing_factors2_canonical.dat", &f_tsallis, 1.0 + LDBL_EPSILON, 1.0, M_PI);
//	squeezing_factors_superposition2_t("squeezing_factors2_harmonious.dat", &f_tsallis, 2.0 - LDBL_EPSILON, 0.9, M_PI);


//	squeezing_factors_superposition_alpha("squeezing_factors_alpha.dat", &f_tsallis, 2.0 - LDBL_EPSILON, LDBL_EPSILON);


	// SECTION 1
//	wigner("wigner_canonical_100.dat", &f_jackson, 1.0 - LDBL_EPSILON, 1.0);
//	wigner("wigner_harmonious_90.dat", &f_jackson, LDBL_EPSILON, 0.9);

//	wigner("wigner_jackson_5_25.dat", &f_jackson, 0.05, 0.25);
//	wigner("wigner_macfarlane_5_25.dat", &f_macfarlane, 0.05, 0.25);

//	wigner("wigner_jackson_30_50.dat", &f_jackson, 0.3, 0.5);
//	wigner("wigner_macfarlane_30_50.dat", &f_macfarlane, 0.3, 0.5);

//	wigner("wigner_jackson_80_75.dat", &f_jackson, 0.8, 0.75);
//	wigner("wigner_macfarlane_80_75.dat", &f_macfarlane, 0.8, 0.75);


//	wigner_superposition("wigner_canonical_100_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 1.0, 0.0);
//	wigner_superposition("wigner_harmonious_90_even.dat", &f_jackson, LDBL_EPSILON, 0.9, 0.0);

//	wigner_superposition("wigner_jackson_5_25_even.dat", &f_jackson, 0.05, 0.25, 0.0);
//	wigner_superposition("wigner_macfarlane_5_25_even.dat", &f_macfarlane, 0.05, 0.25, 0.0);

//	wigner_superposition("wigner_jackson_30_50_even.dat", &f_jackson, 0.3, 0.5, 0.0);
//	wigner_superposition("wigner_macfarlane_30_50_even.dat", &f_macfarlane, 0.3, 0.5, 0.0);

//	wigner_superposition("wigner_jackson_80_75_even.dat", &f_jackson, 0.8, 0.75, 0.0);
//	wigner_superposition("wigner_macfarlane_80_75_even.dat", &f_macfarlane, 0.8, 0.75, 0.0);


//	wigner_superposition("wigner_canonical_100_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 1.0, M_PI);
//	wigner_superposition("wigner_harmonious_90_odd.dat", &f_jackson, LDBL_EPSILON, 0.9, M_PI);

//	wigner_superposition("wigner_jackson_5_25_odd.dat", &f_jackson, 0.05, 0.25, M_PI);
//	wigner_superposition("wigner_macfarlane_5_25_odd.dat", &f_macfarlane, 0.05, 0.25, M_PI);

//	wigner_superposition("wigner_jackson_30_50_odd.dat", &f_jackson, 0.3, 0.5, M_PI);
//	wigner_superposition("wigner_macfarlane_30_50_odd.dat", &f_macfarlane, 0.3, 0.5, M_PI);

//	wigner_superposition("wigner_jackson_80_75_odd.dat", &f_jackson, 0.8, 0.75, M_PI);
//	wigner_superposition("wigner_macfarlane_80_75_odd.dat", &f_macfarlane, 0.8, 0.75, M_PI);


//	fidelity_superposition("fidelity_jackson_0_even.dat", &f_jackson, LDBL_EPSILON, 0.0);
//	fidelity_superposition("fidelity_macfarlane_0_even.dat", &f_macfarlane, LDBL_EPSILON, 0.0);

//	fidelity_superposition("fidelity_jackson_25_even.dat", &f_jackson, 0.25, 0.0);
//	fidelity_superposition("fidelity_macfarlane_25_even.dat", &f_macfarlane, 0.25, 0.0);

//	fidelity_superposition("fidelity_jackson_50_even.dat", &f_jackson, 0.5, 0.0);
//	fidelity_superposition("fidelity_macfarlane_50_even.dat", &f_macfarlane, 0.5, 0.0);

//	fidelity_superposition("fidelity_jackson_100_even.dat", &f_jackson, 1.0, 0.0);
//	fidelity_superposition("fidelity_macfarlane_100_even.dat", &f_macfarlane, 1.0, 0.0);


//	fidelity_superposition("fidelity_jackson_0_odd.dat", &f_jackson, LDBL_EPSILON, M_PI);
//	fidelity_superposition("fidelity_macfarlane_0_odd.dat", &f_macfarlane, LDBL_EPSILON, M_PI);

//	fidelity_superposition("fidelity_jackson_25_odd.dat", &f_jackson, 0.25, M_PI);
//	fidelity_superposition("fidelity_macfarlane_25_odd.dat", &f_macfarlane, 0.25, M_PI);

//	fidelity_superposition("fidelity_jackson_50_odd.dat", &f_jackson, 0.5, M_PI);
//	fidelity_superposition("fidelity_macfarlane_50_odd.dat", &f_macfarlane, 0.5, M_PI);

//	fidelity_superposition("fidelity_jackson_100_odd.dat", &f_jackson, 1.0, M_PI);
//	fidelity_superposition("fidelity_macfarlane_100_odd.dat", &f_macfarlane, 1.0, M_PI);

	// SECTION 2

	system("pause");
	return 0;
}