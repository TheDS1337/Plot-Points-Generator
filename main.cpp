#include "photon_statistics.h"
#include "mandel_factor.h"
#include "second_order_correlation_parameter.h"
#include "squeezing_factors.h"
#include "heisenberg.h"
#include "wigner.h"
#include "visibility.h"
#include "fidelity.h"

int main()
{
/*		FOR WIGNER FUNCTION BENCHMARK PURPOSES
*
*	auto start = std::time(NULL);
*	execute the time-dependent wigner function here
*	auto end = std::time(NULL);
*
*	std::cout << "Started at " << time_stamp(localtime_xp(start), "%T") << std::endl;
*	std::cout << "Finished at " << time_stamp(localtime_xp(end), "%T") << std::endl;
*
*	auto duration = difftime(end, start) / 60.0;
*	std::cout << "Duration: " << duration << " minutes" << std::endl;
*/

	// SECTION 1
//	photon_statistics_n("photon_statistics_n_jackson_25_100.dat", &f_jackson, 0.25, 1.0, -0.25);
//	photon_statistics_n("photon_statistics_n_jackson_75_200.dat", &f_jackson, 0.75, 2.0 - LDBL_EPSILON, -0.25 + 0.5 / 3);
//	photon_statistics_n("photon_statistics_n_jackson_90_300.dat", &f_jackson, 0.9, 3.0, -0.25 + 2 * 0.5 / 3);

//	photon_statistics_n("photon_statistics_n_macfarlane_25_100.dat", &f_macfarlane, 0.25, 1.0, -0.25);
//	photon_statistics_n("photon_statistics_n_macfarlane_75_150.dat", &f_macfarlane, 0.75, 1.5, -0.25 + 0.5 / 3);
//	photon_statistics_n("photon_statistics_n_macfarlane_90_225.dat", &f_macfarlane, 0.9, 2.25, -0.25 + 2 * 0.5 / 3);

//	photon_statistics_n("photon_statistics_n_canonical_600.dat", &f_jackson, 1.0 - LDBL_EPSILON, 6.0);
//	photon_statistics_n("photon_statistics_n_jackson_99_650.dat", &f_jackson, 0.99, 6.5);
//	photon_statistics_n("photon_statistics_n_macfarlane_99_500.dat", &f_macfarlane, 0.99, 5.0);


////	wigner("wigner_canonical_100.dat", &f_jackson, 1.0 - LDBL_EPSILON, 1.0);
////	wigner("wigner_canonical_250.dat", &f_jackson, 1.0 - LDBL_EPSILON, 2.5);
////	wigner("wigner_canonical_500.dat", &f_jackson, 1.0 - LDBL_EPSILON, 5.0);

////	wigner("wigner_harmonious_25.dat", &f_jackson, LDBL_EPSILON, 0.25);
////	wigner("wigner_harmonious_50.dat", &f_jackson, LDBL_EPSILON, 0.5);
////	wigner("wigner_harmonious_75.dat", &f_jackson, LDBL_EPSILON, 0.75);

////	wigner("wigner_jackson_25_100.dat", &f_jackson, 0.25, 1.0);
////	wigner("wigner_jackson_75_200.dat", &f_jackson, 0.75, 2.0 - LDBL_EPSILON);
////	wigner("wigner_jackson_90_300.dat", &f_jackson, 0.9, 3.0);

////	wigner("wigner_macfarlane_25_100.dat", &f_macfarlane, 0.25, 1.0);
////	wigner("wigner_macfarlane_75_150.dat", &f_macfarlane, 0.75, 1.5);
////	wigner("wigner_macfarlane_90_225.dat", &f_macfarlane, 0.9, 2.25);


//	heisenberg_alpha("heisenberg_alpha_canonical.dat", &f_jackson, 1.0 - LDBL_EPSILON);
//	heisenberg_alpha("heisenberg_alpha_harmonious.dat", &f_jackson, LDBL_EPSILON);

//	heisenberg_alpha("heisenberg_alpha_jackson_25.dat", &f_jackson, 0.25);
//	heisenberg_alpha("heisenberg_alpha_macfarlane_25.dat", &f_macfarlane, 0.25);

//	heisenberg_alpha("heisenberg_alpha_jackson_50.dat", &f_jackson, 0.5);
//	heisenberg_alpha("heisenberg_alpha_macfarlane_50.dat", &f_macfarlane, 0.5);

//	heisenberg_alpha("heisenberg_alpha_jackson_75.dat", &f_jackson, 0.75);
//	heisenberg_alpha("heisenberg_alpha_macfarlane_75.dat", &f_macfarlane, 0.75);

//	fidelity_superposition_alpha("fidelity_alpha_canonical_0_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.0, 0.0);
//	fidelity_superposition_alpha("fidelity_alpha_harmonious_0_even.dat", &f_jackson, LDBL_EPSILON, 0.0, 0.0);
//	fidelity_superposition_t("fidelity_t_canonical_10_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.1, 0.0);
//	fidelity_superposition_t("fidelity_t_harmonious_10_even.dat", &f_jackson, LDBL_EPSILON, 0.1, 0.0);
//	fidelity_superposition("fidelity_jackson_0_even.dat", &f_jackson, LDBL_EPSILON, 0.0);
//	fidelity_superposition("fidelity_macfarlane_0_even.dat", &f_macfarlane, LDBL_EPSILON, 0.0);

//	fidelity_superposition_alpha("fidelity_alpha_canonical_25_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.25, 0.0);
//	fidelity_superposition_alpha("fidelity_alpha_harmonious_25_even.dat", &f_jackson, LDBL_EPSILON, 0.25, 0.0);
//	fidelity_superposition_t("fidelity_t_canonical_25_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.25, 0.0);
//	fidelity_superposition_t("fidelity_t_harmonious_25_even.dat", &f_jackson, LDBL_EPSILON, 0.25, 0.0);
//	fidelity_superposition("fidelity_jackson_25_even.dat", &f_jackson, 0.25, 0.0);
//	fidelity_superposition("fidelity_macfarlane_25_even.dat", &f_macfarlane, 0.25, 0.0);

//	fidelity_superposition_alpha("fidelity_alpha_canonical_50_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.5, 0.0);
//	fidelity_superposition_alpha("fidelity_alpha_harmonious_50_even.dat", &f_jackson, LDBL_EPSILON, 0.5, 0.0);
//	fidelity_superposition_t("fidelity_t_canonical_50_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.5, 0.0);
//	fidelity_superposition_t("fidelity_t_harmonious_50_even.dat", &f_jackson, LDBL_EPSILON, 0.5, 0.0);
//	fidelity_superposition("fidelity_jackson_50_even.dat", &f_jackson, 0.5, 0.0);
//	fidelity_superposition("fidelity_macfarlane_50_even.dat", &f_macfarlane, 0.5, 0.0);

//	fidelity_superposition_alpha("fidelity_alpha_canonical_100_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 1.0, 0.0);
//	fidelity_superposition_alpha("fidelity_alpha_harmonious_100_even.dat", &f_jackson, LDBL_EPSILON, 1.0, 0.0);
//	fidelity_superposition_t("fidelity_t_canonical_90_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.9, 0.0);
//	fidelity_superposition_t("fidelity_t_harmonious_90_even.dat", &f_jackson, LDBL_EPSILON, 0.9, 0.0);
//	fidelity_superposition("fidelity_jackson_100_even.dat", &f_jackson, 1.0, 0.0);
//	fidelity_superposition("fidelity_macfarlane_100_even.dat", &f_macfarlane, 1.0, 0.0);


//	fidelity_superposition_alpha("fidelity_alpha_canonical_0_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.0, M_PI);
//	fidelity_superposition_alpha("fidelity_alpha_harmonious_0_odd.dat", &f_jackson, LDBL_EPSILON, 0.0, M_PI);
//	fidelity_superposition_t("fidelity_t_canonical_10_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.1, M_PI);
//	fidelity_superposition_t("fidelity_t_harmonious_10_odd.dat", &f_jackson, LDBL_EPSILON, 0.1, M_PI);
//	fidelity_superposition("fidelity_jackson_0_odd.dat", &f_jackson, LDBL_EPSILON, M_PI);
//	fidelity_superposition("fidelity_macfarlane_0_odd.dat", &f_macfarlane, LDBL_EPSILON, M_PI);

//	fidelity_superposition_alpha("fidelity_alpha_canonical_25_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.25, M_PI);
//	fidelity_superposition_alpha("fidelity_alpha_harmonious_25_odd.dat", &f_jackson, LDBL_EPSILON, 0.25, M_PI);
//	fidelity_superposition_t("fidelity_t_canonical_25_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.25, M_PI);
//	fidelity_superposition_t("fidelity_t_harmonious_25_odd.dat", &f_jackson, LDBL_EPSILON, 0.25, M_PI);
//	fidelity_superposition("fidelity_jackson_25_odd.dat", &f_jackson, 0.25, M_PI);
//	fidelity_superposition("fidelity_macfarlane_25_odd.dat", &f_macfarlane, 0.25, M_PI);

//	fidelity_superposition_alpha("fidelity_alpha_canonical_50_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.5, M_PI);
//	fidelity_superposition_alpha("fidelity_alpha_harmonious_50_odd.dat", &f_jackson, LDBL_EPSILON, 0.5, M_PI);
//	fidelity_superposition_t("fidelity_t_canonical_50_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.5, M_PI);
//	fidelity_superposition_t("fidelity_t_harmonious_50_odd.dat", &f_jackson, LDBL_EPSILON, 0.5, M_PI);
//	fidelity_superposition("fidelity_jackson_50_odd.dat", &f_jackson, 0.5, M_PI);
//	fidelity_superposition("fidelity_macfarlane_50_odd.dat", &f_macfarlane, 0.5, M_PI);

//	fidelity_superposition_alpha("fidelity_alpha_canonical_100_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 1.0, M_PI);
//	fidelity_superposition_alpha("fidelity_alpha_harmonious_100_odd.dat", &f_jackson, LDBL_EPSILON, 1.0, M_PI);
//	fidelity_superposition_t("fidelity_t_canonical_90_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.9, M_PI);
//	fidelity_superposition_t("fidelity_t_harmonious_90_odd.dat", &f_jackson, LDBL_EPSILON, 0.9, M_PI);
//	fidelity_superposition("fidelity_jackson_100_odd.dat", &f_jackson, 1.0, M_PI);
//	fidelity_superposition("fidelity_macfarlane_100_odd.dat", &f_macfarlane, 1.0, M_PI);


////	visibility_superposition_t("visibility_t_canonical_5.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.05, 0);
////	visibility_superposition_t("visibility_t_canonical_20.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.2, 0);
////	visibility_superposition_t("visibility_t_canonical_50.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.5, 0);
////	visibility_superposition_t("visibility_t_canonical_80.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.8, 0);

////	visibility_superposition_t("visibility_t_harmonious_5.dat", &f_jackson, LDBL_EPSILON, 0.05, 0);
////	visibility_superposition_t("visibility_t_harmonious_50.dat", &f_jackson, LDBL_EPSILON, 0.5, 0);
////	visibility_superposition_t("visibility_t_harmonious_80.dat", &f_jackson, LDBL_EPSILON, 0.8, 0);
////	visibility_superposition_t("visibility_t_harmonious_95.dat", &f_jackson, LDBL_EPSILON, 0.95, 0);

////	visibility_superposition_t("visibility_t_harmonious_5_1.dat", &f_jackson, LDBL_EPSILON, 0.05, 1);
////	visibility_superposition_t("visibility_t_harmonious_50_1.dat", &f_jackson, LDBL_EPSILON, 0.5, 1);
////	visibility_superposition_t("visibility_t_harmonious_80_1.dat", &f_jackson, LDBL_EPSILON, 0.8, 1);
////	visibility_superposition_t("visibility_t_harmonious_95_1.dat", &f_jackson, LDBL_EPSILON, 0.95, 1);

//	visibility_superposition_n("visibility_n_harmonious_5_0.dat", &f_jackson, LDBL_EPSILON, 0.05, 0.0);
//	visibility_superposition_n("visibility_n_harmonious_5_25.dat", &f_jackson, LDBL_EPSILON, 0.05, 0.25);
//	visibility_superposition_n("visibility_n_harmonious_5_50.dat", &f_jackson, LDBL_EPSILON, 0.05, 0.5);
//	visibility_superposition_n("visibility_n_harmonious_5_100.dat", &f_jackson, LDBL_EPSILON, 0.05, 1.0);

//	visibility_superposition_n("visibility_n_harmonious_20_0.dat", &f_jackson, LDBL_EPSILON, 0.2, 0.0);
//	visibility_superposition_n("visibility_n_harmonious_20_25.dat", &f_jackson, LDBL_EPSILON, 0.2, 0.25);
//	visibility_superposition_n("visibility_n_harmonious_20_50.dat", &f_jackson, LDBL_EPSILON, 0.2, 0.5);
//	visibility_superposition_n("visibility_n_harmonious_20_100.dat", &f_jackson, LDBL_EPSILON, 0.2, 1.0);

//	visibility_superposition_n("visibility_n_harmonious_50_0.dat", &f_jackson, LDBL_EPSILON, 0.5, 0.0);
//	visibility_superposition_n("visibility_n_harmonious_50_25.dat", &f_jackson, LDBL_EPSILON, 0.5, 0.25);
//	visibility_superposition_n("visibility_n_harmonious_50_50.dat", &f_jackson, LDBL_EPSILON, 0.5, 0.5);
//	visibility_superposition_n("visibility_n_harmonious_50_100.dat", &f_jackson, LDBL_EPSILON, 0.5, 1.0);

//	visibility_superposition_n("visibility_n_harmonious_80_0.dat", &f_jackson, LDBL_EPSILON, 0.8, 0.0);
//	visibility_superposition_n("visibility_n_harmonious_80_25.dat", &f_jackson, LDBL_EPSILON, 0.8, 0.25);	
//	visibility_superposition_n("visibility_n_harmonious_80_50.dat", &f_jackson, LDBL_EPSILON, 0.8, 0.5);
//	visibility_superposition_n("visibility_n_harmonious_80_100.dat", &f_jackson, LDBL_EPSILON, 0.8, 1.0);


////	visibility_superposition_t("visibility_t_jackson_25_5.dat", &f_jackson, 0.25, 0.05, 0);
////	visibility_superposition_t("visibility_t_jackson_25_50.dat", &f_jackson, 0.25, 0.5, 0);
////	visibility_superposition_t("visibility_t_jackson_25_80.dat", &f_jackson, 0.25, 0.8, 0);
////	visibility_superposition_t("visibility_t_jackson_25_95.dat", &f_jackson, 0.25, 0.95, 0);

////	visibility_superposition_t("visibility_t_jackson_75_5.dat", &f_jackson, 0.75, 0.05, 0);
////	visibility_superposition_t("visibility_t_jackson_75_50.dat", &f_jackson, 0.75, 0.5, 0);
////	visibility_superposition_t("visibility_t_jackson_75_80.dat", &f_jackson, 0.75, 0.8, 0);
////	visibility_superposition_t("visibility_t_jackson_75_95.dat", &f_jackson, 0.75, 0.95, 0);

//	visibility_superposition_n("visibility_n_jackson_25_5_0.dat", &f_jackson, 0.25, 0.05, 0.0);
//	visibility_superposition_n("visibility_n_jackson_25_5_25.dat", &f_jackson, 0.25, 0.05, 0.25);
//	visibility_superposition_n("visibility_n_jackson_25_5_50.dat", &f_jackson, 0.25, 0.05, 0.5);
//	visibility_superposition_n("visibility_n_jackson_25_5_100.dat", &f_jackson, 0.25, 0.05, 1.0);

//	visibility_superposition_n("visibility_n_jackson_25_20_0.dat", &f_jackson, 0.25, 0.2, 0.0);
//	visibility_superposition_n("visibility_n_jackson_25_20_25.dat", &f_jackson, 0.25, 0.2, 0.25);
//	visibility_superposition_n("visibility_n_jackson_25_20_50.dat", &f_jackson, 0.25, 0.2, 0.5);
//	visibility_superposition_n("visibility_n_jackson_25_20_100.dat", &f_jackson, 0.25, 0.2, 1.0);

//	visibility_superposition_n("visibility_n_jackson_25_50_0.dat", &f_jackson, 0.25, 0.5, 0.0);
//	visibility_superposition_n("visibility_n_jackson_25_50_25.dat", &f_jackson, 0.25, 0.5, 0.25);
//	visibility_superposition_n("visibility_n_jackson_25_50_50.dat", &f_jackson, 0.25, 0.5, 0.5);
//	visibility_superposition_n("visibility_n_jackson_25_50_100.dat", &f_jackson, 0.25, 0.5, 1.0);

//	visibility_superposition_n("visibility_n_jackson_25_80_0.dat", &f_jackson, 0.25, 0.8, 0.0);
//	visibility_superposition_n("visibility_n_jackson_25_80_25.dat", &f_jackson, 0.25, 0.8, 0.25);
//	visibility_superposition_n("visibility_n_jackson_25_80_50.dat", &f_jackson, 0.25, 0.8, 0.5);
//	visibility_superposition_n("visibility_n_jackson_25_80_100.dat", &f_jackson, 0.25, 0.8, 1.0);

//	visibility_superposition_n("visibility_n_jackson_50_5_0.dat", &f_jackson, 0.5, 0.05, 0.0);
//	visibility_superposition_n("visibility_n_jackson_50_5_25.dat", &f_jackson, 0.5, 0.05, 0.25);
//	visibility_superposition_n("visibility_n_jackson_50_5_50.dat", &f_jackson, 0.5, 0.05, 0.5);
//	visibility_superposition_n("visibility_n_jackson_50_5_100.dat", &f_jackson, 0.5, 0.05, 1.0);

//	visibility_superposition_n("visibility_n_jackson_50_20_0.dat", &f_jackson, 0.5, 0.2, 0.0);
//	visibility_superposition_n("visibility_n_jackson_50_20_25.dat", &f_jackson, 0.5, 0.2, 0.25);
//	visibility_superposition_n("visibility_n_jackson_50_20_50.dat", &f_jackson, 0.5, 0.2, 0.5);
//	visibility_superposition_n("visibility_n_jackson_50_20_100.dat", &f_jackson, 0.5, 0.2, 1.0);

//	visibility_superposition_n("visibility_n_jackson_50_50_0.dat", &f_jackson, 0.5, 0.5, 0.0);
//	visibility_superposition_n("visibility_n_jackson_50_50_25.dat", &f_jackson, 0.5, 0.5, 0.25);
//	visibility_superposition_n("visibility_n_jackson_50_50_50.dat", &f_jackson, 0.5, 0.5, 0.5);
//	visibility_superposition_n("visibility_n_jackson_50_50_100.dat", &f_jackson, 0.5, 0.5, 1.0);

//	visibility_superposition_n("visibility_n_jackson_50_80_0.dat", &f_jackson, 0.5, 0.8, 0.0);
//	visibility_superposition_n("visibility_n_jackson_50_80_25.dat", &f_jackson, 0.5, 0.8, 0.25);
//	visibility_superposition_n("visibility_n_jackson_50_80_50.dat", &f_jackson, 0.5, 0.8, 0.5);
//	visibility_superposition_n("visibility_n_jackson_50_80_100.dat", &f_jackson, 0.5, 0.8, 1.0);

//	visibility_superposition_n("visibility_n_jackson_75_5_0.dat", &f_jackson, 0.75, 0.05, 0.0);
//	visibility_superposition_n("visibility_n_jackson_75_5_25.dat", &f_jackson, 0.75, 0.05, 0.25);
//	visibility_superposition_n("visibility_n_jackson_75_5_50.dat", &f_jackson, 0.75, 0.05, 0.5);
//	visibility_superposition_n("visibility_n_jackson_75_5_100.dat", &f_jackson, 0.75, 0.05, 1.0);

//	visibility_superposition_n("visibility_n_jackson_75_20_0.dat", &f_jackson, 0.75, 0.2, 0.0);
//	visibility_superposition_n("visibility_n_jackson_75_20_25.dat", &f_jackson, 0.75, 0.2, 0.25);
//	visibility_superposition_n("visibility_n_jackson_75_20_50.dat", &f_jackson, 0.75, 0.2, 0.5);
//	visibility_superposition_n("visibility_n_jackson_75_20_100.dat", &f_jackson, 0.75, 0.2, 1.0);

//	visibility_superposition_n("visibility_n_jackson_75_50_0.dat", &f_jackson, 0.75, 0.5, 0.0);
//	visibility_superposition_n("visibility_n_jackson_75_50_25.dat", &f_jackson, 0.75, 0.5, 0.25);
//	visibility_superposition_n("visibility_n_jackson_75_50_50.dat", &f_jackson, 0.75, 0.5, 0.5);
//	visibility_superposition_n("visibility_n_jackson_75_50_100.dat", &f_jackson, 0.75, 0.5, 1.0);

//	visibility_superposition_n("visibility_n_jackson_75_80_0.dat", &f_jackson, 0.75, 0.8, 0.0);
//	visibility_superposition_n("visibility_n_jackson_75_80_25.dat", &f_jackson, 0.75, 0.8, 0.25);
//	visibility_superposition_n("visibility_n_jackson_75_80_50.dat", &f_jackson, 0.75, 0.8, 0.5);
//	visibility_superposition_n("visibility_n_jackson_75_80_100.dat", &f_jackson, 0.75, 0.8, 1.0);


////	visibility_superposition_t("visibility_t_macfarlane_25_5.dat", &f_macfarlane, 0.25, 0.05, 0);
////	visibility_superposition_t("visibility_t_macfarlane_25_50.dat", &f_macfarlane, 0.25, 0.5, 0);
////	visibility_superposition_t("visibility_t_macfarlane_25_80.dat", &f_macfarlane, 0.25, 0.8, 0);
////	visibility_superposition_t("visibility_t_macfarlane_25_95.dat", &f_macfarlane, 0.25, 0.95, 0);

////	visibility_superposition_t("visibility_t_macfarlane_75_5.dat", &f_macfarlane, 0.75, 0.05, 0);
////	visibility_superposition_t("visibility_t_macfarlane_75_50.dat", &f_macfarlane, 0.75, 0.5, 0);
////	visibility_superposition_t("visibility_t_macfarlane_75_80.dat", &f_macfarlane, 0.75, 0.8, 0);
////	visibility_superposition_t("visibility_t_macfarlane_75_95.dat", &f_macfarlane, 0.75, 0.95, 0);

//	visibility_superposition_n("visibility_n_macfarlane_25_5_0.dat", &f_macfarlane, 0.25, 0.05, 0.0);
//	visibility_superposition_n("visibility_n_macfarlane_25_5_25.dat", &f_macfarlane, 0.25, 0.05, 0.25);
//	visibility_superposition_n("visibility_n_macfarlane_25_5_50.dat", &f_macfarlane, 0.25, 0.05, 0.5);
//	visibility_superposition_n("visibility_n_macfarlane_25_5_100.dat", &f_macfarlane, 0.25, 0.05, 1.0);

//	visibility_superposition_n("visibility_n_macfarlane_25_20_0.dat", &f_macfarlane, 0.25, 0.2, 0.0);
//	visibility_superposition_n("visibility_n_macfarlane_25_20_25.dat", &f_macfarlane, 0.25, 0.2, 0.25);
//	visibility_superposition_n("visibility_n_macfarlane_25_20_50.dat", &f_macfarlane, 0.25, 0.2, 0.5);
//	visibility_superposition_n("visibility_n_macfarlane_25_20_100.dat", &f_macfarlane, 0.25, 0.2, 1.0);

//	visibility_superposition_n("visibility_n_macfarlane_25_50_0.dat", &f_macfarlane, 0.25, 0.5, 0.0);
//	visibility_superposition_n("visibility_n_macfarlane_25_50_25.dat", &f_macfarlane, 0.25, 0.5, 0.25);
//	visibility_superposition_n("visibility_n_macfarlane_25_50_50.dat", &f_macfarlane, 0.25, 0.5, 0.5);
//	visibility_superposition_n("visibility_n_macfarlane_25_50_100.dat", &f_macfarlane, 0.25, 0.5, 1.0);

//	visibility_superposition_n("visibility_n_macfarlane_25_80_0.dat", &f_macfarlane, 0.25, 0.8, 0.0);
//	visibility_superposition_n("visibility_n_macfarlane_25_80_25.dat", &f_macfarlane, 0.25, 0.8, 0.25);
//	visibility_superposition_n("visibility_n_macfarlane_25_80_50.dat", &f_macfarlane, 0.25, 0.8, 0.5);
//	visibility_superposition_n("visibility_n_macfarlane_25_80_100.dat", &f_macfarlane, 0.25, 0.8, 1.0);

//	visibility_superposition_n("visibility_n_macfarlane_50_5_0.dat", &f_macfarlane, 0.5, 0.05, 0.0);
//	visibility_superposition_n("visibility_n_macfarlane_50_5_25.dat", &f_macfarlane, 0.5, 0.05, 0.25);
//	visibility_superposition_n("visibility_n_macfarlane_50_5_50.dat", &f_macfarlane, 0.5, 0.05, 0.5);
//	visibility_superposition_n("visibility_n_macfarlane_50_5_100.dat", &f_macfarlane, 0.5, 0.05, 1.0);

//	visibility_superposition_n("visibility_n_macfarlane_50_20_0.dat", &f_macfarlane, 0.5, 0.2, 0.0);
//	visibility_superposition_n("visibility_n_macfarlane_50_20_25.dat", &f_macfarlane, 0.5, 0.2, 0.25);
//	visibility_superposition_n("visibility_n_macfarlane_50_20_50.dat", &f_macfarlane, 0.5, 0.2, 0.5);
//	visibility_superposition_n("visibility_n_macfarlane_50_20_100.dat", &f_macfarlane, 0.5, 0.2, 1.0);

//	visibility_superposition_n("visibility_n_macfarlane_50_50_0.dat", &f_macfarlane, 0.5, 0.5, 0.0);
//	visibility_superposition_n("visibility_n_macfarlane_50_50_25.dat", &f_macfarlane, 0.5, 0.5, 0.25);
//	visibility_superposition_n("visibility_n_macfarlane_50_50_50.dat", &f_macfarlane, 0.5, 0.5, 0.5);
//	visibility_superposition_n("visibility_n_macfarlane_50_50_100.dat", &f_macfarlane, 0.5, 0.5, 1.0);

//	visibility_superposition_n("visibility_n_macfarlane_50_80_0.dat", &f_macfarlane, 0.5, 0.8, 0.0);
//	visibility_superposition_n("visibility_n_macfarlane_50_80_25.dat", &f_macfarlane, 0.5, 0.8, 0.25);
//	visibility_superposition_n("visibility_n_macfarlane_50_80_50.dat", &f_macfarlane, 0.5, 0.8, 0.5);
//	visibility_superposition_n("visibility_n_macfarlane_50_80_100.dat", &f_macfarlane, 0.5, 0.8, 1.0);

//	visibility_superposition_n("visibility_n_macfarlane_75_5_0.dat", &f_macfarlane, 0.75, 0.05, 0.0);
//	visibility_superposition_n("visibility_n_macfarlane_75_5_25.dat", &f_macfarlane, 0.75, 0.05, 0.25);
//	visibility_superposition_n("visibility_n_macfarlane_75_5_50.dat", &f_macfarlane, 0.75, 0.05, 0.5);
//	visibility_superposition_n("visibility_n_macfarlane_75_5_100.dat", &f_macfarlane, 0.75, 0.05, 1.0);

//	visibility_superposition_n("visibility_n_macfarlane_75_20_0.dat", &f_macfarlane, 0.75, 0.2, 0.0);
//	visibility_superposition_n("visibility_n_macfarlane_75_20_25.dat", &f_macfarlane, 0.75, 0.2, 0.25);
//	visibility_superposition_n("visibility_n_macfarlane_75_20_50.dat", &f_macfarlane, 0.75, 0.2, 0.5);
//	visibility_superposition_n("visibility_n_macfarlane_75_20_100.dat", &f_macfarlane, 0.75, 0.2, 1.0);

//	visibility_superposition_n("visibility_n_macfarlane_75_50_0.dat", &f_macfarlane, 0.75, 0.5, 0.0);
//	visibility_superposition_n("visibility_n_macfarlane_75_50_25.dat", &f_macfarlane, 0.75, 0.5, 0.25);
//	visibility_superposition_n("visibility_n_macfarlane_75_50_50.dat", &f_macfarlane, 0.75, 0.5, 0.5);
//	visibility_superposition_n("visibility_n_macfarlane_75_50_100.dat", &f_macfarlane, 0.75, 0.5, 1.0);

//	visibility_superposition_n("visibility_n_macfarlane_75_80_0.dat", &f_macfarlane, 0.75, 0.8, 0.0);
//	visibility_superposition_n("visibility_n_macfarlane_75_80_25.dat", &f_macfarlane, 0.75, 0.8, 0.25);
//	visibility_superposition_n("visibility_n_macfarlane_75_80_50.dat", &f_macfarlane, 0.75, 0.8, 0.5);
//	visibility_superposition_n("visibility_n_macfarlane_75_80_100.dat", &f_macfarlane, 0.75, 0.8, 1.0);


////	photon_statistics_superposition_n("photon_statistics_n_canonical_300_0_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 3.0, 0.0, 0.0, -0.25);
////	photon_statistics_superposition_n("photon_statistics_n_canonical_450_0_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 4.5, 0.0, 0.0, 0.25);

////	photon_statistics_superposition_n("photon_statistics_n_canonical_300_0_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 3.0, 0.0, M_PI, -0.25);
////	photon_statistics_superposition_n("photon_statistics_n_canonical_450_0_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 4.5, 0.0, M_PI, 0.25);

//	photon_statistics_superposition_n("photon_statistics_n_canonical_5_0_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.05, 0.0, 0.0);
//	photon_statistics_superposition_n("photon_statistics_n_canonical_80_0_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.8, 0.0, 0.0);

//	photon_statistics_superposition_n("photon_statistics_n_canonical_5_0_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.05, 0.0, M_PI);
//	photon_statistics_superposition_n("photon_statistics_n_canonical_80_0_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.8, 0.0, M_PI);

////	photon_statistics_superposition_n("photon_statistics_n_harmonious_75_0_even.dat", &f_jackson, LDBL_EPSILON, 0.75, 0.0, 0.0, -0.25);
////	photon_statistics_superposition_n("photon_statistics_n_harmonious_95_0_even.dat", &f_jackson, LDBL_EPSILON, 0.95, 0.0, 0.0, 0.25);

////	photon_statistics_superposition_n("photon_statistics_n_harmonious_75_0_odd.dat", &f_jackson, LDBL_EPSILON, 0.75, 0.0, M_PI, -0.25);
////	photon_statistics_superposition_n("photon_statistics_n_harmonious_95_0_odd.dat", &f_jackson, LDBL_EPSILON, 0.95, 0.0, M_PI, 0.25);

//	photon_statistics_superposition_n("photon_statistics_n_harmonious_5_0_even.dat", &f_jackson, LDBL_EPSILON, 0.05, 0.0, 0.0);
//	photon_statistics_superposition_n("photon_statistics_n_harmonious_80_0_even.dat", &f_jackson, LDBL_EPSILON, 0.8, 0.0, 0.0);

//	photon_statistics_superposition_n("photon_statistics_n_harmonious_5_0_odd.dat", &f_jackson, LDBL_EPSILON, 0.05, 0.0, M_PI);	
//	photon_statistics_superposition_n("photon_statistics_n_harmonious_80_0_odd.dat", &f_jackson, LDBL_EPSILON, 0.8, 0.0, M_PI);


//	photon_statistics_superposition_n("photon_statistics_n_jackson_25_5_0_even.dat", &f_jackson, 0.25, 0.05, 0.0, 0.0);
//	photon_statistics_superposition_n("photon_statistics_n_jackson_25_80_0_even.dat", &f_jackson, 0.25, 0.8, 0.0, 0.0);

//	photon_statistics_superposition_n("photon_statistics_n_jackson_25_5_0_odd.dat", &f_jackson, 0.25, 0.05, 0.0, M_PI);
//	photon_statistics_superposition_n("photon_statistics_n_jackson_25_80_0_odd.dat", &f_jackson, 0.25, 0.8, 0.0, M_PI);

////	photon_statistics_superposition_n("photon_statistics_n_jackson_60_100_0_even.dat", &f_jackson, 0.6, 1.0, 0.0, 0.0, -0.25);
////	photon_statistics_superposition_n("photon_statistics_n_jackson_60_150_0_even.dat", &f_jackson, 0.6, 1.5, 0.0, 0.0, 0.25);

////	photon_statistics_superposition_n("photon_statistics_n_jackson_60_100_0_odd.dat", &f_jackson, 0.6, 1.0, 0.0, M_PI, -0.25);
////	photon_statistics_superposition_n("photon_statistics_n_jackson_60_150_0_odd.dat", &f_jackson, 0.6, 1.5, 0.0, M_PI, 0.25);

//	photon_statistics_superposition_n("photon_statistics_n_jackson_50_5_0_even.dat", &f_jackson, 0.5, 0.05, 0.0, 0.0);
//	photon_statistics_superposition_n("photon_statistics_n_jackson_50_80_0_even.dat", &f_jackson, 0.5, 0.8, 0.0, 0.0);

//	photon_statistics_superposition_n("photon_statistics_n_jackson_50_5_0_odd.dat", &f_jackson, 0.5, 0.05, 0.0, M_PI);
//	photon_statistics_superposition_n("photon_statistics_n_jackson_50_80_0_odd.dat", &f_jackson, 0.5, 0.8, 0.0, M_PI);

////	photon_statistics_superposition_n("photon_statistics_n_jackson_85_150_0_even.dat", &f_jackson, 0.85, 1.5, 0.0, 0.0, -0.25);
////	photon_statistics_superposition_n("photon_statistics_n_jackson_85_250_0_even.dat", &f_jackson, 0.85, 2.5, 0.0, 0.0, 0.25);

////	photon_statistics_superposition_n("photon_statistics_n_jackson_85_150_0_odd.dat", &f_jackson, 0.85, 1.5, 0.0, M_PI, -0.25);
////	photon_statistics_superposition_n("photon_statistics_n_jackson_85_250_0_odd.dat", &f_jackson, 0.85, 2.5, 0.0, M_PI, 0.25);

//	photon_statistics_superposition_n("photon_statistics_n_jackson_75_5_0_even.dat", &f_jackson, 0.75, 0.05, 0.0, 0.0);
//	photon_statistics_superposition_n("photon_statistics_n_jackson_75_80_0_even.dat", &f_jackson, 0.75, 0.8, 0.0, 0.0);

//	photon_statistics_superposition_n("photon_statistics_n_jackson_75_5_0_odd.dat", &f_jackson, 0.75, 0.05, 0.0, M_PI);
//	photon_statistics_superposition_n("photon_statistics_n_jackson_75_80_0_odd.dat", &f_jackson, 0.75, 0.8, 0.0, M_PI);


//	photon_statistics_superposition_n("photon_statistics_n_macfarlane_25_5_0_even.dat", &f_macfarlane, 0.25, 0.05, 0.0, 0.0);
//	photon_statistics_superposition_n("photon_statistics_n_macfarlane_25_80_0_even.dat", &f_macfarlane, 0.25, 0.8, 0.0, 0.0);

//	photon_statistics_superposition_n("photon_statistics_n_macfarlane_25_5_0_odd.dat", &f_macfarlane, 0.25, 0.05, 0.0, M_PI);
//	photon_statistics_superposition_n("photon_statistics_n_macfarlane_25_80_0_odd.dat", &f_macfarlane, 0.25, 0.8, 0.0, M_PI);

////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_60_100_0_even.dat", &f_macfarlane, 0.6, 1.0, 0.0, 0.0, -0.25);
////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_60_125_0_even.dat", &f_macfarlane, 0.6, 1.25 - LDBL_EPSILON, 0.0, 0.0, 0.25);

////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_60_100_0_odd.dat", &f_macfarlane, 0.6, 1.0, 0.0, M_PI, -0.25);
////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_60_125_0_odd.dat", &f_macfarlane, 0.6, 1.25 - LDBL_EPSILON, 0.0, M_PI, 0.25);

//	photon_statistics_superposition_n("photon_statistics_n_macfarlane_50_5_0_even.dat", &f_macfarlane, 0.5, 0.05, 0.0, 0.0);
//	photon_statistics_superposition_n("photon_statistics_n_macfarlane_50_80_0_even.dat", &f_macfarlane, 0.5, 0.8, 0.0, 0.0);

//	photon_statistics_superposition_n("photon_statistics_n_macfarlane_50_5_0_odd.dat", &f_macfarlane, 0.5, 0.05, 0.0, M_PI);
//	photon_statistics_superposition_n("photon_statistics_n_macfarlane_50_80_0_odd.dat", &f_macfarlane, 0.5, 0.8, 0.0, M_PI);


////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_85_125_0_even.dat", &f_macfarlane, 0.85, 1.25, 0.0, 0.0, -0.25);
////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_85_180_0_even.dat", &f_macfarlane, 0.85, 1.8, 0.0, 0.0, 0.25);

////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_85_125_0_odd.dat", &f_macfarlane, 0.85, 1.25, 0.0, M_PI, -0.25);
////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_85_180_0_odd.dat", &f_macfarlane, 0.85, 1.8, 0.0, M_PI, 0.25);

//	photon_statistics_superposition_n("photon_statistics_n_macfarlane_75_5_0_even.dat", &f_macfarlane, 0.75, 0.05, 0.0, 0.0);
//	photon_statistics_superposition_n("photon_statistics_n_macfarlane_75_80_0_even.dat", &f_macfarlane, 0.75, 0.8, 0.0, 0.0);

//	photon_statistics_superposition_n("photon_statistics_n_macfarlane_75_5_0_odd.dat", &f_macfarlane, 0.75, 0.05, 0.0, M_PI);
//	photon_statistics_superposition_n("photon_statistics_n_macfarlane_75_80_0_odd.dat", &f_macfarlane, 0.75, 0.8, 0.0, M_PI);


////	photon_statistics_superposition("photon_statistics_canonical_300_0.dat", &f_jackson, 1.0 - LDBL_EPSILON, 3.0, 0.0);
////	photon_statistics_superposition("photon_statistics_canonical_450_0.dat", &f_jackson, 1.0 - LDBL_EPSILON, 4.5, 0.0);

////	photon_statistics_superposition("photon_statistics_harmonious_75_0.dat", &f_jackson, LDBL_EPSILON, 0.75, 0.0);
////	photon_statistics_superposition("photon_statistics_harmonious_95_0.dat", &f_jackson, LDBL_EPSILON, 0.95, 0.0);

////	photon_statistics_superposition("photon_statistics_jackson_60_100_0.dat", &f_jackson, 0.6, 1.0, 0.0);
////	photon_statistics_superposition("photon_statistics_jackson_60_150_0.dat", &f_jackson, 0.6, 1.5, 0.0);

////	photon_statistics_superposition("photon_statistics_jackson_85_150_0.dat", &f_jackson, 0.85, 1.5, 0.0);
////	photon_statistics_superposition("photon_statistics_jackson_85_250_0.dat", &f_jackson, 0.85, 2.5, 0.0);

////	photon_statistics_superposition("photon_statistics_macfarlane_60_100_0.dat", &f_macfarlane, 0.6, 1.0, 0.0);
////	photon_statistics_superposition("photon_statistics_macfarlane_60_125_0.dat", &f_macfarlane, 0.6, 1.25 - LDBL_EPSILON, 0.0);

////	photon_statistics_superposition("photon_statistics_macfarlane_85_125_0.dat", &f_macfarlane, 0.85, 1.25, 0.0);
////	photon_statistics_superposition("photon_statistics_macfarlane_85_180_0.dat", &f_macfarlane, 0.85, 1.8, 0.0);


//	mandel_parameter_superposition_t("mandel_parameter_t_canonical_5_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.05, 0.0);
//	mandel_parameter_superposition_t("mandel_parameter_t_canonical_5_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.05, M_PI);
//	mandel_parameter_superposition_t("mandel_parameter_t_canonical_20_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.2, 0.0);
//	mandel_parameter_superposition_t("mandel_parameter_t_canonical_20_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.2, M_PI);
//	mandel_parameter_superposition_t("mandel_parameter_t_canonical_50_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.5, 0.0);
//	mandel_parameter_superposition_t("mandel_parameter_t_canonical_50_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.5, M_PI);
//	mandel_parameter_superposition_t("mandel_parameter_t_canonical_80_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.8, 0.0);
//	mandel_parameter_superposition_t("mandel_parameter_t_canonical_80_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.8, M_PI);

//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_canonical_0_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_canonical_0_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.0, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_canonical_25_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.25, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_canonical_25_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.25, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_canonical_50_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.5, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_canonical_50_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.5, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_canonical_100_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 1.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_canonical_100_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 1.0, M_PI);

//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_harmonious_0_even.dat", &f_jackson, LDBL_EPSILON, 0.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_harmonious_0_odd.dat", &f_jackson, LDBL_EPSILON, 0.0, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_harmonious_25_even.dat", &f_jackson, LDBL_EPSILON, 0.25, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_harmonious_25_odd.dat", &f_jackson, LDBL_EPSILON, 0.25, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_harmonious_50_even.dat", &f_jackson, LDBL_EPSILON, 0.5, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_harmonious_50_odd.dat", &f_jackson, LDBL_EPSILON, 0.5, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_harmonious_100_even.dat", &f_jackson, LDBL_EPSILON, 1.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_harmonious_100_odd.dat", &f_jackson, LDBL_EPSILON, 1.0, M_PI);


//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_25_0_even.dat", &f_jackson, 0.25, 0.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_25_0_odd.dat", &f_jackson, 0.25, 0.0, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_25_25_even.dat", &f_jackson, 0.25, 0.25, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_25_25_odd.dat", &f_jackson, 0.25, 0.25, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_25_50_even.dat", &f_jackson, 0.25, 0.5, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_25_50_odd.dat", &f_jackson, 0.25, 0.5, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_25_100_even.dat", &f_jackson, 0.25, 1.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_25_100_odd.dat", &f_jackson, 0.25, 1.0, M_PI);

//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_50_0_even.dat", &f_jackson, 0.5, 0.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_50_0_odd.dat", &f_jackson, 0.5, 0.0, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_50_25_even.dat", &f_jackson, 0.5, 0.25, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_50_25_odd.dat", &f_jackson, 0.5, 0.25, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_50_50_even.dat", &f_jackson, 0.5, 0.5, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_50_50_odd.dat", &f_jackson, 0.5, 0.5, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_50_100_even.dat", &f_jackson, 0.5, 1.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_50_100_odd.dat", &f_jackson, 0.5, 1.0, M_PI);

//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_75_0_even.dat", &f_jackson, 0.75, 0.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_75_0_odd.dat", &f_jackson, 0.75, 0.0, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_75_25_even.dat", &f_jackson, 0.75, 0.25, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_75_25_odd.dat", &f_jackson, 0.75, 0.25, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_75_50_even.dat", &f_jackson, 0.75, 0.5, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_75_50_odd.dat", &f_jackson, 0.75, 0.5, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_75_100_even.dat", &f_jackson, 0.75, 1.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_jackson_75_100_odd.dat", &f_jackson, 0.75, 1.0, M_PI);


//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_25_0_even.dat", &f_macfarlane, 0.25, 0.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_25_0_odd.dat", &f_macfarlane, 0.25, 0.0, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_25_25_even.dat", &f_macfarlane, 0.25, 0.25, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_25_25_odd.dat", &f_macfarlane, 0.25, 0.25, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_25_50_even.dat", &f_macfarlane, 0.25, 0.5, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_25_50_odd.dat", &f_macfarlane, 0.25, 0.5, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_25_100_even.dat", &f_macfarlane, 0.25, 1.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_25_100_odd.dat", &f_macfarlane, 0.25, 1.0, M_PI);

//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_50_0_even.dat", &f_macfarlane, 0.5, 0.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_50_0_odd.dat", &f_macfarlane, 0.5, 0.0, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_50_25_even.dat", &f_macfarlane, 0.5, 0.25, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_50_25_odd.dat", &f_macfarlane, 0.5, 0.25, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_50_50_even.dat", &f_macfarlane, 0.5, 0.5, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_50_50_odd.dat", &f_macfarlane, 0.5, 0.5, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_50_100_even.dat", &f_macfarlane, 0.5, 1.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_50_100_odd.dat", &f_macfarlane, 0.5, 1.0, M_PI);

//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_75_0_even.dat", &f_macfarlane, 0.75, 0.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_75_0_odd.dat", &f_macfarlane, 0.75, 0.0, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_75_25_even.dat", &f_macfarlane, 0.75, 0.25, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_75_25_odd.dat", &f_macfarlane, 0.75, 0.25, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_75_50_even.dat", &f_macfarlane, 0.75, 0.5, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_75_50_odd.dat", &f_macfarlane, 0.75, 0.5, M_PI);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_75_100_even.dat", &f_macfarlane, 0.75, 1.0, 0.0);
//	mandel_parameter_superposition_alpha("mandel_parameter_alpha_macfarlane_75_100_odd.dat", &f_macfarlane, 0.75, 1.0, M_PI);



//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_canonical_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 0.0);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_canonical_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, M_PI);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_canonical_90.dat", &f_jackson, 1.0 - LDBL_EPSILON, M_PI / 2);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_canonical_60.dat", &f_jackson, 1.0 - LDBL_EPSILON, M_PI / 3);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_canonical_45.dat", &f_jackson, 1.0 - LDBL_EPSILON, M_PI / 4);

//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_harmonious_even.dat", &f_jackson, LDBL_EPSILON, 0.0);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_harmonious_odd.dat", &f_jackson, LDBL_EPSILON, M_PI);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_harmonious_90.dat", &f_jackson, LDBL_EPSILON, M_PI / 2);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_harmonious_60.dat", &f_jackson, LDBL_EPSILON, M_PI / 3);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_harmonious_45.dat", &f_jackson, LDBL_EPSILON, M_PI / 4);


//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_jackson_25_even.dat", &f_jackson, 0.25, 0.0);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_jackson_25_odd.dat", &f_jackson, 0.25, M_PI);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_jackson_25_90.dat", &f_jackson, 0.25, M_PI / 2);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_jackson_25_60.dat", &f_jackson, 0.25, M_PI / 3);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_jackson_25_45.dat", &f_jackson, 0.25, M_PI / 4);

//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_jackson_50_even.dat", &f_jackson, 0.5, 0.0);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_jackson_50_odd.dat", &f_jackson, 0.5, M_PI);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_jackson_50_90.dat", &f_jackson, 0.5, M_PI / 2);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_jackson_50_60.dat", &f_jackson, 0.5, M_PI / 3);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_jackson_50_45.dat", &f_jackson, 0.5, M_PI / 4);

//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_jackson_75_even.dat", &f_jackson, 0.75, 0.0);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_jackson_75_odd.dat", &f_jackson, 0.75, M_PI);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_jackson_75_90.dat", &f_jackson, 0.75, M_PI / 2);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_jackson_75_60.dat", &f_jackson, 0.75, M_PI / 3);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_jackson_75_45.dat", &f_jackson, 0.75, M_PI / 4);


//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_macfarlane_25_even.dat", &f_macfarlane, 0.25, 0.0);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_macfarlane_25_odd.dat", &f_macfarlane, 0.25, M_PI);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_macfarlane_25_90.dat", &f_macfarlane, 0.25, M_PI / 2);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_macfarlane_25_60.dat", &f_macfarlane, 0.25, M_PI / 3);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_macfarlane_25_45.dat", &f_macfarlane, 0.25, M_PI / 4);

//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_macfarlane_50_even.dat", &f_macfarlane, 0.5, 0.0);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_macfarlane_50_odd.dat", &f_macfarlane, 0.5, M_PI);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_macfarlane_50_90.dat", &f_macfarlane, 0.5, M_PI / 2);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_macfarlane_50_60.dat", &f_macfarlane, 0.5, M_PI / 3);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_macfarlane_50_45.dat", &f_macfarlane, 0.5, M_PI / 4);

//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_macfarlane_75_even.dat", &f_macfarlane, 0.75, 0.0);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_macfarlane_75_odd.dat", &f_macfarlane, 0.75, M_PI);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_macfarlane_75_90.dat", &f_macfarlane, 0.75, M_PI / 2);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_macfarlane_75_60.dat", &f_macfarlane, 0.75, M_PI / 3);
//	second_order_correlation_parameter_superposition_alpha("second_order_correlation_parameter_alpha_macfarlane_75_45.dat", &f_macfarlane, 0.75, M_PI / 4);



////	wigner_superposition("wigner_canonical_100_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 1.0, 0.0);
////	wigner_superposition("wigner_canonical_100_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 1.0, M_PI);

////	wigner_superposition("wigner_canonical_250_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 2.5, 0.0);
////	wigner_superposition("wigner_canonical_250_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 2.5, M_PI);

////	wigner_superposition("wigner_canonical_500_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 5.0, 0.0);
////	wigner_superposition("wigner_canonical_500_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 5.0, M_PI);


////	wigner_superposition("wigner_harmonious_25_even.dat", &f_jackson, LDBL_EPSILON, 0.25, 0.0);
////	wigner_superposition("wigner_harmonious_25_odd.dat", &f_jackson, LDBL_EPSILON, 0.25, M_PI);

////	wigner_superposition("wigner_harmonious_50_even.dat", &f_jackson, LDBL_EPSILON, 0.5, 0.0);
////	wigner_superposition("wigner_harmonious_50_odd.dat", &f_jackson, LDBL_EPSILON, 0.5, M_PI);

////	wigner_superposition("wigner_harmonious_75_even.dat", &f_jackson, LDBL_EPSILON, 0.75, 0.0);
////	wigner_superposition("wigner_harmonious_75_odd.dat", &f_jackson, LDBL_EPSILON, 0.75, M_PI);

////	wigner_superposition("wigner_jackson_25_100_even.dat", &f_jackson, 0.25, 1.0, 0.0);
////	wigner_superposition("wigner_jackson_25_100_odd.dat", &f_jackson, 0.25, 1.0, M_PI);

////	wigner_superposition("wigner_jackson_75_200_even.dat", &f_jackson, 0.75, 2.0 - LDBL_EPSILON, 0.0);
////	wigner_superposition("wigner_jackson_75_200_odd.dat", &f_jackson, 0.75, 2.0 - LDBL_EPSILON, M_PI);

////	wigner_superposition("wigner_jackson_90_300_even.dat", &f_jackson, 0.9, 3.0, 0.0);
////	wigner_superposition("wigner_jackson_90_300_odd.dat", &f_jackson, 0.9, 3.0, M_PI);


////	wigner_superposition("wigner_macfarlane_25_100_even.dat", &f_macfarlane, 0.25, 1.0, 0.0);
////	wigner_superposition("wigner_macfarlane_25_100_odd.dat", &f_macfarlane, 0.25, 1.0, M_PI);

////	wigner_superposition("wigner_macfarlane_75_150_even.dat", &f_macfarlane, 0.75, 1.5, 0.0);
////	wigner_superposition("wigner_macfarlane_75_150_odd.dat", &f_macfarlane, 0.75, 1.5, M_PI);

////	wigner_superposition("wigner_macfarlane_90_225_even.dat", &f_macfarlane, 0.9, 2.25, 0.0);
////	wigner_superposition("wigner_macfarlane_90_225_odd.dat", &f_macfarlane, 0.9, 2.25, M_PI);



	// SECTION 2

////	photon_statistics_superposition_n("photon_statistics_n_canonical_300_0_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 3.0, 0.0, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_canonical_300_0_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 3.0, 0.0, M_PI, 0.15);

////	photon_statistics_superposition_n("photon_statistics_n_canonical_300_10_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 3.0, 0.1, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_canonical_300_10_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 3.0, 0.1, M_PI, 0.15);

////	photon_statistics_superposition_n("photon_statistics_n_canonical_300_50_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 3.0, 0.5, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_canonical_300_50_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 3.0, 0.5, M_PI, 0.15);

////	photon_statistics_superposition_n("photon_statistics_n_canonical_300_100_even.dat", &f_jackson, 1.0 - LDBL_EPSILON, 3.0, 1.0, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_canonical_300_100_odd.dat", &f_jackson, 1.0 - LDBL_EPSILON, 3.0, 1.0, M_PI, 0.15);


////	photon_statistics_superposition_n("photon_statistics_n_harmonious_95_0_even.dat", &f_jackson, LDBL_EPSILON, 0.95, 0.0, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_harmonious_95_0_odd.dat", &f_jackson, LDBL_EPSILON, 0.95, 0.0, M_PI, 0.15);

////	photon_statistics_superposition_n("photon_statistics_n_harmonious_95_10_even.dat", &f_jackson, LDBL_EPSILON, 0.95, 0.1, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_harmonious_95_10_odd.dat", &f_jackson, LDBL_EPSILON, 0.95, 0.1, M_PI, 0.15);

////	photon_statistics_superposition_n("photon_statistics_n_harmonious_95_50_even.dat", &f_jackson, LDBL_EPSILON, 0.95, 0.5, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_harmonious_95_50_odd.dat", &f_jackson, LDBL_EPSILON, 0.95, 0.5, M_PI, 0.15);

////	photon_statistics_superposition_n("photon_statistics_n_harmonious_95_100_even.dat", &f_jackson, LDBL_EPSILON, 0.95, 1.0, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_harmonious_95_100_odd.dat", &f_jackson, LDBL_EPSILON, 0.95, 1.0, M_PI, 0.15);


////	photon_statistics_superposition_n("photon_statistics_n_jackson_75_200_0_even.dat", &f_jackson, 0.75, 2.0 - LDBL_EPSILON, 0.0, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_jackson_75_200_0_odd.dat", &f_jackson, 0.75, 2.0 - LDBL_EPSILON, 0.0, M_PI, 0.15);

////	photon_statistics_superposition_n("photon_statistics_n_jackson_75_200_10_even.dat", &f_jackson, 0.75, 2.0 - LDBL_EPSILON, 0.1, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_jackson_75_200_10_odd.dat", &f_jackson, 0.75, 2.0 - LDBL_EPSILON, 0.1, M_PI, 0.15);

////	photon_statistics_superposition_n("photon_statistics_n_jackson_75_200_50_even.dat", &f_jackson, 0.75, 2.0 - LDBL_EPSILON, 0.5, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_jackson_75_200_50_odd.dat", &f_jackson, 0.75, 2.0 - LDBL_EPSILON, 0.5, M_PI, 0.15);

////	photon_statistics_superposition_n("photon_statistics_n_jackson_75_200_100_even.dat", &f_jackson, 0.75, 2.0 - LDBL_EPSILON, 1.0, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_jackson_75_200_100_odd.dat", &f_jackson, 0.75, 2.0 - LDBL_EPSILON, 1.0, M_PI, 0.15);


////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_75_150_0_even.dat", &f_macfarlane, 0.75, 1.5, 0.0, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_75_150_0_odd.dat", &f_macfarlane, 0.75, 1.5, 0.0, M_PI, 0.15);

////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_75_150_10_even.dat", &f_macfarlane, 0.75, 1.5, 0.1, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_75_150_10_odd.dat", &f_macfarlane, 0.75, 1.5, 0.1, M_PI, 0.15);

////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_75_150_50_even.dat", &f_macfarlane, 0.75, 1.5, 0.5, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_75_150_50_odd.dat", &f_macfarlane, 0.75, 1.5, 0.5, M_PI, 0.15);

////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_75_150_100_even.dat", &f_macfarlane, 0.75, 1.5, 1.0, 0.0, -0.15);
////	photon_statistics_superposition_n("photon_statistics_n_macfarlane_75_150_100_odd.dat", &f_macfarlane, 0.75, 1.5, 1.0, M_PI, 0.15);

	long double interval = INTERVAL_T / 10.0;
/*
	for( auto t = LDBL_EPSILON; t <= INTERVAL_T; t += interval )
	{
		char filename[128];
		sprintf_s(filename, sizeof(filename), "wigner_canonical_100_%i_even.dat", (int) (t * 100));
		wigner_superposition_batouli(filename, &f_jackson, 1.0 - LDBL_EPSILON, 1.0, t, 0.0);
	}
*/
/*
	for( auto t = LDBL_EPSILON; t <= INTERVAL_T; t += interval )
	{
		char filename[128];
		sprintf_s(filename, sizeof(filename), "wigner_canonical_100_%i_odd.dat", (int) (t * 100));
		wigner_superposition_batouli(filename, &f_jackson, 1.0 - LDBL_EPSILON, 1.0, t, M_PI);
	}
*/

/*
	for( auto t = LDBL_EPSILON; t <= INTERVAL_T; t += interval )
	{
		char filename[128];
		sprintf_s(filename, sizeof(filename), "wigner_harmonious_50_%i_even.dat", (int) (t * 100));
		wigner_superposition_batouli(filename, &f_jackson, LDBL_EPSILON, 0.5, t, 0.0);
	}
*/
/*
	for( auto t = LDBL_EPSILON; t <= INTERVAL_T; t += interval )
	{
		char filename[128];
		sprintf_s(filename, sizeof(filename), "wigner_harmonious_50_%i_odd.dat", (int) (t * 100));
		wigner_superposition_batouli(filename, &f_jackson, LDBL_EPSILON, 0.5, t, M_PI);
	}
*/

/*
	for( auto t = LDBL_EPSILON; t <= INTERVAL_T; t += interval )
	{
		char filename[128];
		sprintf_s(filename, sizeof(filename), "wigner_jackson_75_200_%i_even.dat", (int) (t * 100));
		wigner_superposition_batouli(filename, &f_jackson, 0.75, 2.0 - LDBL_EPSILON, t, 0.0);
	}
*/
/*
	for( auto t = LDBL_EPSILON; t <= INTERVAL_T; t += interval )
	{
		char filename[128];
		sprintf_s(filename, sizeof(filename), "wigner_jackson_75_200_%i_odd.dat", (int) (t * 100));
		wigner_superposition_batouli(filename, &f_jackson, 0.75, 2.0 - LDBL_EPSILON, t, M_PI);
	}
*/

/*
	for( auto t = LDBL_EPSILON; t <= INTERVAL_T; t += interval )
	{
		char filename[128];
		sprintf_s(filename, sizeof(filename), "wigner_macfarlane_75_150_%i_even.dat", (int) (t * 100));
		wigner_superposition_batouli(filename, &f_macfarlane, 0.75, 1.5, t, 0.0);
	}
*/
/*	
	for( auto t = LDBL_EPSILON; t <= INTERVAL_T; t += interval )
	{
		char filename[128];
		sprintf_s(filename, sizeof(filename), "wigner_macfarlane_75_150_%i_odd.dat", (int) (t * 100));
		wigner_superposition_batouli(filename, &f_macfarlane, 0.75, 1.5, t, M_PI);
	}
*/
	system("pause");
	return 0;
}