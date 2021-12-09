#include "insarinverse.h"
#include "okadasurface.h"
#include "okadapointsourcesurface.h"
#include "fitter2d.h"

#include <cmath>
#include <random>
#include <iostream>
#include <chrono>

int main() {

	std::chrono::high_resolution_clock::time_point start =
			std::chrono::high_resolution_clock::now();

//	InSAR::BinnerNepalQuake({"los_T048_C02_detrend.lltnde"});
//	InSAR::UnitConverterNepalQuake("nepal_binned_75x75.txt");

//	InSAR::MultiBinnerNepalQuake({"los_T048_C02_detrend.lltnde","los_T157_C01_detrend.lltnde"});

	Fitter fit("nepalparams.txt", "nepal_binned_75x75_converted.txt");
//	fit.GeneticAlgorithm(true);
	fit.Setup();

//	Fitter fit("iranparams.txt", "iran_displacements_inc.txt");
//	fit.SaveGreensCatalog("greens.txt");
//	fit.MomentGA(true);

//	Fitter2D fit("iranparams.txt", "iran_binned.txt");
//	Fitter2D fit("ridgecrestparams.txt", "ridgecrest_binned.txt");
//	fit.Setup();
//	fit.SaveGreensCatalog("greens.txt");
//	fit.RecordFittedData();
//	fit.MomentGA(true);
//	fit.LocGeneticAlgorithm(true);

	std::chrono::high_resolution_clock::time_point end =
			std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> diff = end - start;

	std::cout << "Elapsed time: " << diff.count() << "s." << "\n";

	return 0;
}
