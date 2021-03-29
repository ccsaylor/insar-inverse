#include "insarinverse.h"

#include <cmath>
#include <random>

int main() {

//	InSAR sar("SanAnd_08519_14004-007_14145-008_0260d_s01_L090HH_01.ann","SanAnd_08519_14004-007_14145-008_0260d_s01_L090HH_01.unw.grd");
//	sar.TranslateToASCII();
//	sar.Binner();
//	sar.UnitConverter("binned_data.txt");

//	Fitter fit;
//	fit.ImportData("converted_data.txt");

//	fit.LocationGeneticAlgorithm(1000);

//	sar.GeneratePoints();
//	fit.OrientationGeneticAlgorithm(1000);

//	std::string filename = "interferogram_1.txt";
//	DataGenerator::InterferogramGenerator(filename, 0.0, 30.0, 0.0, 10.0, -3.0, -10.0, 30, 10, 10);
//	Fitter fit;
//	fit.ImportData(filename);
//	fit.BothGeneticAlgorithm(10000);

//	InSAR::BinnerNepalQuake({"los_T048_POST01_detrend.lltnde","los_T048_C02_detrend.lltnde"});
//	InSAR::UnitConverterNepalQuake("total_displacement_binned_50x50.txt");

//	Fitter fit;
//	fit.ImportParameters("parametertest.txt");
//	fit.ImportData("total_displacement_binned_50x50_converted.txt");
//	fit.GeneticAlgorithm("total_displacement_binned_50x50_converted_model_points.txt");

//	double xmin = 0.0;
//	double xmax = 10.0;
//	double ymin = 0.0;
//	double ymax = 10.0;
//	double zmin = -3.0;
//	double zmax = -10.0;
//	int xpoints = 20;
//	int ypoints = 20;
//	int numpoints = 1;

//	Fitter fit;
//	fit.ImportParameters("parametertest.txt");

//	DataGenerator gen("modelanalysisparams.txt", "interferogram.txt");
//	gen.Interferogram();

	ModelAnalysis modan("modelanalysisparams.txt", "interferogram.txt");
	modan.Run(50);
	modan.RecordStats();

//	System sys;
//	sys.PlaceSource(48.5774, 43.6555, -10.4041, 5.25892, 0.523487, 0, 8.29375e+09);
////	modan.SensAnalysis(sys);
//	modan.TwoDSensAnalysis(sys, 50);

	return 0;
}
