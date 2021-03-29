#ifndef INSARINVERSE_H_
#define INSARINVERSE_H_

#include <vector>
#include <fstream>
//#include <climits>
//#include <cmath>
#include <random>
#include <algorithm>
#include <map>
#include <omp.h>
#include <float.h>
#include "okadapointsourcesurface.h"

class InSAR {

		std::string param_filename_, data_filename_;

	public:

		InSAR(const std::string& param_filename, const std::string& data_filename);

		double GetParameter(std::string parameter);

		void TranslateToASCII ();

		void Binner();

		void UnitConverter(const std::string& filename);

		static void UnitConverterNepalQuake(const std::string& filename);

		static void BinnerNepalQuake(double latstart, double lonstart,
									 double latspacing, double lonspacing,
									 int latlines, int lonsamples);

		static void BinnerNepalQuake(std::vector<std::string> files);

		static void GeneratePointsNepalQuake(const int xpoints, const int ypoints,
				                             const double xmin, const double xmax,
											 const double ymin, const double ymax);

		void GeneratePoints();
};

class PointSource {

		double x_coord_, y_coord_, z_coord_, strike_angle_,
		       dip_angle_, rake_angle_, seismic_moment_;

	public:

		PointSource (const double x, const double y, const double z,
				     const double strike, const double dip,
					 const double rake, const double M0);

		double GetX() {return x_coord_;}
		double GetY() {return y_coord_;}
		double GetZ() {return z_coord_;}

		double GetStrikeAngle() {return strike_angle_;}
		double GetDipAngle() {return dip_angle_;}
		double GetRakeAngle() {return rake_angle_;}
		double GetSeismicMoment() {return seismic_moment_;}

		void SetX(double x) {x_coord_ = x;}
		void SetY(double y) {y_coord_ = y;}
		void SetZ(double z) {z_coord_ = z;}

		void SetStrikeAngle(double strike) {strike_angle_ = strike;}
		void SetDipAngle(double dip) {dip_angle_ = dip;}
		void SetRakeAngle(double rake) {rake_angle_ = rake;}
		void SetSeismicMoment(double M0) {seismic_moment_ = M0;}

		PointSource CopySource();
};

class System {

		int index_ = -1;
		std::vector<PointSource> sources_;

	public:

		System();

		std::vector<PointSource>::const_iterator begin() {
			return sources_.begin();
		}

		std::vector<PointSource>::const_iterator end() {
			return sources_.end();
		}

		void PlaceSource(const double x, const double y, const double z,
						 const double strike, const double dip,
						 const double rake, const double M0) {
			sources_.push_back(PointSource(x, y, z, strike, dip, rake, M0));
		}

		void AddSource(PointSource point) {sources_.push_back(point);}

		std::vector<PointSource>& GetSources() {return sources_;}

		void ClearSources() {sources_.clear();}

		//copies the sources in sys to the calling System object
		void CopySources(System sys);

		void SetIndex(int index) {index_ = index;}
		int GetIndex() {return index_;}

};

struct Parameters {

		double xmin_ = 0; //Minimum x value
		double xmax_ = 0; //Maximum x value
		double ymin_ = 0; //Minimum y value
		double ymax_ = 0; //Maximum y value
		double zmin_ = 0; //Minimum z value
		double zmax_ = 0; //Maximum z value
		double min_moment_ = 0; //Minimum seismic moment value
		double max_moment_ = 0; //Maximum seismic moment value
		double chance_to_mutate_ = 0; //chance for a given member of a population to mutate

		int xdatapoints_ = 0; //Number of data points in the x direction
		int ydatapoints_ = 0; //Number of data points in the y direction
		int gen_num_sources_ = 0; //Number of sources used when generating data
		int fit_num_sources_ = 0; //Number of sources used when fitting data
		int num_generations_ = 0; //Number of generations used when fitting data
		int pop_size_ = 0; //Size of population used when fitting data

		std::string source_type_ = "notset"; //Type of sources used when generating/fitting data

		Parameters();

		Parameters(const std::string& param_filename);
};

class DataGenerator {

		Parameters params_;

		std::string filename_, filename_noext_;

		static std::mt19937 engine;

	public:

		DataGenerator(const std::string& param_filename,
					  const std::string int_filename = "interferogram.txt");

		std::vector<std::array<double, 3> > Interferogram();
};

class ParameterReader {

		std::map<std::string, std::string> params_;

	public:

		ParameterReader(const std::string& param_filename);

		std::string GetParameter(const std::string& param_name);

};

class DataReader {

		std::vector<std::array<double, 3> > data_;
		std::string datafilenoext_;

	public:

		DataReader(const std::string& data_filename);

		std::vector<std::array<double, 3> > Subset(const double xmin,
												   const double xmax,
												   const double ymin,
												   const double ymax);

		std::vector<double> FindBounds();
};

class Fitter {

		std::vector<System> models_;
		std::vector<double> errors_;
		System bestmodel_;
		double bestmodelerror_ = DBL_MAX;
		std::vector<std::array<double, 3> > data_;
		std::vector<std::tuple<int, double, double, double> > gen_error_;
		static std::mt19937 engine;

		Parameters params_;

	public:

		std::string datafilenoext_;

		Fitter();

		Fitter(const std::string& param_filename, const std::string& data_filename);

		std::vector<System>& GetModels() {return models_;}
		std::vector<double>& GetErrors() {return errors_;}

		void ImportParameters(const std::string& param_filename);

		//file format is x y displacement -- separated by tabs, each triplet on a new line
		void ImportData(const std::string& data_filename);

//		int GetDataXMax();
//
//		int GetDataXMin();
//
//		int GetDataYMax();
//
//		int GetDataYMin();

		std::vector<double> GetDataDisplacementField();

		std::vector<double> GetDataUncertainty();

		std::vector<double> CalcDisplacementField(System model);

		std::vector<System> Crossover(System parent1, System parent2);

		std::vector<System> SimulatedBinaryCrossover(System parent1, System parent2);

		void LocationMutation (System &sys);

		void OrientationMutation (System &sys);

		System SelectParent();

		void GeneticAlgorithm(bool save, std::string initial_state = "");

		double TestFitness(std::vector<double> datadisp, std::vector<double> modeldisp);

		void RecordModel(System model);

		void RecordBestModel();

		void RecordFittedData();

		System& GetBestGenModel();

		double& GetBestGenError();

		std::tuple<System, double> GetBestGenModelAndError();

		System& GetBestOverallModel() {return bestmodel_;};

		double& GetBestOverallError() {return bestmodelerror_;};

};

class ModelAnalysis {

		std::vector<System> bestmodels_;
		std::vector<double> modelerrors_;
		std::vector<std::vector<double> > modeldata_;
		Fitter fitter_;

	public:

		ModelAnalysis(const std::string param_filename,
				   	  const std::string data_filename = "");

		void Run(const int iterations = 1);

		void RecordStats();

		void SensAnalysis(System sys, const int steps = 100);

		void TwoDSensAnalysis(System sys, const int steps = 100);
};

#endif
