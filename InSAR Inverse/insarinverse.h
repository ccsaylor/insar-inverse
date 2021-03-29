#ifndef INSARINVERSE_H_
#define INSARINVERSE_H_

#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <map>
#include <omp.h>
#include <float.h>
#include "okadapointsourcesurface.h"


//extra functions used to get the data into a usable format
//conversion from latitude/longitude to x/y coordinates
//two-dimensional binning to reduce computational complexity
//these will change depending on the starting format of the dataset
class InSAR {

		std::string param_filename_, data_filename_;

	public:

//		InSAR(const std::string& param_filename, const std::string& data_filename);

//		double GetParameter(std::string parameter);

//		void TranslateToASCII ();

//		void Binner();

//		void UnitConverter(const std::string& filename);

		static void UnitConverterNepalQuake(const std::string& filename);

//		static void BinnerNepalQuake(double latstart, double lonstart,
//									 double latspacing, double lonspacing,
//									 int latlines, int lonsamples);

		static void BinnerNepalQuake(std::vector<std::string> files);

//		static void GeneratePointsNepalQuake(const int xpoints, const int ypoints,
//				                             const double xmin, const double xmax,
//											 const double ymin, const double ymax);

		void GeneratePoints();
};

//stores the parameters needed to fully describe a seismic point source
//includes both position and orientation
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

//stores a set of PointSource objects
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

//stores the parameters necessary to run the fitting algorithm
struct Parameters {

		double xmin_ = 0; //Minimum x value
		double xmax_ = 0; //Maximum x value
		double ymin_ = 0; //Minimum y value
		double ymax_ = 0; //Maximum y value
		double zmin_ = 0; //Minimum z value
		double zmax_ = 0; //Maximum z value
		double min_moment_ = 0; //Minimum seismic moment value
		double max_moment_ = 0; //Maximum seismic moment value
		double chance_to_mutate_ = 0; //chance for a given member of a System to mutate

		int xdatapoints_ = 0; //Number of data points in the x direction
		int ydatapoints_ = 0; //Number of data points in the y direction
		int gen_num_sources_ = 0; //Number of sources in each System when generating data
		int fit_num_sources_ = 0; //Number of sources in each System when fitting data
		int num_generations_ = 0; //Number of generations used when fitting data
		int pop_size_ = 0; //Number of Systems used when fitting data

		std::string source_type_ = "notset"; //Type of sources used when generating/fitting data

		Parameters();

		Parameters(const std::string& param_filename);
};

//generates a random interferogram based on the Parameters
class DataGenerator {

		Parameters params_;

		std::string filename_, filename_noext_;

		static std::mt19937 engine;

	public:

		DataGenerator(const std::string& param_filename,
					  const std::string int_filename = "interferogram.txt");

		std::vector<std::array<double, 3> > Interferogram();
};

//reads Parameters from a text file
//lines must be of the form: "parametername = value"
//one parameter per line
class ParameterReader {

		std::map<std::string, std::string> params_;

	public:

		ParameterReader(const std::string& param_filename);

		std::string GetParameter(const std::string& param_name);

};

//reads data to be fit from a text file
//lines must be of the form: "x-coordinate	y-coordinate  z-coordinate"
//one data point per line
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

//includes all functions necessary to fit the data
//stores a list of Systems that are the possible models
//as well as their error relative to the dataset
//uses a genetic algorithm to slightly alter the model parameters
//over a specified number of generations to reach a final
//model that most accurately represents the data
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

//		void ImportParameters(const std::string& param_filename);

		//file format is x y displacement -- separated by tabs, each triplet on a new line
//		void ImportData(const std::string& data_filename);

//		int GetDataXMax();
//
//		int GetDataXMin();
//
//		int GetDataYMax();
//
//		int GetDataYMin();

		//gets displacement values from the data
		std::vector<double> GetDataDisplacementField();

//		std::vector<double> GetDataUncertainty();

		//calculates the displacement values of a System based
		//on the parameters of its PointSource objects
		std::vector<double> CalcDisplacementField(System model);

		//an operator that crosses the parameters of two Systems
		//by swapping random PointSource objects between them
		std::vector<System> Crossover(System parent1, System parent2);

		//an operator that crosses the parameters of two Systems
		//using simulated binary crossover
		std::vector<System> SimulatedBinaryCrossover(System parent1, System parent2);

		//mutates the locations of PointSource objects in sys
		//based on input parameters
		void LocationMutation (System &sys);

		//mutates the orientations of PointSource objects in sys
		//based on input parameters
		void OrientationMutation (System &sys);

		//selects a random System from the current population of models
		//probabilities are weighted based on System error value
		System SelectParent();

		//runs the genetic algorithm based on the input parameters
		void GeneticAlgorithm(bool save);

		//calculates the error between the data and model
		double TestFitness(std::vector<double> datadisp, std::vector<double> modeldisp);

		//save the parameters of the PointSource objects in model
		void RecordModel(System model);

		//save the parameters of the PointSource objects
		//in the model with the lowest error
		void RecordBestModel();

		//records the data used in the fitting
		//all of it if no subset is specified
		void RecordFittedData();

		//get the System with lowest error in the current generation
		System& GetBestGenModel();

		//get the error of the System with lowest error in the
		//current generation
		double& GetBestGenError();

		std::tuple<System, double> GetBestGenModelAndError();

		System& GetBestOverallModel() {return bestmodel_;};

		double& GetBestOverallError() {return bestmodelerror_;};

};

//used to explore the sensitivity of the error of a model
//to changes in the values of the parameters
//can change either a single parameter at a time or every
//possible pairing of parameters to see if correlations exist
class ModelAnalysis {

		std::vector<System> bestmodels_;
		std::vector<double> modelerrors_;
		std::vector<std::vector<double> > modeldata_;
		Fitter fitter_;

	public:

		ModelAnalysis(const std::string param_filename,
				   	  const std::string data_filename = "");

		//runs many fits to see the spread in the final values of
		//the model parameters
		void Run(const int iterations = 1);

		void RecordStats();

		void SensAnalysis(System sys, const int steps = 100);

		void TwoDSensAnalysis(System sys, const int steps = 100);
};

#endif
