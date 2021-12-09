#ifndef INSARINVERSE_H_
#define INSARINVERSE_H_

#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <map>
//#include <omp.h>
#include <float.h>
//#include <unordered_set>
#include "okadapointsourcesurface.h"
#include "okadasurface.h"

struct Position {
		double x = 0;
		double y = 0;
		double z = 0;

		Position() {};

		Position(const double x, const double y, const double z);
};

struct Data {
		Position pos;
		double dx = 0;
		double dy = 0;
		double dz = 0;

		Data() {};

		Data(Position& pos, double dx, double dy, double dz);
};

struct Displacement {
		double dx = 0;
		double dy = 0;
		double dz = 0;

		Displacement() {};

		Displacement(double dx, double dy, double dz);
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

//stores the parameters necessary to run the fitting algorithm
struct Parameters {

//		std::unordered_set<std::string> found;

		double xmin = 0; //Minimum data x value
		double xmax = 0; //Maximum data x value
		double ymin = 0; //Minimum data y value
		double ymax = 0; //Maximum data y value
		double min_moment = 0; //Minimum seismic moment value
		double max_moment = 0; //Maximum seismic moment value
		double chance_to_mutate = 0; //chance for a given member of a System to mutate
		double strike_angle = 0;
		double dip_angle = 0;
		double rake_angle = 0;

		int xdatapoints = 0; //Number of data points in the x direction
		int ydatapoints = 0; //Number of data points in the y direction
		int gen_num_sources = 0; //Number of sources in each System when generating data
		int fit_num_sources = 0; //Number of sources in each System when fitting data
		int num_generations = 0; //Number of generations used when fitting data
		int pop_size = 0; //Number of Systems used when fitting data
		int xfitpoints = 0; //Number of sources in the x direction, used in MomentGA function
		int yfitpoints = 0; //Number of sources in the y direction
		int zfitpoints = 0; //Number of sources in the z direction

		double sxmin = 0; //Minimum source x value
		double sxmax = 0; //Maximum source x value
		double symin = 0; //Minimum source y value
		double symax = 0; //Maximum source y value
		double szmin = 0; //Minimum source z value
		double szmax = 0; //Maximum source z value

//		std::string source_type = "notset"; //Type of sources used when generating/fitting data

		Parameters() {};

		Parameters(const std::string& param_filename);

		bool DataLimits();
		bool SourceLimits();
};

//extra functions used to get the data into a usable format
//conversion from latitude/longitude to x/y coordinates
//two-dimensional binning to reduce computational complexity
//these will change depending on the starting format of the dataset
class InSAR {

	public:

		static void UnitConverterNepalQuake(const std::string& filename);

		static void BinnerNepalQuake(std::vector<std::string> files);

		static void MultiBinnerNepalQuake(std::vector<std::string> files);

		static void TwoDBinner(const std::string& filename);

};

//stores the parameters needed to fully describe a seismic point source
//includes both position and orientation
//class PointSource { //TODO: Change to struct?
//
//		Position pos_;
//		double strike_angle_ = 0;
//		double dip_angle_ = 0;
//		double rake_angle_ = 0;
//		double seismic_moment_ = 0;
//
//	public:
//
//		PointSource(const double x, const double y, const double z,
//					const double strike = 0, const double dip = 0,
//					const double rake = 0, const double moment = 0);
//
//		double GetX() {return pos_.x;}
//		double GetY() {return pos_.y;}
//		double GetZ() {return pos_.z;}
//
//		double GetStrikeAngle() {return strike_angle_;}
//		double GetDipAngle() {return dip_angle_;}
//		double GetRakeAngle() {return rake_angle_;}
//		double GetSeismicMoment() {return seismic_moment_;}
//
//		void SetX(double x) {pos_.x = x;}
//		void SetY(double y) {pos_.y = y;}
//		void SetZ(double z) {pos_.z = z;}
//
//		void SetStrikeAngle(double strike) {strike_angle_ = strike;}
//		void SetDipAngle(double dip) {dip_angle_ = dip;}
//		void SetRakeAngle(double rake) {rake_angle_ = rake;}
//		void SetSeismicMoment(double M0) {seismic_moment_ = M0;}
//
//		PointSource CopySource();
//};

struct PointSource {

		Position pos;
		double strike = 0;
		double dip = 0;
		double rake = 0;
		double moment = 0;

		PointSource(const double x, const double y, const double z,
					const double strike = 0, const double dip = 0,
					const double rake = 0, const double moment = 0);

		PointSource CopySource();
};

//stores a vector of PointSource objects
class System {

		int index_ = -1;
		std::vector<PointSource> sources_;

	public:

		System() {};

		System(const std::string& filename);

		System(Parameters params);

		std::vector<PointSource>::iterator begin() {
			return sources_.begin();
		}

		std::vector<PointSource>::iterator end() {
			return sources_.end();
		}

		int size() {return sources_.size();}

		void PlaceSource(const double x, const double y, const double z,
						 const double strike, const double dip,
						 const double rake, const double M0) {
			sources_.emplace_back(x, y, z, strike, dip, rake, M0);
		}

		void AddSource(PointSource point) {sources_.push_back(point);}

		std::vector<PointSource>& GetSources() {return sources_;}

		void ClearSources() {sources_.clear();}

		//copies the sources in sys to the calling System object
		void CopySources(System sys);

		void SetIndex(int index) {index_ = index;}
		int GetIndex() {return index_;}

		void Displacement(const std::string& filename, int numx, int numy,
						  double xmin, double xmax, double ymin, double ymax);

};

//generates a random interferogram based on the Parameters
class DataGenerator {

		Parameters params_;

		std::string filename_, filename_noext_;

		static std::mt19937 engine;

	public:

		DataGenerator(const std::string& param_filename,
					  const std::string int_filename = "interferogram.txt");

		std::vector<Data> Interferogram();
};

//reads data to be fit from a text file
//lines must be of the form: "x-coordinate	y-coordinate  z-coordinate"
//one data point per line
class DataReader {

		std::vector<Data> data_;
		std::string datafilenoext_;

	public:

		DataReader(const std::string& data_filename);

		std::vector<Data> GetData();

		std::vector<Data> Subset(const double xmin,
								 const double xmax,
								 const double ymin,
								 const double ymax);

		std::vector<double> FindBounds();

		int CountColumns(const std::string& filename);
};

struct Greens {

		double u_x = 0;
		double u_y = 0;
		double u_z = 0;

		Greens() {};
};

struct GreensCatalog {

		int numsources = 0;
		int numdatapoints = 0;

//		std::vector<StrikeGreens> strikes;
//		std::vector<DipGreens> dips;
		std::vector<Greens> greens;

		GreensCatalog() {};

		GreensCatalog(System sys, std::vector<Data> data);

//		void SaveCatalog(std::ofstream& fout);
		void SaveCatalog(const std::string& filename);
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
		std::vector<Data> data_;
		std::vector<std::tuple<int, double, double, double> > gen_error_;
		static std::mt19937 engine;

		Parameters params_;

	public:

		std::string datafilenoext_;

		Fitter() {};

		Fitter(const std::string& param_filename, const std::string& data_filename);

		std::vector<System>& GetModels() {return models_;}
		std::vector<double>& GetErrors() {return errors_;}

		//gets displacement values from the data
		std::vector<double> GetDataDisplacementField();

		//calculates the displacement values of a System based
		//on the parameters of its PointSource objects
		std::vector<double> CalcDisplacementField(System model);

		std::vector<double> CalcDisplacementField(System model, GreensCatalog catalog);

		//an operator that crosses the parameters of two Systems
		//by swapping random PointSource objects between them
		std::vector<System> Crossover(System parent1, System parent2);

		//an operator that crosses all parameters of two Systems
		//using simulated binary crossover
		std::vector<System> SimulatedBinaryCrossover(System parent1, System parent2);

		//an operator that crosses the seismic moments of two Systems
		//using simulated binary crossover
		std::vector<System> SBXMoment(System parent1, System parent2);

		//mutates the locations of PointSource objects in sys
		//based on input parameters
		void LocationMutation (System &sys);

		//mutates the orientations & moments of PointSource objects in sys
		//based on input parameters
		void OrientationMutation (System &sys);

		//mutates the seismic moment of PointSource objects in sys
		//based on input parameters
		void MomentMutation (System &sys);

		//selects a random System from the current population of models
		//probabilities are weighted based on System error value
		System SelectParent();

		//runs the genetic algorithm based on the input parameters
		void GeneticAlgorithm(bool save);

		//runs the moment genetic algorithm based on the input parameters
		void MomentGA(bool save, std::string sourcefile = "");

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

//		void SaveGreensCatalog(std::ofstream& fout);
		void SaveGreensCatalog(const std::string& filename);

		void Setup();
};

#endif
