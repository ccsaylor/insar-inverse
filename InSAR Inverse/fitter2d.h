#ifndef FITTER2D_H_
#define FITTER2D_H_

#include "insarinverse.h"

class Fitter2D {

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

		Fitter2D() {};

		Fitter2D(const std::string& param_filename, const std::string& data_filename);

		std::vector<Displacement> GetDataDisplacementField();

		std::vector<Displacement> CalcDisplacementField(System model);
		std::vector<Displacement> CalcDisplacementField(System model, GreensCatalog catalog);

		double TestFitness(std::vector<Displacement> modeldisp);

		System SelectParent();

		std::vector<System> Crossover(System parent1, System parent2);

		std::vector<System> SBXMoment(System parent1, System parent2);

		void LocationMutation(System& sys);

		void XYMutation(System& sys);

		void OrientationMutation(System& sys);

		void MomentMutation (System& sys);

		std::tuple<System, double> GetBestGenModelAndError();

		void GeneticAlgorithm(bool save);

		void LocGeneticAlgorithm(bool save);

		void MomentGA(bool save, std::string sourcefile = "");

		void RecordModel(System& model);

		void RecordFittedData();

		void RecordBestModel() {RecordModel(bestmodel_);}

		void SaveGreensCatalog(const std::string& filename);

		void Setup();
};



#endif /* FITTER2D_H_ */
