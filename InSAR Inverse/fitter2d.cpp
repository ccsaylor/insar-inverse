#include "fitter2d.h"
#include "okadasurface.h"
#include <chrono>
#include <iostream>
#include <cmath>

std::mt19937 Fitter2D::engine(time(0)); //TODO: Find a faster random number generator

Fitter2D::Fitter2D(const std::string& param_filename, const std::string& data_filename) {

	size_t lastindex = data_filename.find_last_of(".");
	datafilenoext_ = data_filename.substr(0, lastindex);

	params_ = Parameters(param_filename);

	DataReader data(data_filename);

	if (params_.DataLimits()) {
		data_ = data.Subset(params_.xmin, params_.xmax,
							params_.ymin, params_.ymax);
	}
	else {
		data_ = data.GetData();

		std::vector<double> bounds = data.FindBounds();
		params_.xmin = bounds[0];
		params_.xmax = bounds[1];
		params_.ymin = bounds[2];
		params_.ymax = bounds[3];
	}

	if (!params_.SourceLimits()) {
		params_.sxmin = params_.xmin;
		params_.sxmax = params_.xmax;
		params_.symin = params_.ymin;
		params_.symax = params_.ymax;
	}
}

std::vector<Displacement> Fitter2D::GetDataDisplacementField() {
	std::vector<Displacement> displacements;
	for (Data& data : data_) {
		displacements.emplace_back(data.dx, data.dy, data.dz);
	}
	return displacements;
}

std::vector<Displacement> Fitter2D::CalcDisplacementField(System model) { //TODO: Update to use new Okada equations
	std::vector<Displacement> displacements(data_.size(), Displacement(0.0, 0.0, 0.0));
	std::vector<PointSource>& sources = model.GetSources();

	int i = 0;
	for (Data& data : data_) {
		double datax = data.pos.x;
		double datay = data.pos.y;

//		displacements.emplace_back(0.0, 0.0, 0.0);

		for (PointSource& source : sources) {

			double sx = source.pos.x;
			double sy = source.pos.y;
			double sz = source.pos.z;
			double strike = source.strike;
			double dip = source.dip;
			double rake = source.rake;
			double moment = source.moment;

			double newx = (datax - sx)*cos_(strike) - (datay - sy)*sin_(strike);
			double newy = (datay - sy)*sin_(strike) + (datay - sy)*cos_(strike);

			double strx = strike_u_x_0(newx, newy, std::abs(sz), dip, moment);
			double dipx = dip_u_x_0(newx, newy, std::abs(sz), dip, moment);
			double strz = strike_u_z_0(newx, newy, std::abs(sz), dip, moment);
			double dipz = dip_u_z_0(newx, newy, std::abs(sz), dip, moment);

			displacements[i].dx += strx*cos_(rake) + dipx*sin_(rake);
			displacements[i].dz += strz*cos_(rake) + dipz*sin_(rake);
		}
		i++;
	}
	return displacements;
}

std::vector<Displacement> Fitter2D::CalcDisplacementField(System model, GreensCatalog catalog) {
	std::vector<Displacement> displacements(data_.size(), Displacement(0.0, 0.0, 0.0));
	std::vector<PointSource>& sources = model.GetSources();

	int i = 0;

	for (Displacement& disp : displacements) {
		for (PointSource& source : sources) {

			double moment = source.moment;

			disp.dx += catalog.greens[i].u_x*moment;
			disp.dz += catalog.greens[i].u_z*moment;

			i++;
		}
	}
	return displacements;
}

System Fitter2D::SelectParent() { //TODO: Update to run only a single time to pick all needed pairings
	System parent;
	std::vector<double> weights;
	std::vector<double> errors = errors_;
	double sum;

	for (std::vector<double>::iterator it = errors.begin(); it != errors.end(); ++it) {
		*it = 1 / *it;
	}

	sum = std::accumulate(errors.begin(), errors.end(), 0.0);

	for (std::vector<double>::iterator it = errors.begin(); it != errors.end(); ++it) {
		weights.push_back(*it/sum);
	}

	std::discrete_distribution<int> dist(weights.begin(), weights.end());
	parent = models_[dist(engine)];

	return parent;
}

std::vector<System> Fitter2D::Crossover(System parent1, System parent2) {

	std::vector<PointSource>& parent1_sources = parent1.GetSources();
	std::vector<PointSource>& parent2_sources = parent2.GetSources();

	std::vector<System> new_sys;
	std::uniform_real_distribution<double> zerotoone(0.0, 1.0);

	for (int i = 0; i < static_cast<int>(parent1_sources.size()); i++) {
		if (zerotoone(engine) < 0.5) {
			PointSource temp = parent1_sources[i];
			parent1_sources[i] = parent2_sources[i];
			parent2_sources[i] = temp;
		}
	}

	new_sys.push_back(parent1);
	new_sys.push_back(parent2);

	return new_sys;
}

std::vector<System> Fitter2D::SBXMoment(System parent1, System parent2) {

	double beta, eta_c, u;
	double p1mom, p2mom;
	double p1mom_new, p2mom_new;

	std::vector<PointSource>& parent1_sources = parent1.GetSources();
	std::vector<PointSource>& parent2_sources = parent2.GetSources();
	if (parent1_sources.size() != parent2_sources.size()) {
		throw std::runtime_error("Parents don't have the same number of sources");
	}
	std::vector<PointSource>::iterator pit, sit;

	std::vector<System> new_sys;

	eta_c = 2;
	std::uniform_real_distribution<double> zerotoone(0.0, 1.0);
	u = zerotoone(engine);

	if (u <= 0.5) {
		beta = 2*std::pow(u, 1/(eta_c + 1));
	}
	else {
		beta = std::pow((0.5/(1 - u)), 1/(eta_c + 1));
	}

	for (pit = parent1_sources.begin(), sit = parent2_sources.begin(); pit != parent1_sources.end() && sit != parent2_sources.end(); ++pit, ++sit) {

		p1mom = pit->moment;
		p2mom = sit->moment;

		p1mom_new = 0.5*((1 + beta)*p1mom + (1 - beta)*p2mom);
		p2mom_new = 0.5*((1 - beta)*p1mom + (1 + beta)*p2mom);

		pit->moment = p1mom_new;
		sit->moment = p2mom_new;
	}

	new_sys.push_back(parent1);
	new_sys.push_back(parent2);

	return new_sys;

}

void Fitter2D::LocationMutation (System &sys) {
	std::vector<PointSource>& sources = sys.GetSources();
	std::uniform_real_distribution<double> zerotoone(0.0, 1.0);
	double x, y, z;

	for (std::vector<PointSource>::iterator it = sources.begin(); it != sources.end(); ++it) {
		x = it->pos.x;
		y = it->pos.y;
		z = it->pos.z;
		std::normal_distribution<double> x_spread(x, 2);
		std::normal_distribution<double> y_spread(y, 2);
		std::normal_distribution<double> z_spread(z, 0.5);

		if (zerotoone(engine) < 0.1) {
			double new_x = x_spread(engine);
			while (new_x > params_.sxmax || new_x < params_.sxmin) {
				new_x = x_spread(engine);
			}
			double new_y = y_spread(engine);
			while (new_y > params_.symax || new_y < params_.symin) {
				new_y = y_spread(engine);
			}
			double new_z = z_spread(engine);
			while (new_z > 0) {
				new_z = z_spread(engine);
			}
			it->pos.x = new_x;
			it->pos.y = new_y;
			it->pos.z = new_z;
		}
	}
}

void Fitter2D::XYMutation (System &sys) {
	std::vector<PointSource>& sources = sys.GetSources();
	std::uniform_real_distribution<double> zerotoone(0.0, 1.0);
	double x, y, moment;

	for (std::vector<PointSource>::iterator it = sources.begin(); it != sources.end(); ++it) {
		x = it->pos.x;
		y = it->pos.y;
		moment = it->moment;
		std::normal_distribution<double> x_spread(x, 2);
		std::normal_distribution<double> y_spread(y, 2);
		std::normal_distribution<double> moment_spread(std::log10(moment), 0.5);

		if (zerotoone(engine) < 0.1) {
			double new_x = x_spread(engine);
			while (new_x > params_.sxmax || new_x < params_.sxmin) {
				new_x = x_spread(engine);
			}
			double new_y = y_spread(engine);
			while (new_y > params_.symax || new_y < params_.symin) {
				new_y = y_spread(engine);
			}
			double new_moment = pow(10.0, moment_spread(engine));
			it->pos.x = new_x;
			it->pos.y = new_y;
			it->moment = new_moment;
		}
	}
}

void Fitter2D::OrientationMutation (System &sys) {
	std::vector<PointSource>& sources = sys.GetSources();
	std::uniform_real_distribution<double> zerotoone(0.0, 1.0);
	double strike, dip, moment;

	for (std::vector<PointSource>::iterator it = sources.begin(); it != sources.end(); ++it) {
		strike = it->strike;
		dip = it->dip;
		moment = it->moment;
		std::normal_distribution<double> strike_spread(strike, M_PI/6.0);
		std::normal_distribution<double> dip_spread(dip, M_PI/24.0);
		std::normal_distribution<double> moment_spread(std::log10(moment), 0.5);

		if (zerotoone(engine) < 0.1) {
			double new_strike = strike_spread(engine);
			double new_dip = dip_spread(engine);
			while (new_dip > M_PI/180*90 || new_dip < 0) {
				new_dip = dip_spread(engine);
			}
			double new_moment = pow(10.0, moment_spread(engine));
			it->strike = new_strike;
			it->dip = new_dip;
			it->moment = new_moment;
		}
	}
}

void Fitter2D::MomentMutation (System &sys) {
	std::vector<PointSource>& sources = sys.GetSources();
	std::uniform_real_distribution<double> zerotoone(0.0, 1.0);

	for (std::vector<PointSource>::iterator it = sources.begin(); it != sources.end(); ++it) {
		double moment = it->moment;
		std::normal_distribution<double> moment_spread(std::log10(moment), 0.25);

		if (zerotoone(engine) < 0.1) {
			double new_moment = pow(10.0, moment_spread(engine));
			it->moment = new_moment;
		}
	}
}

std::tuple<System, double> Fitter2D::GetBestGenModelAndError() {
	std::vector<double>::iterator min = std::min_element(errors_.begin(), errors_.end());
	int index = std::distance(errors_.begin(), min);
	return std::make_tuple(models_[index], errors_[index]);
}

void Fitter2D::GeneticAlgorithm(bool save) { //TODO: Add rake as a parameter

	models_.clear();
	errors_.clear();
	bestmodelerror_ = DBL_MAX;
	gen_error_.clear();

	std::uniform_int_distribution<int> dist_x(params_.sxmin, params_.sxmax);
	std::uniform_int_distribution<int> dist_y(params_.symin, params_.symax);
	std::uniform_int_distribution<int> dist_z(params_.szmin, params_.szmax);
	std::uniform_real_distribution<double> strike(0.0, 2.0*M_PI);
	std::uniform_real_distribution<double> diprake(0.0, M_PI/2.0);
	std::uniform_real_distribution<double> moment(params_.min_moment,
												  params_.max_moment);
	std::uniform_real_distribution<double> mutate(0.0, 1.0);

	System model;
	System bestmodel;

	std::vector<System> next_gen;

	double error;

	std::vector<Displacement> datadisp = this->GetDataDisplacementField();

	for (int i = 0; i < params_.pop_size; i++) { //generate random first generation parents

		model.ClearSources();

		for (int j = 0; j < params_.fit_num_sources; j++) { //number of sources in each parent
			model.PlaceSource(dist_x(engine), dist_y(engine), -dist_z(engine),
							  strike(engine), diprake(engine), diprake(engine), moment(engine));
		}

		std::vector<Displacement> modeldisp = CalcDisplacementField(model);

		error = TestFitness(modeldisp);

		models_.push_back(model);
		errors_.push_back(error);

	}

	int counter = 0;
	for (int generation = 0; generation < params_.num_generations; generation++) {

		std::cout << "Generation " << generation << std::endl;

		for (std::vector<System>::iterator it = models_.begin(); it != models_.end(); ++it) {
			it->SetIndex(std::distance(models_.begin(), it));
		}

		next_gen.clear();

		while (static_cast<int>(next_gen.size()) < params_.pop_size) {
			System parent1 = SelectParent();
			System parent2 = SelectParent();
			if (parent1.GetIndex() != parent2.GetIndex()) {
				std::vector<System> crossed = Crossover(parent1, parent2);
				for (std::vector<System>::iterator it = crossed.begin(); it != crossed.end(); ++it) {
					next_gen.push_back(*it);
				}
			}
			else {
				continue;
			}
		}

		models_.clear();
		errors_.clear();

		#pragma omp parallel
		{
			std::vector<System> private_models;
			std::vector<double> private_errors;
			#pragma omp for nowait
			for (int k = 0; k < static_cast<int>(next_gen.size()); ++k) {
				std::vector<Displacement> modeldisp = CalcDisplacementField(next_gen[k]);
				error = TestFitness(modeldisp);
				private_models.push_back(next_gen[k]);
				private_errors.push_back(error);
			}
			#pragma omp critical
			{
				models_.insert(models_.end(), private_models.begin(), private_models.end());
				errors_.insert(errors_.end(), private_errors.begin(), private_errors.end());
			}
		}

//		std::cout << "Minimum error: " << *std::min_element(errors_.begin(), errors_.end()) << "\n";
		std::cout << "Average error: " << std::accumulate(errors_.begin(), errors_.end(), 0.0)/errors_.size() << "\n";
//		std::cout << "Maximum error: " << *std::max_element(errors_.begin(), errors_.end()) << "\n";

		gen_error_.push_back({generation,
							  *std::min_element(errors_.begin(), errors_.end()),
							  std::accumulate(errors_.begin(), errors_.end(), 0.0)/errors_.size(),
							  *std::max_element(errors_.begin(), errors_.end())});

		System best_gen_model;
		double best_gen_error;

		std::tie(best_gen_model, best_gen_error) = GetBestGenModelAndError();

		if (best_gen_error < bestmodelerror_) {
			bestmodel_ = best_gen_model;
			bestmodelerror_ = best_gen_error;
			counter = 0;
		}

		else {
			counter++;
		}

		if (counter > 500) {
			break;
		}

		for (std::vector<System>::iterator it = models_.begin(); it != models_.end(); ++it) {

			if (mutate(engine) < params_.chance_to_mutate) {
				OrientationMutation(*it);
				LocationMutation(*it);
			}
		}
	}

	if (save) {
		std::ofstream ofile(datafilenoext_ + "_error.txt");
		for (std::tuple<int, double, double, double> tup : gen_error_) {
			ofile << (std::get<0>(tup)) << "\t";
			ofile << (std::get<1>(tup)) << "\t";
			ofile << (std::get<2>(tup)) << "\t";
			ofile << (std::get<3>(tup)) << "\n";
		}
		RecordBestModel();
	}
}

void Fitter2D::LocGeneticAlgorithm(bool save) {

	models_.clear();
	errors_.clear();
	bestmodelerror_ = DBL_MAX;
	gen_error_.clear();

	std::uniform_int_distribution<int> dist_x(params_.sxmin, params_.sxmax);
	std::uniform_int_distribution<int> dist_y(params_.symin, params_.symax);
	std::uniform_int_distribution<int> dist_z(params_.szmin, params_.szmax);
	double strike = params_.strike_angle;
	double dip = params_.dip_angle;
	double rake = params_.rake_angle;
	std::uniform_real_distribution<double> moment(params_.min_moment,
												  params_.max_moment);
	std::uniform_real_distribution<double> mutate(0.0, 1.0);

	System model;
	System bestmodel;

	std::vector<System> next_gen;

	double error;

	std::vector<Displacement> datadisp = this->GetDataDisplacementField();

	for (int i = 0; i < params_.pop_size; i++) { //generate random first generation parents

		model.ClearSources();

		for (int j = 0; j < params_.fit_num_sources; j++) { //number of sources in each parent
			model.PlaceSource(dist_x(engine), dist_y(engine), -dist_z(engine),
							  strike, dip, rake, moment(engine));
		}

		std::vector<Displacement> modeldisp = CalcDisplacementField(model);

		error = TestFitness(modeldisp);

		models_.push_back(model);
		errors_.push_back(error);

	}

	int counter = 0;
	for (int generation = 0; generation < params_.num_generations; generation++) {

		std::cout << "Generation " << generation << std::endl;

		for (std::vector<System>::iterator it = models_.begin(); it != models_.end(); ++it) {
			it->SetIndex(std::distance(models_.begin(), it));
		}

		next_gen.clear();

		while (static_cast<int>(next_gen.size()) < params_.pop_size) {
			System parent1 = SelectParent();
			System parent2 = SelectParent();
			if (parent1.GetIndex() != parent2.GetIndex()) {
				std::vector<System> crossed = Crossover(parent1, parent2);
				for (std::vector<System>::iterator it = crossed.begin(); it != crossed.end(); ++it) {
					next_gen.push_back(*it);
				}
			}
			else {
				continue;
			}
		}

		models_.clear();
		errors_.clear();

		#pragma omp parallel
		{
			std::vector<System> private_models;
			std::vector<double> private_errors;
			#pragma omp for nowait
			for (int k = 0; k < static_cast<int>(next_gen.size()); ++k) {
				std::vector<Displacement> modeldisp = CalcDisplacementField(next_gen[k]);
				error = TestFitness(modeldisp);
				private_models.push_back(next_gen[k]);
				private_errors.push_back(error);
			}
			#pragma omp critical
			{
				models_.insert(models_.end(), private_models.begin(), private_models.end());
				errors_.insert(errors_.end(), private_errors.begin(), private_errors.end());
			}
		}

//		std::cout << "Minimum error: " << *std::min_element(errors_.begin(), errors_.end()) << "\n";
		std::cout << "Average error: " << std::accumulate(errors_.begin(), errors_.end(), 0.0)/errors_.size() << "\n";
//		std::cout << "Maximum error: " << *std::max_element(errors_.begin(), errors_.end()) << "\n";

		gen_error_.push_back({generation,
							  *std::min_element(errors_.begin(), errors_.end()),
							  std::accumulate(errors_.begin(), errors_.end(), 0.0)/errors_.size(),
							  *std::max_element(errors_.begin(), errors_.end())});

		System best_gen_model;
		double best_gen_error;

		std::tie(best_gen_model, best_gen_error) = GetBestGenModelAndError();

		if (best_gen_error < bestmodelerror_) {
			bestmodel_ = best_gen_model;
			bestmodelerror_ = best_gen_error;
			counter = 0;
		}

		else {
			counter++;
		}

		if (counter > 500) {
			break;
		}

		for (std::vector<System>::iterator it = models_.begin(); it != models_.end(); ++it) {

			if (mutate(engine) < params_.chance_to_mutate) {
				XYMutation(*it);
			}
		}
	}

	if (save) {
		std::ofstream ofile(datafilenoext_ + "_error.txt");
		for (std::tuple<int, double, double, double> tup : gen_error_) {
			ofile << (std::get<0>(tup)) << "\t";
			ofile << (std::get<1>(tup)) << "\t";
			ofile << (std::get<2>(tup)) << "\t";
			ofile << (std::get<3>(tup)) << "\n";
		}
		RecordBestModel();
	}
}

void Fitter2D::MomentGA(bool save, std::string sourcefile) {

	models_.clear();
	errors_.clear();
	bestmodelerror_ = DBL_MAX;
	gen_error_.clear();

	std::chrono::high_resolution_clock::time_point start =
			std::chrono::high_resolution_clock::now();

	std::uniform_real_distribution<double> moment(params_.min_moment,
												  params_.max_moment);
	std::uniform_real_distribution<double> mutate(0.0, 1.0);

	System model;
	if (!sourcefile.empty()) {
		std::ifstream fin(sourcefile);
		double x, y, z, strike, dip, rake, moment;
		while (fin >> x >> y >> z >> strike >> dip >> rake >> moment) {
			model.AddSource(PointSource(x, y, z, strike, dip));
		}
	}
	else {
		model = System(params_);
	}

	GreensCatalog catalog(model, data_);
	catalog.SaveCatalog("greens.txt");

	std::vector<System> next_gen;

	double error;

	for (int i = 0; i < params_.pop_size; i++) { //generate random first generation parents

		for (PointSource& source : model) {
			source.moment = moment(engine);
		}

		std::vector<Displacement> modeldisp = CalcDisplacementField(model, catalog);

		error = TestFitness(modeldisp);

		models_.push_back(model);
		errors_.push_back(error);

	}

	int counter = 0;
	for (int generation = 0; generation < params_.num_generations; generation++) {

		std::cout << "Generation " << generation << std::endl;

		for (std::vector<System>::iterator it = models_.begin(); it != models_.end(); ++it) { //TODO: Move outside of loop
			it->SetIndex(std::distance(models_.begin(), it));
		}

		next_gen.clear();

		while (static_cast<int>(next_gen.size()) < params_.pop_size) {
			System parent1 = SelectParent();
			System parent2 = SelectParent();
			if (parent1.GetIndex() != parent2.GetIndex()) {
//				std::vector<System> crossed = SBXMoment(parent1, parent2);
				std::vector<System> crossed = Crossover(parent1, parent2);
				for (std::vector<System>::iterator it = crossed.begin(); it != crossed.end(); ++it) {
					next_gen.push_back(*it);
				}
			}
			else {
				continue;
			}
		}

		models_.clear();
		errors_.clear();

		#pragma omp parallel
		{
			std::vector<System> private_models;
			std::vector<double> private_errors;
			#pragma omp for nowait
			for (int k = 0; k < static_cast<int>(next_gen.size()); ++k) {
				std::vector<Displacement> modeldisp = CalcDisplacementField(next_gen[k], catalog);
				error = TestFitness(modeldisp);
				private_models.push_back(next_gen[k]);
				private_errors.push_back(error);
			}
			#pragma omp critical
			{
				models_.insert(models_.end(), private_models.begin(), private_models.end());
				errors_.insert(errors_.end(), private_errors.begin(), private_errors.end());
			}
		}

//		std::cout << "Minimum error: " << *std::min_element(errors_.begin(), errors_.end()) << "\n";
		std::cout << "Average error: " << std::accumulate(errors_.begin(), errors_.end(), 0.0)/errors_.size() << "\n";
//		std::cout << "Maximum error: " << *std::max_element(errors_.begin(), errors_.end()) << "\n";

		gen_error_.push_back({generation,
							  *std::min_element(errors_.begin(), errors_.end()),
							  std::accumulate(errors_.begin(), errors_.end(), 0.0)/errors_.size(),
							  *std::max_element(errors_.begin(), errors_.end())});

		System best_gen_model;
		double best_gen_error;

		std::tie(best_gen_model, best_gen_error) = GetBestGenModelAndError();

		if (best_gen_error < bestmodelerror_) {
			bestmodel_ = best_gen_model;
			bestmodelerror_ = best_gen_error;
			counter = 0;
		}

		else {
			counter++;
		}

		if (counter > 1000) {
			break;
		}

		for (std::vector<System>::iterator it = models_.begin(); it != models_.end(); ++it) {

			if (mutate(engine) < params_.chance_to_mutate) {
				MomentMutation(*it);
			}
		}
	}

	if (save) {
		std::ofstream ofile(datafilenoext_ + "_error.txt");
		for (std::tuple<int, double, double, double> tup : gen_error_) {
			ofile << (std::get<0>(tup)) << "\t";
			ofile << (std::get<1>(tup)) << "\t";
			ofile << (std::get<2>(tup)) << "\t";
			ofile << (std::get<3>(tup)) << "\n";
		}
		RecordBestModel();
	}
	std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_time = end - start;
	std::cout << elapsed_time.count() << "s\n";
}

double Fitter2D::TestFitness(std::vector<Displacement> modeldisp) {
	double error = 0;
	std::vector<Displacement>::iterator it;

	if (data_.size() != modeldisp.size()) {throw std::runtime_error("Vectors are different lengths");}

	for (int i = 0; i < static_cast<int>(data_.size()); i++) {
//		error += std::abs(data_[i].dx - modeldisp[i].dx);
		error += std::abs(data_[i].dz - modeldisp[i].dz);
	}

	return error;
}

void Fitter2D::RecordModel(System& model) {

	std::vector<Displacement> displacements = CalcDisplacementField(model);
	std::ofstream ofile;
	ofile.open(datafilenoext_ + "_model_data.txt");

	int i = 0;
	for (std::vector<Data>::iterator dit = data_.begin();
										 dit != data_.end(); ++dit, ++i) {
		ofile << dit->pos.x << "\t";
		ofile << dit->pos.y << "\t";
		ofile << displacements[i].dx << "\t";
		ofile << displacements[i].dz << "\n";
	}

	ofile.close();

	ofile.open(datafilenoext_ + "_model_points.txt");
	std::vector<PointSource> sources = model.GetSources();
	for (std::vector<PointSource>::iterator it = sources.begin(); it != sources.end(); ++it) {
		ofile << it->pos.x << "\t";
		ofile << it->pos.y << "\t";
		ofile << it->pos.z << "\t";
		ofile << it->strike << "\t";
		ofile << it->dip << "\t";
		ofile << it->rake << "\t";
		ofile << it->moment << "\n";
	}
	ofile.close();
}

void Fitter2D::RecordFittedData() {

	std::ofstream dout(datafilenoext_ + "_fitted_data.txt");

	for (Data& data : data_) {
		dout << data.pos.x << "\t";
		dout << data.pos.y << "\t";
		dout << data.dx << "\t";
		dout << data.dz << "\n";
	}
}

void Fitter2D::SaveGreensCatalog(const std::string& filename) {
	System model = System(params_);

	GreensCatalog catalog(model, data_);
	catalog.SaveCatalog(filename);
}

void Fitter2D::Setup() {

	System model = System(params_);
	RecordModel(model);
	RecordFittedData();

	GreensCatalog catalog(model, data_);
	catalog.SaveCatalog("greens.txt");
}
