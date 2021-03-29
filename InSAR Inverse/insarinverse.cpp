#include <iomanip>
#include <iostream>
#include <chrono>
#include "insarinverse.h"

const double M_PI = 3.14159265358979323846;

//InSAR::InSAR(const std::string& param_filename, const std::string& data_filename) :
//		param_filename_(param_filename),
//		data_filename_(data_filename) {
//}

//double InSAR::GetParameter(std::string parameter) {
//	std::stringstream counter(parameter);
//	int count = 0;
//	std::string line, buffer;
//	std::ifstream fin(param_filename_);
//	if (!fin) {std::cout << "nope" << "\n";}
//	std::vector<std::string> tokens;
//	while (counter >> buffer) {++count;}
//	while(std::getline(fin, line)) {
//		if (line.find(parameter) != std::string::npos) {
//			std::stringstream ss(line);
//			while(ss >> buffer) {
//				tokens.push_back(buffer);
//			}
//		}
//	}
//	return std::stod(tokens[count + 2]);
//}

//void InSAR::TranslateToASCII () {
//
//	double latstart = GetParameter("Ground Range Data Starting Latitude");
//	double lonstart = GetParameter("Ground Range Data Starting Longitude");
//
//	double latspacing = GetParameter("Ground Range Data Latitude Spacing");
//	double lonspacing = GetParameter("Ground Range Data Longitude Spacing");
//
//	int latlines = GetParameter("Ground Range Data Latitude Lines");
//	int lonsamples = GetParameter("Ground Range Data Longitude Samples");
//
//    float f;
//    std::ifstream fin(data_filename_, std::ios::binary);
//    std::ofstream unwout("unwrapped.txt");
//    if (!fin) {std::cout << "nope" << "\n";}
//    while (fin.read(reinterpret_cast<char*>(&f), sizeof(float))) {
//        unwout << f << "\n";
//    }
//	fin.close();
//	unwout.close();
//
//	std::ofstream locout("locations.txt");
//	locout << std::fixed << std::setprecision(8);
//	for (int j = 0; j < latlines; j++) {
//		for (int i = 0; i < lonsamples; i++) {
//			locout << lonstart + lonspacing*i << "\t";
//			locout << latstart + latspacing*j << "\n";
//		}
//	}
//	locout.close();
//
//}

//void InSAR::Binner() {
//
//	const int latbins = 30;
//	const int lonbins = 30;
//
//	std::array<std::array<double, latbins>, lonbins> hist {0};
//	std::array<std::array<int, latbins>, lonbins> counts {0};
//	std::array<double, latbins + 1> latbounds;
//	std::array<double, lonbins + 1> lonbounds;
//
//	double latstart = GetParameter("Ground Range Data Starting Latitude");
//	double lonstart = GetParameter("Ground Range Data Starting Longitude");
//
//	double latspacing = GetParameter("Ground Range Data Latitude Spacing");
//	double lonspacing = GetParameter("Ground Range Data Longitude Spacing");
//
//	int latlines = GetParameter("Ground Range Data Latitude Lines");
//	int lonsamples = GetParameter("Ground Range Data Longitude Samples");
//
//	double latitude, longitude, phase;
//
//	std::set<double> lats, lons;
//	lats.insert(latstart);
//	lats.insert(latstart + latspacing*latlines);
//	lons.insert(lonstart);
//	lons.insert(lonstart + lonspacing*lonsamples);
//
//	double latmin, latmax, lonmin, lonmax;
//
//	latmin = *lats.begin();
//	latmax = *lats.rbegin();
//	lonmin = *lons.begin();
//	lonmax = *lons.rbegin();
//
//	std::cout << std::fixed << std::setprecision(8);
//	for (int i = 0; i < latbins + 1; i++) {
//		latbounds[i] = latmin + (latmax - latmin)/latbins*i;
//		std::cout << latbounds[i] << "\n";
//	}
//
//	for (int j = 0; j < lonbins + 1; j++) {
//		lonbounds[j] = lonmin + (lonmax - lonmin)/lonbins*j;
//		std::cout << lonbounds[j] << "\n";
//	}
//
//	std::ifstream locin("locations.txt");
//	if (!locin) {std::cout << "nope" << "\n";}
//	std::ifstream unwin("unwrapped.txt");
//	if (!unwin) {std::cout << "nope" << "\n";}
//
//	while (locin >> longitude >> latitude && unwin >> phase) {
//		for (int i = 0; i < latbins; i++) {
//			for (int j = 0; j < lonbins; j++) {
//				if (latitude <= latbounds[i + 1] && longitude <= lonbounds[j + 1]) {
//					hist[j][i] += phase;
//					counts[j][i]++;
//					goto finish;
//				}
//			}
//		}
//		finish:;
//	}
//
//	std::ofstream smout("binned_data.txt");
//	smout << std::fixed << std::setprecision(8);
//	for (int i = 0; i < latbins; i++) {
//		for (int j = 0; j < lonbins; j++) {
//			smout << (lonbounds[j] + lonbounds[j + 1])/2.0 << "\t";
//			smout << (latbounds[i] + latbounds[i + 1])/2.0 << "\t";
//			smout << hist[j][i]/counts[j][i] << "\n";
//		}
//	}
//}

//void InSAR::UnitConverter(const std::string& filename) {
//	double longitude, latitude, phase;
//	double wavelength = GetParameter("Center Wavelength")/100000.0; //km
////			double wavelength = GetParameter("Center Wavelength"); //cm
//	double lonorigin, latorigin;
//	double newx, newy, newz;
//
//	std::ifstream fin(filename);
//	if (!fin) {std::cout << "nope" << "\n";}
//
//	fin >> longitude >> latitude >> phase;
//	lonorigin = longitude;
//	latorigin = latitude;
//
//	std::ofstream fout("converted_data.txt");
//	fout << std::fixed << std::setprecision(8);
//	fout << 0.00000000 << "\t";
//	fout << 0.00000000 << "\t";
//	fout << phase*wavelength/(2.0*M_PI) << "\n";
//
//	while (fin >> longitude >> latitude >> phase) {
//		newx = (longitude - lonorigin)*111.32*std::cos(latitude/180.0*M_PI);
//		newy = (latitude - latorigin)*110.574;
//		newz = phase*wavelength/(2.0*M_PI);
//		fout << newx << "\t";
//		fout << newy << "\t";
//		fout << newz << "\n";
//	}
//	fin.close();
//	fout.close();
//}

void InSAR::UnitConverterNepalQuake(const std::string& filename) {

	size_t lastindex = filename.find_last_of(".");
	std::string file = filename.substr(0, lastindex);

	double longitude, latitude;
	std::string displacement;
	double lonorigin, latorigin;
	double x, y;

	std::ifstream fin(filename);
	if (!fin) {std::cout << "nope" << "\n";}

	fin >> longitude >> latitude >> displacement;
	lonorigin = longitude;
	latorigin = latitude;

	std::ofstream fout(file + "_converted.txt");
	fout << std::fixed << std::setprecision(8);
	fout << 0.00000000 << "\t";
	fout << 0.00000000 << "\t";
	fout << std::stod(displacement)/1.0e6 << "\n";


	while (fin >> longitude >> latitude >> displacement) {
		x = (longitude - lonorigin)*111.32*std::cos(latitude/180.0*M_PI);
		y = (latitude - latorigin)*110.574;
		fout << x << "\t";
		fout << y << "\t";
		fout << std::stod(displacement)/1.0e6 << "\n";
	}
	fin.close();
	fout.close();
}

//void InSAR::BinnerNepalQuake(double latstart, double lonstart, double latspacing, double lonspacing, int latlines, int lonsamples) {
//
//	std::string noval = "--";
//
//	const int latbins = 30;
//	const int lonbins = 30;
//
//	std::array<std::array<double, latbins>, lonbins> hist {0};
//	std::array<std::array<int, latbins>, lonbins> counts {0};
//	std::array<double, latbins + 1> latbounds;
//	std::array<double, lonbins + 1> lonbounds;
//
//	double latitude, longitude;
//	std::string phase;
//
//	std::set<double> lats, lons;
//	lats.insert(latstart);
//	lats.insert(latstart + latspacing*latlines);
//	lons.insert(lonstart);
//	lons.insert(lonstart + lonspacing*lonsamples);
//
//	double latmin, latmax, lonmin, lonmax;
//
//	latmin = *lats.begin();
//	latmax = *lats.rbegin();
//	lonmin = *lons.begin();
//	lonmax = *lons.rbegin();
//
//	std::cout << std::fixed << std::setprecision(8);
//	for (int i = 0; i < latbins + 1; i++) {
//		latbounds[i] = latmin + (latmax - latmin)/latbins*i;
////				std::cout << latbounds[i] << "\n";
//	}
//
//	for (int j = 0; j < lonbins + 1; j++) {
//		lonbounds[j] = lonmin + (lonmax - lonmin)/lonbins*j;
////				std::cout << lonbounds[j] << "\n";
//	}
//
//	std::ifstream datain("testfile.txt");
//	if (!datain) {std::cout << "nope" << "\n";}
//
//	while (datain >> longitude >> latitude >> phase) {
////				if (phase == noval) {phase = "0";}
//		if (phase == noval) {continue;}
//		for (int i = 0; i < latbins; i++) {
//			for (int j = 0; j < lonbins; j++) {
//				if (latitude <= latbounds[i + 1] && longitude <= lonbounds[j + 1]) {
//					hist[j][i] += std::stod(phase);
//					counts[j][i]++;
//					goto finish;
//				}
//			}
//		}
//		finish:;
//	}
//
//	std::ofstream smout("testfile_binned_novalskip.txt");
//	smout << std::fixed << std::setprecision(8);
//	for (int i = 0; i < latbins; i++) {
//		for (int j = 0; j < lonbins; j++) {
//			smout << (lonbounds[j] + lonbounds[j + 1])/2.0 << "\t";
//			smout << (latbounds[i] + latbounds[i + 1])/2.0 << "\t";
//			smout << hist[j][i]/counts[j][i] << "\n";
//			if (counts[j][i] == 0) {
//				std::cout << "zero" << "\n";
//			}
//		}
//	}
//}

void InSAR::BinnerNepalQuake(std::vector<std::string> files) {

	const int latbins = 50; //TODO: Fix so that bins can differ in each direction
	const int lonbins = 50;

	std::array<std::array<double, latbins>, lonbins> hist {0};
	std::array<std::array<int, latbins>, lonbins> counts {0};
	std::array<double, latbins + 1> latbounds;
	std::array<double, lonbins + 1> lonbounds;

	double longitude, latitude, topo, east, north, up, los, sig;

	double lonmin = DBL_MAX;
	double lonmax = -DBL_MAX;
	double latmin = DBL_MAX;
	double latmax = -DBL_MAX;

	for (std::vector<std::string>::iterator it = files.begin(); it != files.end(); ++it) {

		std::ifstream dataread(*it);
		if (!dataread) {std::cout << "nope" << "\n";}

		while (dataread >> longitude >> latitude >> topo >> east >> north >> up >> los >> sig) {
			if (longitude > lonmax) {lonmax = longitude;}
			if (longitude < lonmin) {lonmin = longitude;}
			if (latitude > latmax) {latmax = latitude;}
			if (latitude < latmin) {latmin = latitude;}
		}
	}

	std::cout << std::fixed << std::setprecision(8);
	for (int i = 0; i < latbins + 1; i++) {
		latbounds[i] = latmin + (latmax - latmin)/latbins*i;
//				std::cout << latbounds[i] << "\n";
	}

	for (int j = 0; j < lonbins + 1; j++) {
		lonbounds[j] = lonmin + (lonmax - lonmin)/lonbins*j;
//				std::cout << lonbounds[j] << "\n";
	}

	for (std::vector<std::string>::iterator it = files.begin(); it != files.end(); ++it) {

		std::ifstream datain(*it);
		if (!datain) {std::cout << "nope" << "\n";}

		while (datain >> longitude >> latitude >> topo >> east >> north >> up >> los >> sig) {
			for (int i = 0; i < latbins; i++) {
				for (int j = 0; j < lonbins; j++) {
					if (latitude <= latbounds[i + 1] && longitude <= lonbounds[j + 1]) {
						hist[j][i] += los*up;
//								hist[j][i] += los;
						counts[j][i]++;
						goto finish;
					}
				}
			}
			finish:;
		}
	}

	std::ofstream smout("total_displacement_binned_50x50.txt");
	smout << std::fixed << std::setprecision(8);
	for (int i = 0; i < latbins; i++) {
		for (int j = 0; j < lonbins; j++) {
			smout << (lonbounds[j] + lonbounds[j + 1])/2.0 << "\t";
			smout << (latbounds[i] + latbounds[i + 1])/2.0 << "\t";
			smout << hist[j][i]/counts[j][i] << "\n";
			std::cout << hist[j][i] << " : " << counts[j][i] << "\n";
			if (counts[j][i] == 0) {
				std::cout << "zero" << "\n";
			}
		}
	}
}

//void InSAR::GeneratePointsNepalQuake(const int xpoints, const int ypoints,
//		                                    const double xmin, const double xmax,
//											const double ymin, const double ymax) {
//
//	std::array<double, 2> topleft, topright, botleft, botright;
//	double x, y, depth;
//
//	std::ofstream fout("testfile_binned_converted_sourcepoints.txt");
//	if (!fout) {std::cout << "nope" << "\n";}
//
//	topleft[0] = xmin;
//	topleft[1] = ymax;
//	topright[0] = xmax;
//	topright[1] = ymax;
//	botleft[0] = xmin;
//	botleft[1] = ymin;
//	botright[0] = xmax;
//	botright[1] = ymin;
//
//	double dx = (xmax - xmin)/xpoints;
//	double dy = (ymax - ymin)/ypoints;
//
//	for (depth = 2; depth <= 2; ++ depth) {
//		for (int i = 0; i < ypoints; ++i) {
//			for (int j = 0; j < xpoints; ++j) {
//				x = xmin + dx/2 + dx*j;
//				y = ymin + dy/2 + dy*i;
//				fout << x << "\t" << y << "\t" << -depth << "\n";
//			}
//		}
//	}
//}

//void InSAR::GeneratePoints() {
//	std::array<double, 2> topleft, topright, botleft, botright;
//	double lonorigin, latorigin, longitude, latitude, phase;
//	double leftslope, leftintercept, botslope, botintercept, x, y;
//	int leftpoints = 20;
//	int botpoints = 10;
//	double depth; //km
//
//	std::ifstream fin("binned_data.txt");
//	if (!fin) {std::cout << "nope" << "\n";}
//	fin >> longitude >> latitude >> phase;
//	lonorigin = longitude;
//	latorigin = latitude;
//
//	std::ofstream fout("sourcepoints.txt");
//	if (!fout) {std::cout << "nope" << "\n";}
//
//	topleft[0] = GetParameter("Approximate Upper Left Longitude");
//	topleft[1] = GetParameter("Approximate Upper Left Latitude");
//	topright[0] = GetParameter("Approximate Upper Right Longitude");
//	topright[1] = GetParameter("Approximate Upper Right Latitude");
//	botleft[0] = GetParameter("Approximate Lower Left Longitude");
//	botleft[1] = GetParameter("Approximate Lower Left Latitude");
//	botright[0] = GetParameter("Approximate Lower Right Longitude");
//	botright[1] = GetParameter("Approximate Lower Right Latitude");
//
//	topleft[0] = (topleft[0] - lonorigin)*111.32*std::cos(topleft[1]/180.0*M_PI);
//	topleft[1] = (topleft[1] - latorigin)*110.574;
//	topright[0] = (topright[0] - lonorigin)*111.32*std::cos(topright[1]/180.0*M_PI);
//	topright[1] = (topright[1] - latorigin)*110.574;
//	botleft[0] = (botleft[0] - lonorigin)*111.32*std::cos(botleft[1]/180.0*M_PI);
//	botleft[1] = (botleft[1] - latorigin)*110.574;
//	botright[0] = (botright[0] - lonorigin)*111.32*std::cos(botright[1]/180.0*M_PI);
//	botright[1] = (botright[1] - latorigin)*110.574;
//
//	leftslope = (topleft[1] - botleft[1])/(topleft[0] - botleft[0]);
//	leftintercept = botleft[1] - leftslope*botleft[0];
//
//	botslope = (botright[1] - botleft[1])/(botright[0] - botleft[0]);
//	botintercept = botleft[1] - botslope*botleft[0];
//
//	for (depth = 2; depth <= 4; ++ depth) {
//		for (int i = 0; i < botpoints; ++i) {
//			x = botleft[0] + i*(botright[0] - botleft[0])/(botpoints - 1);
//			y = botslope*x + botintercept;
//			fout << x << "\t" << y << "\t" << depth << "\n";
//			for (int j = 1; j < leftpoints; ++j) {
//				leftintercept = y - leftslope*x;
//				x += (topleft[0] - botleft[0])/(leftpoints - 1);
//				y = leftslope*x + leftintercept;
//				fout << x << "\t" << y << "\t" << depth << "\n";
//			}
//		}
//	}
//}

PointSource::PointSource (const double x, const double y, const double z,
						  const double strike, const double dip,
						  const double rake, const double M0) :
		x_coord_(x),
		y_coord_(y),
		z_coord_(z),
		strike_angle_(strike),
		dip_angle_(dip),
		rake_angle_(rake),
		seismic_moment_(M0) {

}

PointSource PointSource::CopySource() {
	return PointSource(this->GetX(), this->GetY(), this->GetZ(), this->GetStrikeAngle(), this->GetDipAngle(), this->GetRakeAngle(), this->GetSeismicMoment());
}

System::System() {

}

void System::CopySources(System sys) { //copy the sources from sys into the calling object
	this->ClearSources();
	for (std::vector<PointSource>::iterator it = sys.sources_.begin(); it != sys.sources_.end(); ++it) {
		this->sources_.push_back(it->CopySource());
	}
}

Parameters::Parameters() {

}

Parameters::Parameters(const std::string& param_filename) {

	ParameterReader params(param_filename);

	xmin_ = std::stod(params.GetParameter("Minimum x value"));
	xmax_ = std::stod(params.GetParameter("Maximum x value"));
	ymin_ = std::stod(params.GetParameter("Minimum y value"));
	ymax_ = std::stod(params.GetParameter("Maximum y value"));
	zmin_ = std::stod(params.GetParameter("Minimum z value"));
	zmax_ = std::stod(params.GetParameter("Maximum z value"));
	min_moment_ = std::stod(params.GetParameter("Minimum seismic moment"));
	max_moment_ = std::stod(params.GetParameter("Maximum seismic moment"));
	chance_to_mutate_ = std::stod(params.GetParameter("Chance to mutate"));

	xdatapoints_ = std::stoi(params.GetParameter("Num x data points"));
	ydatapoints_ = std::stoi(params.GetParameter("Num y data points"));
	gen_num_sources_ = std::stoi(params.GetParameter("Num sources (gen)"));
	fit_num_sources_ = std::stoi(params.GetParameter("Num sources (fit)"));
	num_generations_ = std::stoi(params.GetParameter("Number of generations"));
	pop_size_ = std::stoi(params.GetParameter("Population size"));

	source_type_ = params.GetParameter("Source type");

}

std::mt19937 DataGenerator::engine(time(0));

DataGenerator::DataGenerator(const std::string& param_filename, std::string const int_filename) {

	filename_ = int_filename;

	size_t lastindex = int_filename.find_last_of(".");
	filename_noext_ = int_filename.substr(0, lastindex);

	params_ = Parameters(param_filename);

}

std::vector<std::array<double, 3> > DataGenerator::Interferogram() {

	double x, y, datadx, datady;
	std::vector<std::array<double, 3> > interferogram;
	std::vector<PointSource> sources;

	double xmin = params_.xmin_;
	double xmax = params_.xmax_;
	double ymin = params_.ymin_;
	double ymax = params_.ymax_;
	double zmin = params_.zmin_;
	double zmax = params_.zmax_;
	double moment_order_min = std::log10(params_.min_moment_);
	double moment_order_max = std::log10(params_.max_moment_);

	int xdatapoints = params_.xdatapoints_;
	int ydatapoints = params_.ydatapoints_;
	int numsources = params_.gen_num_sources_;

	std::string source_type = params_.source_type_;

	std::uniform_real_distribution<double> dist_x(xmin, xmax);
	std::uniform_real_distribution<double> dist_y(ymin, ymax);
	std::uniform_real_distribution<double> dist_z(zmin, zmax);
	std::uniform_real_distribution<double> strike(0.0, 2.0*M_PI);
	std::uniform_real_distribution<double> dip(0.0, M_PI/2.0);
	std::uniform_real_distribution<double> order(moment_order_min,
												 moment_order_max);

	datadx = (xmax - xmin)/xdatapoints;
	datady = (ymax - ymin)/ydatapoints;

	for (int i = 0; i < numsources; ++i) {
		sources.push_back(PointSource(dist_x(engine), dist_y(engine), -dist_z(engine), strike(engine), dip(engine), 0.0, pow(10, order(engine))));
	}


	for (y = ymin + datady/2; y < ymax; y += datady) {
		for (x = xmin + datadx/2; x < xmax; x += datadx) {
			interferogram.push_back({x, y, 0.0});
		}
	}

	if (source_type == "strike") {
		for (std::array<double, 3>& datapoint: interferogram) {
			for (PointSource& source : sources) {
				datapoint[2] += strike_disp_z(datapoint[0], datapoint[1], 0, source.GetX(), source.GetY(), -source.GetZ(), source.GetStrikeAngle(), source.GetDipAngle(), source.GetSeismicMoment());
			}
		}
	}

	else if (source_type == "dip") {
		for (std::array<double, 3>& datapoint: interferogram) {
			for (PointSource& source : sources) {
				datapoint[2] += dip_disp_z(datapoint[0], datapoint[1], 0, source.GetX(), source.GetY(), -source.GetZ(), source.GetStrikeAngle(), source.GetDipAngle(), source.GetSeismicMoment());
			}
		}
	}

	else {
		std::cout << "Source type not recognized or set" << std::endl;
		throw;
	}

	std::ofstream pout(filename_noext_ + "_sources.txt");
	if (!pout) {std::cout << "Couldn't open file" << "\n";}
	for (PointSource& source : sources) {
		pout << source.GetX() << "\t";
		pout << source.GetY() << "\t";
		pout << source.GetZ() << "\t";
		pout << source.GetStrikeAngle() << "\t";
		pout << source.GetDipAngle() << "\t";
		pout << source.GetSeismicMoment() << "\n";
	}
	pout.close();

	std::ofstream iout(filename_);
	if (!iout) {std::cout << "Couldn't open file" << "\n";}
	for (std::array<double, 3>& datapoint : interferogram) {
		iout << datapoint[0] << "\t";
		iout << datapoint[1] << "\t";
		iout << datapoint[2] << "\n";
	}
	iout.close();

	return interferogram;
}

std::mt19937 Fitter::engine(time(0));

ParameterReader::ParameterReader(const std::string& param_filename) {
	std::ifstream file(param_filename);
	std::string line;
	std::vector<std::string> tokens;

	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string token;
		while (std::getline(iss, token, '\t')) {
			if (!token.empty() && token != "=") {
				tokens.push_back(token);
			}
		}

		std::string key = tokens[0];
		std::string value = tokens[1];

		params_.insert({key, value});
		tokens.clear();
	}
}

std::string ParameterReader::GetParameter(const std::string& param_name) {

	if (params_.find(param_name) != params_.end()) {
		return params_[param_name];
	}

	else {
		std::cout << "Parameter \"" << param_name << "\" not found" << std::endl;
		throw;
	}
}

DataReader::DataReader(const std::string& data_filename) {

	size_t lastindex = data_filename.find_last_of(".");
	datafilenoext_ = data_filename.substr(0, lastindex);

	std::ifstream fin(data_filename);
	if (!fin) {std::cout << "nope" << "\n";}
	double x, y, z;
	std::string zin;
	std::array<double, 3> triplet;

	while (fin >> x >> y >> zin) {
		z = std::stod(zin);
		triplet[0] = x;
		triplet[1] = y;
		triplet[2] = z;

		if (!std::isnan(triplet[2])) {
			data_.push_back(triplet);
	    }
	}
}

std::vector<std::array<double, 3> > DataReader::Subset(const double xmin,
													   const double xmax,
													   const double ymin,
													   const double ymax) {

	double x, y;
	std::vector<std::array<double, 3> > subset;

	for (std::array<double, 3> triplet : data_) {
		x = triplet[0];
		y = triplet[1];

		if (x <= xmax && x >= xmin && y <= ymax && y >= ymin) {
			subset.push_back(triplet);
	    }
	}

	return subset;
}

std::vector<double> DataReader::FindBounds() {

	double xmin = DBL_MAX;
	double xmax = -DBL_MAX;
	double ymin = DBL_MAX;
	double ymax = -DBL_MAX;

	for (std::array<double, 3> triplet : data_) {
		if (triplet[0] > xmax) {
			xmax = triplet[0];
		}
		if (triplet[0] < xmin) {
			xmin = triplet[0];
		}
		if (triplet[1] < ymin) {
			ymin = triplet[1];
		}
		if (triplet[1] > ymax) {
			ymax = triplet[1];
		}
	}

	return {xmin, xmax, ymin, ymax};
}

Fitter::Fitter() {

}

Fitter::Fitter(const std::string& param_filename, const std::string& data_filename) {

	size_t lastindex = data_filename.find_last_of(".");
	datafilenoext_ = data_filename.substr(0, lastindex);

	params_ = Parameters(param_filename);

	DataReader data(data_filename);
	data_ = data.Subset(params_.xmin_, params_.xmax_,
						params_.ymin_, params_.ymax_);
}

//void Fitter::ImportParameters(const std::string& param_filename) {
//
//	std::ifstream file(param_filename);
//	std::string line;
//	std::vector<std::string> tokens;
//	std::map<std::string, std::string> params;
//
//	while (std::getline(file, line)) {
//		std::istringstream iss(line);
//		std::string token;
//		while (std::getline(iss, token, '\t')) {
//			if (!token.empty() && token != "=")
//				tokens.push_back(token);
//		}
//		params.insert({tokens[0], tokens[1]});
//		tokens.clear();
//	}
//
//	xmin_ = std::stod(params.at("Minimum x value"));
//	xmax_ = std::stod(params.at("Maximum x value"));
//	ymin_ = std::stod(params.at("Minimum y value"));
//	ymax_ = std::stod(params.at("Maximum y value"));
//	zmin_ = std::stod(params.at("Minimum z value"));
//	zmax_ = std::stod(params.at("Maximum z value"));
//	min_moment_ = std::stod(params.at("Minimum seismic moment"));
//	max_moment_ = std::stod(params.at("Maximum seismic moment"));
//	chance_to_mutate_ = std::stod(params.at("Chance to mutate"));
//	num_generations_ = std::stoi(params.at("Number of generations"));
//	num_sources_ = std::stoi(params.at("Number of sources"));
//	pop_size_ = std::stoi(params.at("Population size"));
//	source_type_ = params.at("Source type");
//}

//void Fitter::ImportData(const std::string& data_filename) {
//
//	size_t lastindex = data_filename.find_last_of(".");
//	datafilenoext_ = data_filename.substr(0, lastindex);
//
//	std::ifstream fin(data_filename);
//	if (!fin) {std::cout << "nope" << "\n";}
//	double x, y, z;
//	std::string zin;
//	std::array<double, 3> triplet;
//
//	while (fin >> x >> y >> zin) {
//		z = std::stod(zin);
//		triplet[0] = x;
//		triplet[1] = y;
//		triplet[2] = z;
//
//		if (!std::isnan(triplet[2]) && x <= xmax_ && x >= xmin_ && y <= ymax_ && y >= ymin_) {
//			data_.push_back(triplet);
//	    }
//	}
//}

//int Fitter::GetDataXMax() {
//	int x_max = -INT_MAX;
//	for (std::vector<std::array<double, 3> >::iterator it = data_.begin(); it != data_.end(); ++it) {
//		if (*(it->data()) > x_max) {
//			x_max = *(it->data());
//		}
//	}
//	return x_max;
//}
//
//int Fitter::GetDataXMin() {
//	int x_min = INT_MAX;
//	for (std::vector<std::array<double, 3> >::iterator it = data_.begin(); it != data_.end(); ++it) {
//		if (*(it->data()) < x_min) {
//			x_min = *(it->data());
//		}
//	}
//	return x_min;
//}
//
//int Fitter::GetDataYMax() {
//	int y_max = -INT_MAX;
//	for (std::vector<std::array<double, 3> >::iterator it = data_.begin(); it != data_.end(); ++it) {
//		if (*(it->data() + 1) > y_max) {
//			y_max = *(it->data() + 1);
//		}
//	}
//	return y_max;
//}
//
//int Fitter::GetDataYMin() {
//	int y_min = INT_MAX;
//	for (std::vector<std::array<double, 3> >::iterator it = data_.begin(); it != data_.end(); ++it) {
//		if (*(it->data() + 1) < y_min) {
//			y_min = *(it->data() + 1);
//		}
//	}
//	return y_min;
//}

std::vector<double> Fitter::GetDataDisplacementField() {
	std::vector<double> displacements;
	for (std::vector<std::array<double, 3> >::iterator it = data_.begin(); it != data_.end(); ++it) {
		displacements.push_back(*(it->data() + 2));
	}
	return displacements;
}

//std::vector<double> Fitter::GetDataUncertainty() {
//	std::vector<double> uncertainty;
//	for (std::vector<std::array<double, 3> >::iterator it = data_.begin(); it != data_.end(); ++it) {
//		uncertainty.push_back(*(it->data() + 3));
//	}
//	return uncertainty;
//}

std::vector<double> Fitter::CalcDisplacementField(System model) {
	std::vector<double> displacements(data_.size(), 0.0);
	std::vector<PointSource> sources = model.GetSources();
	for (std::vector<PointSource>::iterator it = sources.begin(); it != sources.end(); ++it) {
		int i = 0;
		for (std::vector<std::array<double, 3> >::iterator dit = data_.begin(); dit != data_.end(); ++dit, ++i) {
			if (params_.source_type_ == "strike") {
				displacements[i] += strike_disp_z(*(dit->data()), *(dit->data() + 1), 0, it->GetX(), it->GetY(), fabs(it->GetZ()), it->GetStrikeAngle(), it->GetDipAngle(), it->GetSeismicMoment());
			}

			else if (params_.source_type_ == "dip") {
				displacements[i] += dip_disp_z(*(dit->data()), *(dit->data() + 1), 0, it->GetX(), it->GetY(), fabs(it->GetZ()), it->GetStrikeAngle(), it->GetDipAngle(), it->GetSeismicMoment());
			}

			else {
				std::cout << "Source type not recognized or not set" << "\n";
				throw;
			}
		}
	}
	return displacements;
}

std::vector<System> Fitter::Crossover(System parent1, System parent2) {

	std::vector<PointSource> parent1_sources = parent1.GetSources();
	std::vector<PointSource> parent2_sources = parent2.GetSources();

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

std::vector<System> Fitter::SimulatedBinaryCrossover(System parent1, System parent2) {

	double beta, eta_c, u;
	double p1strike, p2strike, p1dip, p2dip, p1mom, p2mom;
	double p1strike_new, p2strike_new, p1dip_new, p2dip_new, p1mom_new, p2mom_new;
	double p1x, p1y, p1z, p2x, p2y, p2z;
	double p1x_new, p1y_new, p1z_new, p2x_new, p2y_new, p2z_new;

	std::vector<PointSource> parent1_sources = parent1.GetSources();
	std::vector<PointSource> parent2_sources = parent2.GetSources();
	if (parent1_sources.size() != parent2_sources.size()) {
		std::cout << "Parents don't have the same number of sources" << std::endl;
		throw;
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
		p1x = pit->GetX();
		p2x = sit->GetX();
		p1y = pit->GetY();
		p2y = sit->GetY();
		p1z = pit->GetZ();
		p2z = pit->GetZ();

		p1strike = pit->GetStrikeAngle();
		p2strike = sit->GetStrikeAngle();
		p1dip = pit->GetDipAngle();
		p2dip = sit->GetDipAngle();
		p1mom = pit->GetSeismicMoment();
		p2mom = sit->GetSeismicMoment();

		p1x_new = 0.5*((1 + beta)*p1x + (1 - beta)*p2x);
		p2x_new = 0.5*((1 - beta)*p1x + (1 + beta)*p2x);

		p1y_new = 0.5*((1 + beta)*p1y + (1 - beta)*p2y);
		p2y_new = 0.5*((1 - beta)*p1y + (1 + beta)*p2y);

		p1z_new = 0.5*((1 + beta)*p1z + (1 - beta)*p2z);
		p2z_new = 0.5*((1 - beta)*p1z + (1 + beta)*p2z);

		if (p1z_new > 0) {
			p1z_new = -p1z_new;
		}

		if (p2z_new > 0) {
			p2z_new = -p2z_new;
		}

		p1strike_new = 0.5*((1 + beta)*p1strike + (1 - beta)*p2strike);
		p2strike_new = 0.5*((1 - beta)*p1strike + (1 + beta)*p2strike);

		p1dip_new = 0.5*((1 + beta)*p1dip + (1 - beta)*p2dip);
		p2dip_new = 0.5*((1 - beta)*p1dip + (1 + beta)*p2dip);

		p1mom_new = 0.5*((1 + beta)*p1mom + (1 - beta)*p2mom);
		p2mom_new = 0.5*((1 - beta)*p1mom + (1 + beta)*p2mom);

		pit->SetX(p1x_new);
		sit->SetX(p2x_new);
		pit->SetY(p1y_new);
		sit->SetY(p2y_new);
		pit->SetZ(p1z_new);
		sit->SetZ(p2z_new);

		pit->SetStrikeAngle(p1strike_new);
		sit->SetStrikeAngle(p2strike_new);
		pit->SetDipAngle(p1dip_new);
		sit->SetDipAngle(p2dip_new);
		pit->SetSeismicMoment(p1mom_new);
		sit->SetSeismicMoment(p2mom_new);
	}

	new_sys.push_back(parent1);
	new_sys.push_back(parent2);

	return new_sys;

}

void Fitter::LocationMutation (System &sys) {
	std::vector<PointSource>& sources = sys.GetSources();
	std::uniform_real_distribution<double> zerotoone(0.0, 1.0);
	double x, y, z;

	for (std::vector<PointSource>::iterator it = sources.begin(); it != sources.end(); ++it) {
		x = it->GetX();
		y = it->GetY();
		z = it->GetZ();
		std::normal_distribution<double> x_spread(x, 2);
		std::normal_distribution<double> y_spread(y, 2);
		std::normal_distribution<double> z_spread(z, 0.5);

		if (zerotoone(engine) < 0.1) { //TODO: Mess around with this chance
			double new_x = x_spread(engine);
			while (new_x > params_.xmax_ || new_x < params_.xmin_) {
				new_x = x_spread(engine);
			}
			double new_y = y_spread(engine);
			while (new_y > params_.ymax_ || new_y < params_.ymin_) {
				new_y = y_spread(engine);
			}
			double new_z = z_spread(engine);
			while (new_z > 0) {
				new_z = z_spread(engine);
			}
			it->SetX(new_x);
			it->SetY(new_y);
			it->SetZ(new_z);
		}
	}
}

void Fitter::OrientationMutation (System &sys) {
	std::vector<PointSource>& sources = sys.GetSources();
	std::uniform_real_distribution<double> zerotoone(0.0, 1.0);
	double strike, dip, moment;

	for (std::vector<PointSource>::iterator it = sources.begin(); it != sources.end(); ++it) {
		strike = it->GetStrikeAngle();
		dip = it->GetDipAngle();
		moment = it->GetSeismicMoment();
		std::normal_distribution<double> strike_spread(strike, M_PI/6.0);
		std::normal_distribution<double> dip_spread(dip, M_PI/24.0);
		std::normal_distribution<double> moment_spread(std::log10(moment), 0.5);

		if (zerotoone(engine) < 0.1) { //TODO: Mess around with this chance
			double new_strike = strike_spread(engine);
			double new_dip = dip_spread(engine);
			while (new_dip > M_PI/180*90 || new_dip < 0) {
				new_dip = dip_spread(engine);
			}
			double new_moment = pow(10.0, moment_spread(engine));
			it->SetStrikeAngle(new_strike);
			it->SetDipAngle(new_dip);
			it->SetSeismicMoment(new_moment);
		}
	}
}

System Fitter::SelectParent() {
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

void Fitter::GeneticAlgorithm(bool save) {

	models_.clear();
	errors_.clear();
	bestmodelerror_ = DBL_MAX;
	gen_error_.clear();

	std::chrono::high_resolution_clock::time_point start =
			std::chrono::high_resolution_clock::now();

	std::uniform_int_distribution<int> dist_x(params_.xmin_, params_.xmax_);
	std::uniform_int_distribution<int> dist_y(params_.ymin_, params_.ymax_);
	std::uniform_int_distribution<int> dist_z(params_.zmin_, params_.zmax_);
	std::uniform_real_distribution<double> strike(0.0, 2.0*M_PI);
	std::uniform_real_distribution<double> dip(0.0, M_PI/2.0);
	std::uniform_real_distribution<double> moment(params_.min_moment_,
												  params_.max_moment_);
	std::uniform_real_distribution<double> mutate(0.0, 1.0);

	std::vector<std::array<double, 3> > data = data_;

	System model;
	System bestmodel;

	std::vector<System> next_gen;

	double error;

	std::vector<double> datadisp = this->GetDataDisplacementField();

//	if (!initial_state.empty()) {
//		double x_pos, y_pos, z_pos, strike_angle, dip_angle, rake_angle, seis_moment;
//		std::vector<std::array<double, 7> > initial_points;
//		std::ifstream incon(initial_state);
//		while (incon >> x_pos >> y_pos >> z_pos >> strike_angle >> dip_angle
//			   >> rake_angle >> seis_moment) {
//			initial_points.push_back({x_pos, y_pos, z_pos, strike_angle,
//				                      dip_angle, rake_angle, seis_moment});
//		}
//
//		std::normal_distribution<double> xy_spread(0, 2);
//		std::normal_distribution<double> z_spread(0, 0.5);
//		std::normal_distribution<double> strike_spread(0, M_PI/6.0);
//		std::normal_distribution<double> dip_spread(0, M_PI/24.0);
////		std::normal_distribution rake_spread(0, 0);
////		std::normal_distribution seism_spread(0, 1.0e8);
//
//		for (int i = 0; i < params_.pop_size_; i++) { //generate random first generation parents
//
//			model.ClearSources();
//
//			for (std::array<double, 7> point : initial_points) { //number of sources in each parent
//				x_pos = point[0] + xy_spread(engine);
//				y_pos = point[1] + xy_spread(engine);
//				z_pos = point[2] + z_spread(engine);
//				strike_angle = point[3] + strike_spread(engine);
//				dip_angle = point[4] + dip_spread(engine);
////				rake_angle = point[5] + rake_spread(engine);
////				seis_moment = point[6] + seism_spread(engine);
//				model.PlaceSource(x_pos, y_pos, z_pos, strike_angle,
//	                      	  	  dip_angle, rake_angle, seis_moment);
//			}
//
//			std::vector<double> modeldisp = CalcDisplacementField(model);
//
//			error = TestFitness(datadisp, modeldisp);
//
//			models_.push_back(model);
//			errors_.push_back(error);
//
//		}
//	}


	for (int i = 0; i < params_.pop_size_; i++) { //generate random first generation parents

		model.ClearSources();

		for (int j = 0; j < params_.fit_num_sources_; j++) { //number of sources in each parent
			model.PlaceSource(dist_x(engine), dist_y(engine), -dist_z(engine),
							  strike(engine), dip(engine), 0, moment(engine));
		}

		std::vector<double> modeldisp = CalcDisplacementField(model); //TODO: Must update to include rake angle

		error = TestFitness(datadisp, modeldisp);

		models_.push_back(model);
		errors_.push_back(error);

	}


	int counter = 0;
	for (int generation = 0; generation < params_.num_generations_; generation++) {

		std::cout << "Generation " << generation << std::endl;

		for (std::vector<System>::iterator it = models_.begin(); it != models_.end(); ++it) {
			it->SetIndex(std::distance(models_.begin(), it));
		}

		next_gen.clear();

		while (static_cast<int>(next_gen.size()) < params_.pop_size_) {
			System parent1 = SelectParent();
			System parent2 = SelectParent();
			if (parent1.GetIndex() != parent2.GetIndex()) {
//				std::vector<System> crossed = Crossover(parent1, parent2);
				std::vector<System> crossed = SimulatedBinaryCrossover(parent1, parent2);
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
				std::vector<double> modeldisp = CalcDisplacementField(next_gen[k]);
				error = TestFitness(datadisp, modeldisp);
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

			if (mutate(engine) < params_.chance_to_mutate_) {
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
	std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_time = end - start;
	std::cout << elapsed_time.count() << "s\n";
}

//void Fitter::BothGeneticAlgorithm(const int generations, std::string initial_state) {
//
//	std::chrono::high_resolution_clock::time_point start =
//			std::chrono::high_resolution_clock::now();
//
//	std::uniform_int_distribution<int> dist_x(GetDataXMin(), GetDataXMax());
//	std::uniform_int_distribution<int> dist_y(GetDataYMin(), GetDataYMax());
//	std::uniform_int_distribution<int> dist_z(3.0, 10.0);
//	std::uniform_real_distribution<double> strike(0.0, 2.0*M_PI);
//	std::uniform_real_distribution<double> dip(0.0, M_PI/2.0);
//	std::uniform_real_distribution<double> moment(1.0e9, 1.0e12);
//	std::uniform_real_distribution<double> mutate(0.0, 1.0);
//	std::vector<std::array<double, 3> > data = data_;
//
//	System model;
//	System bestmodel;
//
//	std::vector<System> next_gen;
//	int source_count = 1;
//	int population_size = 500;
//	int number_of_generations = generations;
//
//	double error;
//
//	std::vector<double> datadisp = this->GetDataDisplacementField();
//
//	if (!initial_state.empty()) {
//		double x_pos, y_pos, z_pos, strike_angle, dip_angle, rake_angle, seis_moment;
//		std::vector<std::array<double, 7> > initial_points;
//		std::ifstream incon(initial_state);
//		while (incon >> x_pos >> y_pos >> z_pos >> strike_angle >> dip_angle
//			   >> rake_angle >> seis_moment) {
//			initial_points.push_back({x_pos, y_pos, z_pos, strike_angle,
//				                      dip_angle, rake_angle, seis_moment});
//		}
//
//		std::normal_distribution<double> xy_spread(0, 2);
//		std::normal_distribution<double> z_spread(0, 0.5);
//		std::normal_distribution<double> strike_spread(0, M_PI/6.0);
//		std::normal_distribution<double> dip_spread(0, M_PI/24.0);
////		std::normal_distribution rake_spread(0, 0);
////		std::normal_distribution seism_spread(0, 1.0e8);
//
//		for (int i = 0; i < population_size; i++) { //generate random first generation parents
//
//			model.ClearSources();
//
//			for (std::array<double, 7> point : initial_points) { //number of sources in each parent
//				x_pos = point[0] + xy_spread(engine);
//				y_pos = point[1] + xy_spread(engine);
//				z_pos = point[2] + z_spread(engine);
//				strike_angle = point[3] + strike_spread(engine);
//				dip_angle = point[4] + dip_spread(engine);
////				rake_angle = point[5] + rake_spread(engine);
////				seis_moment = point[6] + seism_spread(engine);
//				model.PlaceSource(x_pos, y_pos, z_pos, strike_angle,
//	                      	  	  dip_angle, rake_angle, seis_moment);
//			}
//
//			std::vector<double> modeldisp = CalcDisplacementField(model); //TODO: Must update to include rake angle
//
//			error = TestFitness(datadisp, modeldisp);
//
//			models_.push_back(model);
//			errors_.push_back(error);
//
//		}
//	}
//
//	else {
//		for (int i = 0; i < population_size; i++) { //generate random first generation parents
//
//			model.ClearSources();
//
//			for (int j = 0; j < source_count; j++) { //number of sources in each parent
//				model.PlaceSource(dist_x(engine), dist_y(engine), -dist_z(engine),
//								  strike(engine), dip(engine), 0, moment(engine));
//			}
//
//			std::vector<double> modeldisp = CalcDisplacementField(model); //TODO: Must update to include rake angle
//
//			error = TestFitness(datadisp, modeldisp);
//
//			models_.push_back(model);
//			errors_.push_back(error);
//
//		}
//	}
//
//	std::ofstream ofile;
//	ofile.open(datafile_ + "_error.txt");
//
//	for (int generation = 0; generation < number_of_generations; generation++) {
//
//		std::cout << "Generation " << generation << std::endl;
//
//		for (std::vector<System>::iterator it = models_.begin(); it != models_.end(); ++it) {
//			it->SetIndex(std::distance(models_.begin(), it));
//		}
//
//		next_gen.clear();
//
//		while (static_cast<int>(next_gen.size()) < population_size) {
//			System parent1 = SelectParent();
//			System parent2 = SelectParent();
//			if (parent1.GetIndex() != parent2.GetIndex()) {
////				std::vector<System> crossed = Crossover(parent1, parent2); //TODO: Testing new crossover
//				std::vector<System> crossed = SimulatedBinaryCrossover(parent1, parent2);
//				for (std::vector<System>::iterator it = crossed.begin(); it != crossed.end(); ++it) {
//					next_gen.push_back(*it);
//				}
//			}
//			else {
//				continue;
//			}
//		}
//
//		models_.clear();
//		errors_.clear();
//
//		#pragma omp parallel
//		{
//			std::vector<System> private_models;
//			std::vector<double> private_errors;
//			#pragma omp for nowait
//			for (int k = 0; k < static_cast<int>(next_gen.size()); ++k) {
//				std::vector<double> modeldisp = CalcDisplacementField(next_gen[k]);
//				error = TestFitness(datadisp, modeldisp);
//				private_models.push_back(next_gen[k]);
//				private_errors.push_back(error);
//			}
//			#pragma omp critical
//			{
//				models_.insert(models_.end(), private_models.begin(), private_models.end());
//				errors_.insert(errors_.end(), private_errors.begin(), private_errors.end());
//			}
//		}
//
//		//std::cout << "Minimum error: " << *std::min_element(errors_.begin(), errors_.end()) << "\n";
//		std::cout << "Average error: " << std::accumulate(errors_.begin(), errors_.end(), 0.0)/errors_.size() << "\n";
//		//std::cout << "Maximum error: " << *std::max_element(errors_.begin(), errors_.end()) << "\n";
//
//		ofile << generation << "\t";
//		ofile << *std::min_element(errors_.begin(), errors_.end()) << "\t";
//		ofile << std::accumulate(errors_.begin(), errors_.end(), 0.0)/errors_.size() << "\t";
//		ofile << *std::max_element(errors_.begin(), errors_.end()) << "\n";
//
//		for (std::vector<System>::iterator it = models_.begin(); it != models_.end(); ++it) {
//
//			if (mutate(engine) < 0.2) { //TODO: Mess around with this chance
//				OrientationMutation(*it);
//				LocationMutation(*it);
////				AdaptiveOrientationMutation(*it, generation, number_of_generations);
//			}
//		}
//	}
//	ofile.close();
//	RecordBestModel();
//	std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
//	std::chrono::duration<double> elapsed_time = end - start;
//	std::cout << elapsed_time.count() << "s\n";
//}

double Fitter::TestFitness(std::vector<double> datadisp, std::vector<double> modeldisp) {
	double error;
	double temp = 0;
	double n = 2.0;

	if (datadisp.size() != modeldisp.size()) {throw "Vectors are different lengths";}

	for (int i = 0; i < static_cast<int>(datadisp.size()); i++) {
		temp += pow(datadisp[i] - modeldisp[i], n);
	}

	error = pow(temp, 1/n);
	return error;
}

void Fitter::RecordModel(System model) {

	std::vector<double> displacements = CalcDisplacementField(model);
	std::ofstream ofile;
	ofile.open(datafilenoext_ + "_model_data.txt");

	int i = 0;
	for (std::vector<std::array<double, 3> >::iterator dit = data_.begin();
													   dit != data_.end();
													   ++dit, ++i) {
		ofile << *(dit->data()) << "\t";
		ofile << *(dit->data() + 1) << "\t";
		ofile << displacements[i] << "\n";
	}

	ofile.close();

	ofile.open(datafilenoext_ + "_residuals.txt");

	int j = 0;
	for (std::vector<std::array<double, 3> >::iterator dit = data_.begin(); dit != data_.end(); ++dit, ++j) {
		ofile << *(dit->data()) << "\t";
		ofile << *(dit->data() + 1) << "\t";
		ofile << *(dit->data() + 2) - displacements[j] << "\n";
	}

	ofile.close();

	ofile.open(datafilenoext_ + "_model_points.txt");
	std::vector<PointSource> sources = model.GetSources();
	for (std::vector<PointSource>::iterator it = sources.begin(); it != sources.end(); ++it) {
		ofile << it->GetX() << "\t";
		ofile << it->GetY() << "\t";
		ofile << it->GetZ() << "\t";
		ofile << it->GetStrikeAngle() << "\t";
		ofile << it->GetDipAngle() << "\t";
		ofile << it->GetRakeAngle() << "\t";
		ofile << it->GetSeismicMoment() << "\n";
	}
	ofile.close();
}

void Fitter::RecordBestModel() {
//	std::vector<double>::iterator min = std::min_element(errors_.begin(), errors_.end());
//	int index = std::distance(errors_.begin(), min);
//	bestmodel_ = models_[index];
	RecordModel(bestmodel_);
}

void Fitter::RecordFittedData() {
	std::ofstream dout(datafilenoext_ + "_fitted_data.txt");
	for (std::array<double, 3> datapoint : data_) {
		dout << datapoint[0] << "\t";
		dout << datapoint[1] << "\t";
		dout << datapoint[2] << "\n";
	}
}

System& Fitter::GetBestGenModel() {
	std::vector<double>::iterator min = std::min_element(errors_.begin(), errors_.end());
	int index = std::distance(errors_.begin(), min);
	return models_[index];
}

double& Fitter::GetBestGenError() {
	std::vector<double>::iterator min = std::min_element(errors_.begin(), errors_.end());
	int index = std::distance(errors_.begin(), min);
	return errors_[index];
}

std::tuple<System, double> Fitter::GetBestGenModelAndError() {
	std::vector<double>::iterator min = std::min_element(errors_.begin(), errors_.end());
	int index = std::distance(errors_.begin(), min);
	return std::make_tuple(models_[index], errors_[index]);
}

ModelAnalysis::ModelAnalysis(const std::string param_filename,
							 const std::string data_filename) {

	if (data_filename.empty()) {
		DataGenerator gen(param_filename, "interferogram.txt");
		gen.Interferogram();
		fitter_ = Fitter(param_filename, "interferogram.txt");
	}
	else {
		fitter_ = Fitter(param_filename, data_filename);
	}
}

void ModelAnalysis::Run(const int iterations) {

	for (int i = 0; i < iterations; ++i) {
		fitter_.GeneticAlgorithm(false);
		System& bestmodel = fitter_.GetBestOverallModel();
		double& bestmodelerror = fitter_.GetBestOverallError();
		std::vector<double> bestmodeldata = fitter_.CalcDisplacementField(bestmodel);
		bestmodels_.push_back(bestmodel);
		modelerrors_.push_back(bestmodelerror);
		modeldata_.push_back(bestmodeldata);
	}
}

void ModelAnalysis::RecordStats() {

	fitter_.RecordFittedData();

	std::ofstream mout(fitter_.datafilenoext_ + "_modelsources.txt", std::ofstream::app);

	for (System model : bestmodels_) {
		for (PointSource source : model) {
			mout << source.GetX() << "\t";
			mout << source.GetY() << "\t";
			mout << source.GetZ() << "\t";
			mout << source.GetStrikeAngle() << "\t";
			mout << source.GetDipAngle() << "\t";
//			mout << source.GetRakeAngle() << "\t";
			mout << source.GetSeismicMoment() << "\n";
		}
	}
	mout.close();

	std::ofstream eout(fitter_.datafilenoext_ + "_errors.txt", std::ofstream::app);

	for (double error : modelerrors_) {
		eout << error << "\n";
	}

	eout.close();

	std::ofstream dout(fitter_.datafilenoext_ + "_modeldata.txt", std::ofstream::app);

	for (std::vector<double> data : modeldata_) {
		for (double displacement : data) {
			dout << displacement << "\t";
		}
		dout << "\n";
	}

	dout.close();
}

void ModelAnalysis::SensAnalysis(System sys, const int steps) {
	std::vector<System> systems;
	std::vector<double> errors;

	std::vector<std::string> param_names = {"x", "y", "z", "strike",
											"dip", "seismic moment"};
	std::vector<double> ranges = {4.0, 4.0, 4.0, M_PI/2.0, M_PI/8.0, 5.0e9};
	std::vector<double> truevals;

	std::vector<std::vector<double> > allvals;

	for (PointSource source : sys) {
		truevals.push_back(source.GetX());
		truevals.push_back(source.GetY());
		truevals.push_back(source.GetZ());
		truevals.push_back(source.GetStrikeAngle());
		truevals.push_back(source.GetDipAngle());
		truevals.push_back(source.GetSeismicMoment());
	}

	for (int j = 0; j < 6; ++j) {
		double truevalue = truevals[j];
		double range = ranges[j];
		double stepsize = ranges[j]*2/steps;
		std::vector<double> vals;
		for (int i = 0; i < steps; ++i) {
			vals.push_back(truevalue - range + stepsize*i);
		}
		allvals.push_back(vals);
	}

	std::vector<double> xvals = allvals[0];
	std::vector<double> yvals = allvals[1];
	std::vector<double> zvals = allvals[2];
	std::vector<double> strikevals = allvals[3];
	std::vector<double> dipvals = allvals[4];
	std::vector<double> momvals = allvals[5];

	for (double x : xvals) {
		System temp;
		temp.PlaceSource(x, truevals[1], truevals[2], truevals[3],
						 truevals[4], 0, truevals[5]);
		systems.push_back(temp);
	}

	for (double y : yvals) {
		System temp;
		temp.PlaceSource(truevals[0], y, truevals[2], truevals[3],
						 truevals[4], 0, truevals[5]);
		systems.push_back(temp);
	}

	for (double z : zvals) {
		System temp;
		temp.PlaceSource(truevals[0], truevals[1], z, truevals[3],
						 truevals[4], 0, truevals[5]);
		systems.push_back(temp);
	}

	for (double strike : strikevals) {
		System temp;
		temp.PlaceSource(truevals[0], truevals[1], truevals[2],
						 strike, truevals[4], 0, truevals[5]);
		systems.push_back(temp);
	}

	for (double dip : dipvals) {
		System temp;
		temp.PlaceSource(truevals[0], truevals[1], truevals[2],
						 truevals[3], dip, 0, truevals[5]);
		systems.push_back(temp);
	}

	for (double mom : momvals) {
		System temp;
		temp.PlaceSource(truevals[0], truevals[1], truevals[2],
						 truevals[3], truevals[4], 0, mom);
		systems.push_back(temp);
	}

	std::vector<double> datadisp = fitter_.GetDataDisplacementField();
	for (System system : systems) {
		std::vector<double> modeldisp = fitter_.CalcDisplacementField(system);
		double error = fitter_.TestFitness(datadisp, modeldisp);
		errors.push_back(error);
	}

	std::ofstream sout("sensanalysis.txt");

	int i = 0;
	for (std::vector<double> paramvals : allvals) {
		for (double paramval : paramvals) {
			sout << paramval << "\t";
			sout << errors[i] << "\n";
			i++;
		}
	}
}

void ModelAnalysis::TwoDSensAnalysis(System sys, const int steps) {
	std::vector<System> systems;
	std::vector<double> errors;

	std::vector<std::string> param_names = {"x", "y", "z", "strike",
											"dip", "seismic moment"};
	std::vector<double> ranges = {4.0, 4.0, 4.0, M_PI/2.0, M_PI/8.0, 5.0e9};
	std::vector<double> truevals;

	std::vector<std::vector<double> > allvals;

	std::vector<std::vector<double> > param_sets;

	for (PointSource source : sys) {
		truevals.push_back(source.GetX());
		truevals.push_back(source.GetY());
		truevals.push_back(source.GetZ());
		truevals.push_back(source.GetStrikeAngle());
		truevals.push_back(source.GetDipAngle());
		truevals.push_back(source.GetSeismicMoment());
	}

	for (int j = 0; j < 6; ++j) {
		double truevalue = truevals[j];
		double range = ranges[j];
		double stepsize = ranges[j]*2/steps;
		std::vector<double> vals;
		for (int i = 0; i < steps; ++i) {
			vals.push_back(truevalue - range + stepsize*i);
		}
		allvals.push_back(vals);
	}

	std::vector<double> xvals = allvals[0];
	std::vector<double> yvals = allvals[1];
	std::vector<double> zvals = allvals[2];
	std::vector<double> strikevals = allvals[3];
	std::vector<double> dipvals = allvals[4];
	std::vector<double> momvals = allvals[5];

	for (int i = 0; i < steps*steps*36; ++i) {
		param_sets.push_back({truevals[0], truevals[1], truevals[2],
								  truevals[3], truevals[4], truevals[5]});
	}

	int k = 0;
	for (int j = 0; j < 6; ++j) {
		for (int i = 0; i < 6; ++i) {
//			if (i != j) {
				for (double valj : allvals[j]) {
					for (double vali : allvals[i]) {

						param_sets[k][i] = vali;
						param_sets[k][j] = valj;
						k++;

					}
				}
//			}
		}
	}

	for (std::vector<double> set: param_sets) {
		System temp;
		temp.PlaceSource(set[0], set[1], set[2],
						 set[3], set[4], 0, set[5]);
		systems.push_back(temp);
	}

	std::vector<double> datadisp = fitter_.GetDataDisplacementField();
	for (System system : systems) {
		std::vector<double> modeldisp = fitter_.CalcDisplacementField(system);
		double error = fitter_.TestFitness(datadisp, modeldisp);
		errors.push_back(error);
	}

	std::ofstream sout("twodsensanalysis.txt");

	for (double error : errors) {
		sout << error << "\n";
	}

//	for (std::vector set: param_sets) {
//		twout << set[0] << "\t";
//		twout << set[1] << "\t";
//		twout << set[2] << "\t";
//		twout << set[3] << "\t";
//		twout << set[4] << "\t";
//		twout << set[5] << "\n";
//	}
}
