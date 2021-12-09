#include <iomanip>
#include <iostream>
#include <chrono>
#include "insarinverse.h"

int CountColumns(const std::string& filename) {
    std::ifstream fin(filename);
    std::string line;
    std::getline(fin, line);
    std::istringstream iss(line);
    std::string value;

    int count = 0;
    while (iss >> value)
        count++;
    return count;
}

Position::Position(const double x, const double y, const double z) :
		x(x), y(y), z(z) {

}

Data::Data(Position& pos, double dx, double dy, double dz) :
		pos(pos), dx(dx), dy(dy), dz(dz) {

}

Displacement::Displacement(double dx, double dy, double dz) :
		dx(dx), dy(dy), dz(dz) {

}

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
		return "0";
	}
}

Parameters::Parameters(const std::string& param_filename) {

	ParameterReader params(param_filename);
//	found = params.GetFoundParameters();

	xmin = std::stod(params.GetParameter("Minimum data x value"));
	xmax = std::stod(params.GetParameter("Maximum data x value"));
	ymin = std::stod(params.GetParameter("Minimum data y value"));
	ymax = std::stod(params.GetParameter("Maximum data y value"));
	min_moment = std::stod(params.GetParameter("Minimum seismic moment"));
	max_moment = std::stod(params.GetParameter("Maximum seismic moment"));
	chance_to_mutate = std::stod(params.GetParameter("Chance to mutate"));
	strike_angle = std::stod(params.GetParameter("Strike angle"));
	dip_angle = std::stod(params.GetParameter("Dip angle"));
	rake_angle = std::stod(params.GetParameter("Rake angle"));

	xdatapoints = std::stoi(params.GetParameter("Num x data points"));
	ydatapoints = std::stoi(params.GetParameter("Num y data points"));
	gen_num_sources = std::stoi(params.GetParameter("Num sources (gen)"));
	fit_num_sources = std::stoi(params.GetParameter("Num sources (fit)"));
	num_generations = std::stoi(params.GetParameter("Number of generations"));
	pop_size = std::stoi(params.GetParameter("Population size"));
	xfitpoints = std::stoi(params.GetParameter("Sources in x direction"));
	yfitpoints = std::stoi(params.GetParameter("Sources in y direction"));
	zfitpoints = std::stoi(params.GetParameter("Sources in z direction"));

	sxmin = std::stod(params.GetParameter("Minimum source x value"));
	sxmax = std::stod(params.GetParameter("Maximum source x value"));
	symin = std::stod(params.GetParameter("Minimum source y value"));
	symax = std::stod(params.GetParameter("Maximum source y value"));
	szmin = std::stod(params.GetParameter("Minimum source z value"));
	szmax = std::stod(params.GetParameter("Maximum source z value"));

//	source_type = params.GetParameter("Source type");

}

bool Parameters::DataLimits() {
	if (xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0) {
		std::cout << "No data limits set" << std::endl;
		return false;
	}
	else {
		return true;
	}
}

bool Parameters::SourceLimits() {
	if (sxmin == 0 && sxmax == 0 && symin == 0 && symax == 0) {
		std::cout << "No source limits set" << std::endl;
		return false;
	}
	else {
		return true;
	}
}

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

void InSAR::BinnerNepalQuake(std::vector<std::string> files) {

	const int latbins = 75; //TODO: Fix so that bins can differ in each direction
	const int lonbins = 75;

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
	}

	for (int j = 0; j < lonbins + 1; j++) {
		lonbounds[j] = lonmin + (lonmax - lonmin)/lonbins*j;
	}

	for (std::vector<std::string>::iterator it = files.begin(); it != files.end(); ++it) {

		std::ifstream datain(*it);
		if (!datain) {std::cout << "nope" << "\n";}

		while (datain >> longitude >> latitude >> topo >> east >> north >> up >> los >> sig) {
			for (int i = 0; i < latbins; i++) {
				for (int j = 0; j < lonbins; j++) {
					if (latitude <= latbounds[i + 1] && longitude <= lonbounds[j + 1]) {
						hist[j][i] += los/up;
//						hist[j][i] += los;
						counts[j][i]++;
						goto finish;
					}
				}
			}
			finish:;
		}
	}

	std::ofstream smout("mainshock_displacement_binned_100x100.txt");
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

void InSAR::MultiBinnerNepalQuake(std::vector<std::string> files) {

	const int latbins = 50; //TODO: Fix so that bins can differ in each direction
	const int lonbins = 50;

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

	for (int i = 0; i < latbins + 1; i++) {
		latbounds[i] = latmin + (latmax - latmin)/latbins*i;
	}

	for (int j = 0; j < lonbins + 1; j++) {
		lonbounds[j] = lonmin + (lonmax - lonmin)/lonbins*j;
	}

	for (std::vector<std::string>::iterator it = files.begin(); it != files.end(); ++it) {

		std::array<std::array<double, latbins>, lonbins> hist {0};
		std::array<std::array<int, latbins>, lonbins> counts {0};

		std::array<std::array<double, latbins>, lonbins> easts {0};
		std::array<std::array<double, latbins>, lonbins> norths {0};
		std::array<std::array<double, latbins>, lonbins> ups {0};

		const std::string filename = *it;
		size_t lastindex = filename.find_last_of(".");
		std::string filenoext = filename.substr(0, lastindex);

		std::ifstream datain(*it);
		if (!datain) {std::cout << "nope" << "\n";}

		while (datain >> longitude >> latitude >> topo >> east >> north >> up >> los >> sig) {
			for (int i = 0; i < latbins; i++) {
				for (int j = 0; j < lonbins; j++) {
					if (latitude <= latbounds[i + 1] && longitude <= lonbounds[j + 1]) {
						hist[j][i] += los;
						easts[j][i] += east;
						norths[j][i] += north;
						ups[j][i] += up;
						counts[j][i]++;
						goto finish;
					}
				}
			}
			finish:;
		}

		std::ofstream smout(filenoext + "_binned.txt");
		smout << std::fixed << std::setprecision(8);
		for (int i = 0; i < latbins; i++) {
			for (int j = 0; j < lonbins; j++) {
				smout << (lonbounds[j] + lonbounds[j + 1])/2.0 << "\t";
				smout << (latbounds[i] + latbounds[i + 1])/2.0 << "\t";
				smout << hist[j][i]/counts[j][i] << "\t";
				smout << easts[j][i]/counts[j][i] << "\t";
				smout << norths[j][i]/counts[j][i] << "\t";
				smout << ups[j][i]/counts[j][i] << "\n";
			}
		}
	}
}

void InSAR::TwoDBinner(const std::string& filename) {

	const int xbins = 50;
	const int ybins = 50;

	std::array<double, xbins + 1> xbounds;
	std::array<double, ybins + 1> ybounds;

	double x, y, east, up;

	double xmin = DBL_MAX;
	double xmax = -DBL_MAX;
	double ymin = DBL_MAX;
	double ymax = -DBL_MAX;

	std::ifstream dataread(filename);
	if (!dataread) {std::cout << "nope" << "\n";}

	while (dataread >> x >> y >> east >> up) {
		if (x > xmax) {xmax = x;}
		if (x < xmin) {xmin = x;}
		if (y > ymax) {ymax = y;}
		if (y < ymin) {ymin = y;}
	}

	for (int i = 0; i < xbins + 1; i++) {
		xbounds[i] = xmin + (xmax - xmin)/xbins*i;
	}

	for (int j = 0; j < ybins + 1; j++) {
		ybounds[j] = ymin + (ymax - ymin)/ybins*j;
	}

	std::array<std::array<int, ybins>, xbins> counts {0};

	std::array<std::array<double, ybins>, xbins> easts {0};
	std::array<std::array<double, ybins>, xbins> ups {0};

	size_t lastindex = filename.find_last_of(".");
	std::string filenoext = filename.substr(0, lastindex);

	std::ifstream datain(filename);
	if (!datain) {std::cout << "nope" << "\n";}

	while (datain >> x >> y >> east >> up) {
		for (int i = 0; i < ybins; i++) {
			for (int j = 0; j < xbins; j++) {
				if (y <= ybounds[i + 1] && x <= xbounds[j + 1]) {
					easts[j][i] += east;
					ups[j][i] += up;
					counts[j][i]++;
					goto finish;
				}
			}
		}
		finish:;
	}

	std::ofstream smout(filenoext + "_binned.txt");
	smout << std::fixed << std::setprecision(8);
	for (int i = 0; i < ybins; i++) {
		for (int j = 0; j < xbins; j++) {
			smout << (xbounds[j] + xbounds[j + 1])/2.0 << "\t";
			smout << (ybounds[i] + ybounds[i + 1])/2.0 << "\t";
			smout << easts[j][i]/counts[j][i] << "\t";
			smout << ups[j][i]/counts[j][i] << "\n";
		}
	}

}

PointSource::PointSource(const double x, const double y, const double z,
						 const double strike, const double dip,
						 const double rake, const double moment) :
		strike(strike),
		dip(dip),
		rake(rake),
		moment(moment) {
	pos.x = x;
	pos.y = y;
	pos.z = z;
}

PointSource PointSource::CopySource() {
	return PointSource(this->pos.x, this->pos.y, this->pos.z, this->strike, this->dip, this->rake, this->moment);
}

System::System(const std::string& filename) {

	int numcols = CountColumns(filename);

	double x, y, z, strike, dip, rake, moment;

	std::ifstream fin(filename);
	if (!fin) {std::cout << "Couldn't find file " << filename << "\n";}

	switch(numcols) {
		case 2:
			{
				while (fin >> x >> y) {
					sources_.emplace_back(x, y, 0, 0, 0, 1);
				}
				break;
			}
		case 7:
			{
				while (fin >> x >> y >> z >> strike >> dip >> rake >> moment) {
					sources_.emplace_back(x, y, z, strike, dip, rake, moment);
				}
				break;
			}
		default:
			{
				throw std::runtime_error("I don't know what to do with this file");
				break;
			}
	}
}

System::System(Parameters params) {

	double dz;
	double dx = (params.sxmax - params.sxmin)/(params.xfitpoints - 1);
	double dy = (params.symax - params.symin)/(params.yfitpoints - 1);

	if (params.zfitpoints == 1) {
		dz = 0;
	}
	else {
		dz = (params.szmax - params.szmin)/(params.zfitpoints - 1);
	}

	for (int k = 0; k < params.zfitpoints; ++k) {
		for (int j = 0; j < params.yfitpoints; ++j) {
			for (int i = 0; i < params.xfitpoints; ++i) {
				sources_.emplace_back(params.sxmin + i*dx, params.symin + j*dy, -(params.szmin + k*dz),
									  params.strike_angle, params.dip_angle, params.rake_angle);
			}
		}
	}
}

void System::CopySources(System sys) { //copy the sources from sys into the calling object
	this->ClearSources();
	for (std::vector<PointSource>::iterator it = sys.sources_.begin(); it != sys.sources_.end(); ++it) {
		this->sources_.push_back(it->CopySource());
	}
}

void System::Displacement(const std::string& filename, int numx, int numy,
				  double xmin, double xmax, double ymin, double ymax) {

	double dx = (xmax - xmin)/numx;
	double dy = (ymax - ymin)/numy;

	std::vector<double> xs, ys, disp;

	for (int i = 0; i <= numx; ++i) {
		xs.push_back(xmin + dx*i);
	}
	for (int j = 0; j <= numy; ++j) {
		ys.push_back(ymin + dy*j);
	}

	for (double y : ys) {
		for (double x : xs) {
			double temp = 0;
			for (PointSource src : sources_) {
				temp += z_disp(x, y, src.pos.x, src.pos.y, -src.pos.z, src.strike,
							   src.dip, src.rake, src.moment);
			}
			disp.push_back(temp);
		}
	}

	std::ofstream fout(filename);

	for (double dz : disp) {
		fout << dz << "\n";
	}
}

std::mt19937 DataGenerator::engine(time(0));

DataGenerator::DataGenerator(const std::string& param_filename, std::string const int_filename) {

	filename_ = int_filename;

	size_t lastindex = int_filename.find_last_of(".");
	filename_noext_ = int_filename.substr(0, lastindex);

	params_ = Parameters(param_filename);

}

std::vector<Data> DataGenerator::Interferogram() {

	double x, y, datadx, datady;
	std::vector<Data> interferogram;
	std::vector<PointSource> sources;

	double xmin = params_.xmin;
	double xmax = params_.xmax;
	double ymin = params_.ymin;
	double ymax = params_.ymax;

	double moment_order_min = std::log10(params_.min_moment);
	double moment_order_max = std::log10(params_.max_moment);

	int xdatapoints = params_.xdatapoints;
	int ydatapoints = params_.ydatapoints;
	int numsources = params_.gen_num_sources;

	std::uniform_real_distribution<double> dist_x(params_.sxmin, params_.sxmax);
	std::uniform_real_distribution<double> dist_y(params_.symin, params_.symax);
	std::uniform_real_distribution<double> dist_z(params_.szmin, params_.szmax);
	std::uniform_real_distribution<double> strike(0.0, 2.0*M_PI); //also used for rake
	std::uniform_real_distribution<double> dip(0.0, M_PI/2.0);
	std::uniform_real_distribution<double> order(moment_order_min,
												 moment_order_max);

	datadx = (xmax - xmin)/xdatapoints;
	datady = (ymax - ymin)/ydatapoints;

	for (int i = 0; i < numsources; ++i) {
		sources.push_back(PointSource(dist_x(engine), dist_y(engine), -dist_z(engine), strike(engine), dip(engine), strike(engine), pow(10, order(engine))));
	}

	for (y = ymin + datady/2; y < ymax; y += datady) {
		for (x = xmin + datadx/2; x < xmax; x += datadx) {
			Position pos(x, y, 0.0);
			interferogram.emplace_back(pos, 0.0, 0.0, 0.0);
		}
	}

	for (Data& data : interferogram) {

		double datax = data.pos.x;
		double datay = data.pos.y;

		for (PointSource& source : sources) {
			double sx = source.pos.x;
			double sy = source.pos.y;
			double sz = source.pos.z;
			double strike = source.strike;
			double dip = source.dip;
			double rake = source.rake;
			double moment = source.moment;

			double xstr = strike_disp_x(datax, datay, 0.0, sx, sy, -sz, strike, dip, moment);
			double xdip = dip_disp_x(datax, datay, 0.0, sx, sy, -sz, strike, dip, moment);

			double zstr = strike_disp_z(datax, datay, 0.0, sx, sy, -sz, strike, dip, moment);
			double zdip = dip_disp_z(datax, datay, 0.0, sx, sy, -sz, strike, dip, moment);

			data.dx += xstr*cos_(rake) + xdip*sin_(rake);
			data.dz += zstr*cos_(rake) + zdip*sin_(rake);
		}
	}

	std::ofstream pout(filename_noext_ + "_sources.txt");
	if (!pout) {std::cout << "Couldn't open file" << "\n";}
	for (PointSource& source : sources) {
		pout << source.pos.x << "\t";
		pout << source.pos.y << "\t";
		pout << source.pos.z << "\t";
		pout << source.strike << "\t";
		pout << source.dip << "\t";
		pout << source.rake << "\t";
		pout << source.moment << "\n";
	}
	pout.close();

	std::ofstream iout(filename_);
	if (!iout) {std::cout << "Couldn't open file" << "\n";}
	for (Data& data : interferogram) {
		iout << data.pos.x << "\t";
		iout << data.pos.y << "\t";
//		iout << data.pos.z << "\t";
		iout << data.dx << "\t";
//		iout << data.dy << "\t";
		iout << data.dz << "\n";
	}
	iout.close();

	return interferogram;
}

std::mt19937 Fitter::engine(time(0));

DataReader::DataReader(const std::string& data_filename) {

	size_t lastindex = data_filename.find_last_of(".");
	datafilenoext_ = data_filename.substr(0, lastindex);

	std::ifstream fin(data_filename);
	if (!fin) {std::cout << "nope" << "\n";}
	double x, y, dx, dz;
	std::string xin, zin;

	int ncols = CountColumns(data_filename);

	switch(ncols) {
		case 3:
			while (fin >> x >> y >> zin) {
				dz = std::stod(zin);

				if (!std::isnan(dz) && dz != 0.0) {
					Position pos(x, y, 0.0);
					data_.emplace_back(pos, 0.0, 0.0, dz);
				}
			}
			break;
		case 4:
			while (fin >> x >> y >> xin >> zin) {
				dx = std::stod(xin);
				dz = std::stod(zin);

				if (!std::isnan(dx) && !std::isnan(dz) && dx != 0.0 && dz != 0.0) {
					Position pos(x, y, 0.0);
					data_.emplace_back(pos, dx, 0.0, dz);
				}
			}
//			std::cout << data_.size() << " data points." << std::endl;
			break;
	}
}

std::vector<Data> DataReader::GetData() {
	return data_;
}

std::vector<Data> DataReader::Subset(const double xmin,
									 const double xmax,
									 const double ymin,
									 const double ymax) {
	double x, y;
	std::vector<Data> subset;

	for (Data& data : data_) {
		x = data.pos.x;
		y = data.pos.y;

		if (x <= xmax && x >= xmin && y <= ymax && y >= ymin) {
			subset.push_back(data);
	    }
	}

	return subset;
}

std::vector<double> DataReader::FindBounds() {

	double xmin = DBL_MAX;
	double xmax = -DBL_MAX;
	double ymin = DBL_MAX;
	double ymax = -DBL_MAX;

	for (Data& data : data_) {
		if (data.pos.x > xmax) {
			xmax = data.pos.x;
		}
		if (data.pos.x < xmin) {
			xmin = data.pos.x;
		}
		if (data.pos.y < ymin) {
			ymin = data.pos.y;
		}
		if (data.pos.y > ymax) {
			ymax = data.pos.y;
		}
	}

	return {xmin, xmax, ymin, ymax};
}

int DataReader::CountColumns(const std::string& filename) {
    std::ifstream fin(filename);
    std::string line;
    std::getline(fin, line);
    std::istringstream iss(line);
    std::string value;

    int count = 0;
    while (iss >> value)
        count++;
    return count;
}

GreensCatalog::GreensCatalog(System sys, std::vector<Data> data) {

	numsources = sys.size();
	numdatapoints = data.size();

	for (Data& datap : data) {

		double datax = datap.pos.x;
		double datay = datap.pos.y;

		for (PointSource& src : sys) {
//			StrikeGreens strgr;
//			DipGreens dipgr;

			Greens temp;

			double sx = src.pos.x;
			double sy = src.pos.y;
			double sz = -src.pos.z;
			double strike = src.strike;
			double dip = src.dip;
			double rake = src.rake;

			double strgrx = strike_disp_x(datax, datay, 0.0, sx, sy, sz, strike, dip, 1.0);
			double strgrz = strike_disp_z(datax, datay, 0.0, sx, sy, sz, strike, dip, 1.0);

			double dipgrx = dip_disp_x(datax, datay, 0.0, sx, sy, sz, strike, dip, 1.0);
			double dipgrz = dip_disp_z(datax, datay, 0.0, sx, sy, sz, strike, dip, 1.0);

			temp.u_x = strgrx*cos_(rake) + dipgrx*sin_(rake);
			temp.u_z = strgrz*cos_(rake) + dipgrz*sin_(rake);

			greens.push_back(temp);
		}
	}
}

void GreensCatalog::SaveCatalog(const std::string& filename) {

	size_t lastindex = filename.find_last_of(".");
	std::string filenoext = filename.substr(0, lastindex);

	std::vector<Greens>::iterator it;

	std::ofstream xout(filenoext + "_x.txt");
	std::ofstream zout(filenoext + "_z.txt");

	int i = 1;
	for (it = greens.begin(); it != greens.end(); ++it, ++i) {
		xout << it->u_x;
		zout << it->u_z;

		if (i % numsources == 0) {
			xout << "\n";
			zout << "\n";
		}
		else {
			xout << "\t";
			zout << "\t";
		}
	}
}

Fitter::Fitter(const std::string& param_filename, const std::string& data_filename) {

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

std::vector<double> Fitter::GetDataDisplacementField() {
	std::vector<double> displacements;
	for (Data& data : data_) {
		displacements.push_back(data.dz);
	}
	return displacements;
}

std::vector<double> Fitter::CalcDisplacementField(System model) {
	std::vector<double> displacements(data_.size(), 0.0);
	std::vector<PointSource> sources = model.GetSources();

	int i = 0;
	for (Data& data : data_) {
		for (std::vector<PointSource>::iterator it = sources.begin(); it != sources.end(); ++it) {
			double strz = strike_disp_z(data.pos.x, data.pos.y, 0.0, it->pos.x, it->pos.y, fabs(it->pos.z), it->strike, it->dip, it->moment);
			double dipz = dip_disp_z(data.pos.x, data.pos.y, 0.0, it->pos.x, it->pos.y, fabs(it->pos.z), it->strike, it->dip, it->moment);
			displacements[i] += strz*cos_(it->rake) + dipz*sin_(it->rake);
		}
		i++;
	}
	return displacements;
}

std::vector<double> Fitter::CalcDisplacementField(System model, GreensCatalog catalog) {
	std::vector<double> displacements(data_.size(), 0.0);
	std::vector<PointSource> sources = model.GetSources();

	int i = 0;

	for (double& disp : displacements) {
		for (std::vector<PointSource>::iterator it = sources.begin(); it != sources.end(); ++it, ++i) {
			disp += catalog.greens[i].u_z*it->moment;
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
		p1x = pit->pos.x;
		p2x = sit->pos.x;
		p1y = pit->pos.y;
		p2y = sit->pos.y;
		p1z = pit->pos.z;
		p2z = pit->pos.z;

		p1strike = pit->strike;
		p2strike = sit->strike;
		p1dip = pit->dip;
		p2dip = sit->dip;
		p1mom = pit->moment;
		p2mom = sit->moment;

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

		pit->pos.x = p1x_new;
		sit->pos.x = p2x_new;
		pit->pos.y = p1y_new;
		sit->pos.y = p2y_new;
		pit->pos.z = p1z_new;
		sit->pos.z = p2z_new;

		pit->strike = p1strike_new;
		sit->strike = p2strike_new;
		pit->dip = p1dip_new;
		sit->dip = p2dip_new;
		pit->moment = p1mom_new;
		sit->moment = p2mom_new;
	}

	new_sys.push_back(parent1);
	new_sys.push_back(parent2);

	return new_sys;

}

std::vector<System> Fitter::SBXMoment(System parent1, System parent2) {

	double beta, eta_c, u;
	double p1mom, p2mom;
	double p1mom_new, p2mom_new;

	std::vector<PointSource> parent1_sources = parent1.GetSources();
	std::vector<PointSource> parent2_sources = parent2.GetSources();
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

void Fitter::LocationMutation (System &sys) {
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

void Fitter::OrientationMutation (System &sys) {
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

void Fitter::MomentMutation (System &sys) {
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

	std::uniform_int_distribution<int> dist_x(params_.sxmin, params_.sxmax);
	std::uniform_int_distribution<int> dist_y(params_.symin, params_.symax);
	std::uniform_int_distribution<int> dist_z(params_.szmin, params_.szmax);
	std::uniform_real_distribution<double> strike(0.0, 2.0*M_PI);
	std::uniform_real_distribution<double> dip(0.0, M_PI/2.0);
	std::uniform_real_distribution<double> moment(params_.min_moment,
												  params_.max_moment);
	std::uniform_real_distribution<double> mutate(0.0, 1.0);

	System model;
	System bestmodel;

	std::vector<System> next_gen;

	double error;

	std::vector<double> datadisp = this->GetDataDisplacementField();

	for (int i = 0; i < params_.pop_size; i++) { //generate random first generation parents

		model.ClearSources();

		for (int j = 0; j < params_.fit_num_sources; j++) { //number of sources in each parent
			model.PlaceSource(dist_x(engine), dist_y(engine), -dist_z(engine),
							  strike(engine), dip(engine), 0, moment(engine));
		}

		std::vector<double> modeldisp = CalcDisplacementField(model);

		error = TestFitness(datadisp, modeldisp);

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
	std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_time = end - start;
	std::cout << elapsed_time.count() << "s\n";
}

void Fitter::MomentGA(bool save, std::string sourcefile) {

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

	std::vector<double> datadisp = this->GetDataDisplacementField();

	for (int i = 0; i < params_.pop_size; i++) { //generate random first generation parents

		for (PointSource& source : model) {
			source.moment = moment(engine);
		}

//		std::vector<double> modeldisp = CalcDisplacementField(model);
		std::vector<double> modeldisp = CalcDisplacementField(model, catalog);

		error = TestFitness(datadisp, modeldisp);

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
				std::vector<System> crossed = SBXMoment(parent1, parent2);
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
//				std::vector<double> modeldisp = CalcDisplacementField(next_gen[k]);
				std::vector<double> modeldisp = CalcDisplacementField(next_gen[k], catalog);
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

double Fitter::TestFitness(std::vector<double> datadisp, std::vector<double> modeldisp) {
	double error;
	double temp = 0;
	double n = 2.0;

	if (datadisp.size() != modeldisp.size()) {throw std::runtime_error("Vectors are different lengths");}

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
	for (std::vector<Data>::iterator dit = data_.begin();
										 dit != data_.end(); ++dit, ++i) {
		ofile << dit->pos.x << "\t";
		ofile << dit->pos.y << "\t";
		ofile << displacements[i] << "\n";
	}

	ofile.close();

//	ofile.open(datafilenoext_ + "_residuals.txt");
//
//	int j = 0;
//	for (std::vector<Data>::iterator dit = data_.begin(); dit != data_.end(); ++dit, ++j) {
//		ofile << dit->pos.x << "\t";
//		ofile << dit->pos.y << "\t";
//		ofile << dit->dz - displacements[j] << "\n";
//	}
//
//	ofile.close();

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

void Fitter::RecordBestModel() {
	RecordModel(bestmodel_);
}

void Fitter::RecordFittedData() {
	std::ofstream dout(datafilenoext_ + "_fitted_data.txt");

	for (Data& data : data_) {
		dout << data.pos.x << "\t";
		dout << data.pos.y << "\t";
//		dout << data.pos.z << "\t";
//		dout << data.dx << "\t";
//		dout << data.dy << "\t";
		dout << data.dz << "\n";
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

void Fitter::SaveGreensCatalog(const std::string& filename) {
	System model = System(params_);

	GreensCatalog catalog(model, data_);
	catalog.SaveCatalog(filename);
}

void Fitter::Setup() {

	System model = System(params_);
	RecordModel(model);
	RecordFittedData();

	GreensCatalog catalog(model, data_);
	catalog.SaveCatalog("greens.txt");

}
