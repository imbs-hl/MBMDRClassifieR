
#include <thread>
#include "Mbmdr.h"
#include "ModelClassification.h"

Mbmdr::Mbmdr() :
data(0),
order(0),
n(0),
alpha(0),
max_models(0),
num_models(0),
feature_combination(0),
predictions(0),
v_levels(0) {
	this->mode = 1;
	this->models = std::priority_queue<Model*, std::vector<Model*>, CompareModelPointers>();

	this->threads = 1;
}

Mbmdr::Mbmdr(Data* data,
		size_t order,
		double alpha,
		size_t max_models,
		size_t mode,
		size_t num_threads,
		std::vector<std::ostream*> v_levels) {
	this->mode = mode;
	this->data = data;
	this->order = order;
	this->alpha = alpha;
	this->max_models = max_models;
	this->num_models = 0;
	this->n = data->getNumObservations();
	this->models = std::priority_queue<Model*, std::vector<Model*>, CompareModelPointers>();
	this->feature_combination = std::vector<size_t>(order);
	std::fill(this->feature_combination.begin(), this->feature_combination.end(), 0);
	this->predictions = std::vector<double>(0);

	if(num_threads == DEFAULT_NUM_THREADS) {
		// Detect cores
		this->threads = std::thread::hardware_concurrency();
	} else {
		this->threads = num_threads;
	}

	this->v_levels = v_levels;
}

Mbmdr::Mbmdr(Data* data, Rcpp::List saved_mbmdr,
		size_t num_threads,
		std::vector<std::ostream*> v_levels) {
	// Check if input is of correct type
	try {
		if(!saved_mbmdr.inherits("mbmdr")) {
			throw std::runtime_error("");
		}
	} catch (...) {
		throw std::runtime_error("Object must be a MB-MDR object.");
	}

	this->data = data;
	this->n = data->getNumObservations();
	this->predictions = std::vector<double>(n);
	this->models = std::priority_queue<Model*, std::vector<Model*>, CompareModelPointers>(CompareModelPointers(true));

	// Construct new models from saved MB-MDR object
	try {
		this->mode = saved_mbmdr["mode"];
		this->order = saved_mbmdr["order"];
		this->alpha = saved_mbmdr["alpha"];
		this->max_models = saved_mbmdr["max_models"];

		Rcpp::List models = saved_mbmdr["models"];

		for(int i = 0; i < models.size(); ++i) {
			Rcpp::List rcpp_model = models[i];

			double order = rcpp_model["order"];
			std::vector<size_t> features = rcpp_model["features"];
			double alpha = rcpp_model["alpha"];
			std::vector<uint> in_cell = rcpp_model["in_cell"];
			std::vector<uint> out_cell = rcpp_model["out_cell"];
			std::vector<double> cell_predictions = rcpp_model["cell_predictions"];
			std::vector<double> cell_statistics = rcpp_model["cell_statistics"];
			std::vector<double> cell_pvalues = rcpp_model["cell_pvalues"];
			std::vector<int> cell_labels = rcpp_model["cell_labels"];
			double statistic = rcpp_model["statistic"];
			double pvalue = rcpp_model["pvalue"];

			Model* model;

			if(mode == 1) {
				// Create classification model
				*v_levels[2] << "Creating classification model..." << std::endl;
				model = new ModelClassification(data,
						order,
						features,
						alpha,
						v_levels);
			} else if(mode == 2) {
				// Create regression model
				*v_levels[2] << "Creating regression model..." << std::endl;
			}

			model->loadModel(in_cell, out_cell, cell_predictions, cell_statistics, cell_pvalues, cell_labels, statistic, pvalue);

			possiblyAdd(model);

		}

		this->num_models = this->models.size();

	} catch(...) {
		throw std::runtime_error("Failed loading MB-MDR object.");
	}

	this->feature_combination = std::vector<size_t>(order);
	std::fill(this->feature_combination.begin(), this->feature_combination.end(), 0);

	if(num_threads == DEFAULT_NUM_THREADS) {
		// Detect cores
		this->threads = std::thread::hardware_concurrency();
	} else {
		this->threads = num_threads;
	}

	this->v_levels = v_levels;
}

Mbmdr::~Mbmdr() {
	while(!this->models.empty()) {
		Model* model = models.top();
		models.pop();
		delete model;
	}
}

void Mbmdr::possiblyAdd(Model* new_model) {

	*v_levels[2] << "Comparing model statistics..." << std::endl;

	// Get model statistic of new model
	double new_model_statistic = new_model->getModelStatistic();

	// Insert in queue if not full yet
	*v_levels[2] << "Statistic of new model: " << new_model_statistic << std::endl;
	if(models.size() < max_models) {
		*v_levels[2] << "Queue not full yet. Adding new model..." << std::endl;
		models.push(new_model);
		return;
	}

	// Get model statistic of model with lowest statistic in queue
	double top_model_statistic = models.top()->getModelStatistic();
	*v_levels[2] << "Statistic of model at top position: " << top_model_statistic << std::endl;
	if(new_model_statistic > top_model_statistic) {
		// get rid of the root, i.e. feature model with lowest statistic
		*v_levels[2] << "Adding new model..." << std::endl;
		Model* model = models.top();
		models.pop();
		delete model;

		// Add new model
		models.push(new_model);
	} else {
		// new model has lower statistic than all other models
		*v_levels[2] << "Discarding new model..." << std::endl;
		delete new_model;
	}
}

// [[Rcpp::plugins(cpp11)]]
bool Mbmdr::getNextFeatureCombination(size_t j) {

	size_t k = data->getNumFeatures();

	// Check if last feature combination is reached
	int last_feature_combination_sum = 0;
	int feature_combination_sum = 0;
	for(size_t i = 0; i < order; ++i) {
		last_feature_combination_sum += (k-i-1);
	}
	for(auto& n : feature_combination) {
		feature_combination_sum += n;
	}
	if (feature_combination_sum == last_feature_combination_sum) {
		return(false);
	}

	if (feature_combination[j] < (k-1) - ((order-1) - j)) {
		// Current index to increase is less than maximum at this index
		feature_combination[j] = feature_combination[j] + 1;
	} else {
		// Current index is at its maximum
		getNextFeatureCombination(j - 1);

		// Restart iteration in current index
		feature_combination[j] = feature_combination[j - 1] + 1;

	}

	return(true);

}

// [[Rcpp::plugins(cpp11)]]
void Mbmdr::fit() {

	// Create thread pool
	*v_levels[1] << "Creating thread pool..." << std::endl;
	std::vector<std::thread> thread_pool;
	thread_pool.reserve(threads);

	// Add threads to thread pool
	for(size_t t=0; t<threads; ++t) {
		*v_levels[2] << "Adding thread..." << std::endl;
		thread_pool.push_back(std::thread(&Mbmdr::fitModelInThread, this));
	}

	// Wait for completion
	*v_levels[1] << "Waiting for threads to complete..." << std::endl;
	for(auto &thread : thread_pool) {
		thread.join();
	}


}

void Mbmdr::fitModelInThread() {

	while(true) {

		// Feed thread with next feature combination
		std::unique_lock<std::mutex> lock(mutex);
		if(getNextFeatureCombination(order-1)) {
			std::vector<size_t> feature_combination = this->feature_combination;

			*v_levels[2] << "Calculating feature combination ";
			for(auto &col : feature_combination) {
				*v_levels[2] << col << " ";
			}
			*v_levels[2] << "..." << std::endl;

			// Increase model counter
			++num_models;
		} else {
			*v_levels[2] << "No further feature combinations available!" << std::endl;
			break;
		}
		lock.unlock();

		Model* model;

		if(mode == 1) {
			// Create classification model
			*v_levels[2] << "Creating classification model..." << std::endl;
			model = new ModelClassification(data,
					order,
					feature_combination,
					alpha,
					v_levels);
		} else if(mode == 2) {
			// Create regression model
			*v_levels[2] << "Creating regression model..." << std::endl;
		}

		// Fit the model
		*v_levels[2] << "Model fit in progress..." << std::endl;
		model->fit();

		// (Possibly) save model
		lock.lock();
		*v_levels[2] << "Saving model..." << std::endl;
		possiblyAdd(model);
		lock.unlock();
	}

}

Rcpp::List Mbmdr::exportModels() {

	Rcpp::List export_object;

	while(!models.empty()) {
		// Get next model
		Model* model = models.top();

		// Initialize model export object
		Rcpp::List export_model_object;

		// Fill model export object
		export_model_object.push_back(model->getFeatures(), "features");
		export_model_object.push_back(model->getObservationsInCell(), "in_cell");
		export_model_object.push_back(model->getObservationsOutCell(), "out_cell");
		export_model_object.push_back(model->getCellPredictions(), "cell_predictions");
		export_model_object.push_back(model->getCellStatistics(), "cell_statistics");
		export_model_object.push_back(model->getCellPValues(), "cell_pvalues");
		export_model_object.push_back(model->getCellLabels(), "cell_labels");
		export_model_object.push_back(model->getModelStatistic(), "statistic");
		export_model_object.push_back(model->getModelPValue(), "pvalue");
		export_model_object.push_back(model->getOrder(), "order");
		export_model_object.push_back(model->getAlpha(), "alpha");

		// Append model export object to list of models
		export_object.push_front(export_model_object);

		// Release model from queue
		models.pop();

		// Delete model
		delete model;
	}

	return export_object;

}


// [[Rcpp::plugins(cpp11)]]
std::vector<double> Mbmdr::predict() {

	// Create thread pool
	*v_levels[1] << "Creating thread pool..." << std::endl;
	std::vector<std::thread> thread_pool;
	thread_pool.reserve(threads);

	// Add threads to thread pool
	for(size_t t=0; t<threads; ++t) {
		*v_levels[2] << "Adding thread..." << std::endl;
		thread_pool.push_back(std::thread(&Mbmdr::predictInThread, this));
	}

	// Wait for completion
	*v_levels[1] << "Waiting for threads to complete..." << std::endl;
	for(auto &thread : thread_pool) {
		thread.join();
	}

	for(uint i = 0; i < predictions.size(); ++i) {
		predictions[i] /= num_models;
	}

	return predictions;

}

void Mbmdr::predictInThread() {

	while(!models.empty()) {

		// Feed thread with next model
		std::unique_lock<std::mutex> lock(mutex);
		Model* model = models.top();
		models.pop();
		lock.unlock();

		// Predict sample outcomes in model
		std::vector<double> model_predictions = model->predict();

		// Add predictions
		lock.lock();
		*v_levels[2] << "Saving predictions..." << std::endl;
		for(uint i = 0; i < predictions.size(); ++i) {
			predictions[i] += model_predictions[i];
		}
		lock.unlock();
		delete model;
	}

}

// Getters
uint Mbmdr::getMode() {
	return this->mode;
}
size_t Mbmdr::getOrder() {
	return this->order;
}
size_t Mbmdr::getN() {
	return this->n;
}
double Mbmdr::getAlpha() {
	return this->alpha;
}
size_t Mbmdr::getMaxModels() {
	return this->max_models;
}
std::priority_queue<Model*, std::vector<Model*>, CompareModelPointers> Mbmdr::getModels() {
	return this->models;
}
