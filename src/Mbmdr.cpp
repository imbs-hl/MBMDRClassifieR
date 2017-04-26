
#include <thread>
#include "Mbmdr.h"
#include "ModelClassification.h"

Mbmdr::Mbmdr() :
data(0),
order(0),
n(0),
alpha(0),
max_models(0),
feature_combination(0),
v_levels(0) {
	this->mode = 1;
	this->models = std::priority_queue<Model*, std::vector<Model*>, CmpModelPtrs>();

	this->threads = 1;
}

Mbmdr::Mbmdr(Data* data,
		size_t order,
		double alpha,
		size_t max_models,
		size_t mode,
		size_t num_threads,
		std::vector<std::ostream*> v_levels,
		Rcpp::List saved_mbmdr) {
	this->mode = mode;
	this->data = data;
	this->order = order;
	this->alpha = alpha;
	this->max_models = max_models;
	this->n = data->getNumObservations();
	this->models = std::priority_queue<Model*, std::vector<Model*>, CmpModelPtrs>();
	this->feature_combination = std::vector<size_t>(order);
	std::fill(this->feature_combination.begin(), this->feature_combination.end(), 0);

	if(num_threads == DEFAULT_NUM_THREADS) {
		// Detect cores
		this->threads = std::thread::hardware_concurrency();
	} else {
		this->threads = num_threads;
	}

	this->v_levels = v_levels;

	// Load saved MB-MDR object
	if(saved_mbmdr.size() > 0) {

	}
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

	// Add threds to thread pool
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
