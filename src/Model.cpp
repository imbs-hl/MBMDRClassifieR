
#include <math.h>
#include "Model.h"

Model::Model() : data(0), order(0), model_index(0), features(0), feature_names(0), n(0), alpha(0), in_cell(0), out_cell(0), cell_predictions(0), cell_statistics(0), cell_pvalues(0), cell_labels(0), statistic(0), pvalue(0), logger(0) {
}

Model::Model(Data* data,
		size_t order,
		size_t model_index,
		std::vector<size_t> features,
		double alpha,
		Logger* logger) : data(data), order(order), model_index(model_index), features(features), n(0), alpha(alpha), statistic(0), pvalue(0) {
	size_t idxs = pow(3, order);
	for(size_t i = 0; i < idxs; ++i) {
		this->in_cell.push_back(0);
		this->out_cell.push_back(0);
		this->cell_predictions.push_back(0.5);
		this->cell_statistics.push_back(0);
		this->cell_pvalues.push_back(0);
		this->cell_labels.push_back(0);
	}
	for(size_t o = 0; o < order; ++o) {
		std::string feature_name = data->getFeatureName(features[o]);
		this->feature_names.push_back(feature_name);
	}
	this->logger = logger;
}

Model::Model(Data* data,
		size_t order,
		size_t model_index,
		std::vector<size_t> features,
		std::vector<std::string> feature_names,
		double alpha,
		Logger* logger) : data(data), order(order), model_index(model_index), features(features), feature_names(feature_names), n(0), alpha(alpha), statistic(0), pvalue(0) {
	size_t idxs = pow(3, order);
	for(size_t i = 0; i < idxs; ++i) {
		this->in_cell.push_back(0);
		this->out_cell.push_back(0);
		this->cell_predictions.push_back(0.5);
		this->cell_statistics.push_back(0);
		this->cell_pvalues.push_back(0);
		this->cell_labels.push_back(0);
	}
	this->logger = logger;
}

Model::~Model() {
	// Empty on purpose
}

std::vector<uint> Model::getObservationsInCell() const {
	return this->in_cell;
}
std::vector<uint> Model::getObservationsOutCell() const {
	return this->out_cell;
}
std::vector<double> Model::getCellPredictions() const {
	return this->cell_predictions;
}
std::vector<double> Model::getCellStatistics() const {
	return this->cell_statistics;
}
std::vector<double> Model::getCellPValues() const {
	return this->cell_pvalues;
}
std::vector<int> Model::getCellLabels() const {
	return this->cell_labels;
}
double Model::getModelStatistic() const {
	return this->statistic;
}
double Model::getModelPValue() const {
	return this->pvalue;
}
size_t Model::getOrder() const {
	return this->order;
}
size_t Model::getModelIndex() const {
	return this->model_index;
}
std::vector<size_t> Model::getFeatures() const {
	return this->features;
}
std::vector<std::string> Model::getFeatureNames() const {
	return this->feature_names;
}
size_t Model::getNumObservations() const {
	return this->n;
}
double Model::getAlpha() const {
	return this->alpha;
}

void Model::loadModel(std::vector<uint> in_cell,
		std::vector<uint> out_cell,
		std::vector<double> cell_predictions,
		std::vector<double> cell_statistics,
		std::vector<double> cell_pvalues,
		std::vector<int> cell_labels,
		double statistic,
		double pvalue) {
	this->in_cell = in_cell;
	this->out_cell = out_cell;
	this->cell_predictions = cell_predictions;
	this->cell_statistics = cell_statistics;
	this->cell_pvalues = cell_pvalues;
	this->cell_labels = cell_labels;
	this->statistic = statistic;
	this->pvalue = pvalue;
}

std::vector<double> Model::predict() {
	// Get total number of observations
	size_t n_obs = data->getNumObservations();

	// Initialize output
	std::vector<double> predictions(n_obs);

	// Initialize base vector
	std::vector<uint> bases(order);
	for(size_t i = 0; i < order; ++i) {
		bases[i] = pow(3, i);
	}

	// Iterate through all samples in dataset
	logger->log(Info, "Iterating through samples...", 3);
	for(size_t i = 0; i < n_obs; ++i) {

		// Get genotype combination
		std::vector<uint> genotype(order);
		for(size_t j = 0; j < order; ++j) {
			genotype[j] = data->getFeature(i, feature_names[j]);
		}
		// Convert genotype combination to index
		int idx = std::inner_product(genotype.begin(), genotype.end(), bases.begin(), 0);

		if(cell_labels[idx] == 0) {
			predictions[i] = NA_REAL;
		} else {
			predictions[i] = cell_predictions[idx];
		}

	}

	return predictions;
}
