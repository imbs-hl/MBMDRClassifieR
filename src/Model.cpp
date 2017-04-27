
#include "Model.h"
#include <math.h>

Model::Model() : data(0), order(0), features(0), n(0), alpha(0), in_cell(0), out_cell(0), cell_predictions(0), cell_statistics(0), cell_pvalues(0), cell_labels(0), statistic(0), pvalue(0), v_levels(0) {
}

Model::Model(Data* data,
		size_t order,
		std::vector<size_t> features,
		double alpha,
		std::vector<std::ostream*> v_levels) : data(data), order(order), features(features), n(0), alpha(alpha), statistic(0), pvalue(0) {
	size_t idxs = pow(3, order);
	for(size_t i = 0; i < idxs; ++i) {
		this->in_cell.push_back(0);
		this->out_cell.push_back(0);
		this->cell_predictions.push_back(0);
		this->cell_statistics.push_back(0);
		this->cell_pvalues.push_back(0);
		this->cell_labels.push_back(0);
	}
	this->v_levels = v_levels;
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
std::vector<size_t> Model::getFeatures() const {
	return this->features;
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
	*v_levels[2] << "Iterating through samples..." << std::endl;
	for(size_t i = 0; i < n_obs; ++i) {

		// Get genotype combination
		std::vector<uint> genotype(order);
		for(size_t j = 0; j < order; ++j) {
			genotype[j] = data->getFeature(i, features[j]);
		}
		// Convert genotype combination to index
		int idx = std::inner_product(genotype.begin(), genotype.end(), bases.begin(), 0);

		predictions[i] = cell_predictions[idx];

	}

	return predictions;
}
