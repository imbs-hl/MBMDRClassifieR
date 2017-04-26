
#include "Model.h"

Model::Model() : data(0), order(0), features(0), n(0), alpha(0), in_cell(0), out_cell(0), cell_statistics(0), cell_pvalues(0), cell_labels(0), statistic(0), pvalue(0), v_levels(0) {
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
double Model::getAlpha() const{
	return this->alpha;
}
