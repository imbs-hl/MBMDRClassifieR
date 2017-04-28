
#ifndef MBMDR_H_
#define MBMDR_H_

#include <Rcpp.h>
#include <queue>
#include <mutex>
#include <functional>
#include "globals.h"
#include "Model.h"
#include "Data.h"

class Mbmdr {

public:

	// Constructors
	Mbmdr();
	Mbmdr(Data* data,
			size_t order,
			double alpha,
			size_t max_models,
			size_t mode,
			size_t num_threads,
			std::vector<std::ostream*> v_levels);
	Mbmdr(Data* data,
			Rcpp::List saved_mbmdr,
			size_t num_threads,
			std::vector<std::ostream*> v_levels);

	// Destructor
	virtual ~Mbmdr();

	// Function to fit MB-MDR models
	void fit();

	// Function to predict MB-MDR models
	std::vector<double> predict();

	Rcpp::List exportModels();

	// Getters
	uint getMode();
	size_t getOrder();
	size_t getN();
	double getAlpha();
	size_t getMaxModels();
	std::priority_queue<Model*, std::vector<Model*>, CompareModelPointers> getModels();
	std::unordered_set<std::string> getModelFeatureNames();

protected:

	// Pointer to original data
	Data* data;

	// Prediction mode
	uint mode;

	// Interaction order
	size_t order;

	// Total observations
	size_t n;

	// Significance threshold for cell labels
	double alpha;

	// For multithreading
	std::mutex mutex;
	size_t threads;

	// Vector of feature models
	size_t max_models;
	size_t num_models;
	std::priority_queue<Model*, std::vector<Model*>, CompareModelPointers> models;
	std::unordered_set<std::string> model_feature_names;

	// Current feature combination
	std::vector<size_t> feature_combination;

	// Predictions
	std::vector<double> predictions;

	// Add feature model to to priority queue, eventually replace model with lowest statistic
	void possiblyAdd(Model* new_model);

	// Iterate through all possible feature combinations
	bool getNextFeatureCombination(size_t j);

	// Verbose streams
	std::vector<std::ostream*> v_levels;

private:
	void fitModelInThread();
	void predictInThread();

	DISALLOW_COPY_AND_ASSIGN(Mbmdr);

};

#endif /* MBMDR_H_ */
