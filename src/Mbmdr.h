
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

	Mbmdr();

	Mbmdr(Data* data,
			size_t order,
			double alpha,
			size_t max_models,
			size_t mode,
			size_t num_threads,
			std::vector<std::ostream*> v_levels);

	virtual ~Mbmdr();

	void fit();

	void predict();

	Rcpp::List exportModels();

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
	std::priority_queue<Model*, std::vector<Model*>, CmpModelPtrs> models;

	// Current feature combination
	std::vector<size_t> feature_combination;

	// Add feature model to to priority queue, eventually replace model with lowest statistic
	void possiblyAdd(Model* new_model);

	// Iterate through all possible feature combinations
	bool getNextFeatureCombination(size_t j);

	// Verbose streams
	std::vector<std::ostream*> v_levels;

private:
	void fitModelInThread();

	DISALLOW_COPY_AND_ASSIGN(Mbmdr);

};

#endif /* MBMDR_H_ */
