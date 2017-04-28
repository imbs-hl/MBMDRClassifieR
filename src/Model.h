
#ifndef MODEL_H_
#define MODEL_H_

#include "globals.h"
#include "Data.h"

class Model {

public:

	// Constructors
	Model();
	Model(Data* data,
			size_t order,
			std::vector<size_t> features,
			double alpha,
			std::vector<std::ostream*> v_levels);

	// Destructor
	virtual ~Model();

	// Function to fit feature model
	virtual void fit() = 0;

	// Function to predict outcome given a fitted feature model
	std::vector<double> predict();

	// Construct model from saved object
	void loadModel(std::vector<uint> in_cell,
			std::vector<uint> out_cell,
			std::vector<double> cell_predictions,
			std::vector<double> cell_statistics,
			std::vector<double> cell_pvalues,
			std::vector<int> cell_labels,
			double statistic,
			double pvalue);

	// Getters
	std::vector<uint> getObservationsInCell() const;
	std::vector<uint> getObservationsOutCell() const;
	std::vector<double> getCellPredictions() const;
	std::vector<double> getCellStatistics() const;
	std::vector<double> getCellPValues() const;
	std::vector<int> getCellLabels() const;
	double getModelStatistic() const;
	double getModelPValue() const;
	size_t getOrder() const;
	std::vector<size_t> getFeatures() const;
	size_t getNumObservations() const;
	double getAlpha() const;

protected:

	// Pointer to original data
	Data* data;

	// Interaction order
	size_t order;

	// Features used in this feature model
	std::vector<size_t> features;

	// Total observations
	size_t n;

	// Significance threshold for cell labels
	double alpha;

	// Total observations per genotype combination
	std::vector<uint> in_cell;

	// Total observations in all other genotype combinations per genotype combination
	std::vector<uint> out_cell;

	// Cell predictions
	std::vector<double> cell_predictions;

	// Test statistic per cell
	std::vector<double> cell_statistics;

	// P value per cell
	std::vector<double> cell_pvalues;

	// Cell label
	std::vector<int> cell_labels;

	// Test statistic of model
	double statistic;

	// P value of model
	double pvalue;

	// Get counts per cell
	virtual void getCounts() = 0;

	// Classify cells
	virtual void classifyCells() = 0;

	// Calculate model test statistic
	virtual void calculateModelTestStatistic() = 0;

	// Verbose streams
	std::vector<std::ostream*> v_levels;

private:
	DISALLOW_COPY_AND_ASSIGN(Model);

};

// Comparisons of feature models
class CompareModelPointers {
	bool reverse;
public:
	CompareModelPointers(const bool& revparam = false) {
		reverse = revparam;
	}
	bool operator() (const Model* lhs, const Model* rhs) const {
		if (reverse) {
			return (lhs->getModelStatistic() < rhs->getModelStatistic());
		} else {
			return (lhs->getModelStatistic() > rhs->getModelStatistic());
		}
	}
};

#endif /* MODEL_H_ */
