#ifndef RUNNING_STATS_H
#define RUNNING_STATS_H

#include <cmath>

class RunningStats
{
public:

	RunningStats();
	void UpdateStats(double);

	double GetMean();
	double GetVar();
	double GetStdDev();
	double GetMax();
	double GetMin();
	int GetNumberOfSamples();

private:
	
	double mean, var, max, min, old_mean;
	int number_of_samples=0;
};
#endif
