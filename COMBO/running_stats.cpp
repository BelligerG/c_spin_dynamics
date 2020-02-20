#include <iostream>
#include "running_stats.h"

RunningStats::RunningStats(){}

void RunningStats::UpdateStats(double value){
	//Take in a value and update the stats
	number_of_samples++;

	if (number_of_samples==1){
		old_mean = mean = value;
		var = 0;
		max = value;
		min = value;
	} else{
		mean = old_mean + (value - old_mean)/number_of_samples;
		var = var + (value - old_mean)*(value - mean);

		max = (value > max) ? value : max;
		min = (value < min) ? value : min;

		old_mean = mean;
	}
}


double RunningStats::GetMean(){ return (number_of_samples > 0) ? mean : 0.f; }
double RunningStats::GetVar(){ return (number_of_samples > 0) ? var/(number_of_samples-1): 0.f; }
double RunningStats::GetStdDev(){ return sqrt(GetVar()); }
double RunningStats::GetMax(){ return (number_of_samples > 0) ? max : 0.f; }
double RunningStats::GetMin(){ return (number_of_samples > 0) ? min : 0.f; }
int RunningStats::GetNumberOfSamples(){ return number_of_samples; }
