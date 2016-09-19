#include "randomFunctions.hpp"


namespace EMC_rnd{
	std::mt19937 rng;
	std::uniform_int_distribution<>  rand_int_1(0, 1);
	std::uniform_real_distribution<> rand_real_1(0,1);
};
