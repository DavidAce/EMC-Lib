
#include "DNA.hpp"



using namespace std;
using namespace EMC_constants;
using namespace Eigen;


std::ostream &operator<<(std::ostream &os, DNA const &genome) {
	for (int i = 0; i < nGenes; i++) {
//		os << genome.chromosomes[i] << endl;
        for (size_t nIndex = 0; nIndex < genome.chromosomes[i].size (); ++ nIndex) {
            os << genome.chromosomes[i][nIndex];
        }
    }
	return os;
}

bool DNA::operator==(const DNA &target) { //Compare DNA and return true if equal
	bool isequal = true;
	for (int i = 0; i < EMC_constants::nGenes; i++) {
		isequal = isequal && (chromosomes[i] == target.chromosomes[i]);
	}
	return isequal;
}
int DNA::operator()(int a) { 				//Return the bit at a.
		int gene = a / geneLength;
		int loci = a%geneLength;
		return chromosomes[gene][geneLength - loci - 1];
	}

void DNA::flip_loci(const int a) {
	//(loci)
	int gene = a / geneLength;
	int locus = a%geneLength;
    chromosomes[gene][geneLength - locus - 1]  =  !chromosomes[gene][geneLength - locus - 1];
//    chromosomes[gene].flip(geneLength - locus - 1);
}

//void DNA::flip_loci(Ref<ArrayXi> loci) {
//	//Flip all bits in array loci
//	int gene;// = a / geneLength;
//	int locus;// = a%geneLength;
//	for (int i = 0; i < loci.size(); i++) {
//		gene = loci(i)/geneLength;
//		locus = loci(i) % geneLength;
//		chromosomes[gene].flip(geneLength - locus - 1);
//	}
//}

void DNA::flip_loci(ArrayXi &loci) {
	//Flip all bits in array loci
	int gene;// = a / geneLength;
	int locus;// = a%geneLength;
	for (int i = 0; i < loci.size(); i++) {
		gene = loci(i)/geneLength;
		locus = loci(i) % geneLength;
//		chromosomes[gene].flip(geneLength - locus - 1);
        chromosomes[gene][geneLength - locus - 1]  =  !chromosomes[gene][geneLength - locus - 1];

    }
}

void DNA::copy_loci(const int a, const int bit) {
	//(loci,bit)
	bool b = bit == 1;
	int gene = a / geneLength;
	int locus = a%geneLength;
//	chromosomes[gene].set(geneLength - locus - 1, b);
    chromosomes[gene][geneLength - locus - 1] = b;
}

double DNA::bin2dec(const int i) {
    auto value = chromosomes[i].begin()._M_p;
//    auto value = chromosomes[i].to_ulong();
	return *value / (pow(2.0, geneLength) - 1) * (obj_fun.upper_bound(i) - obj_fun.lower_bound(i)) + obj_fun.lower_bound(i);
}

//bitset<geneLength> DNA::dec2bin(const int i) {
//	double value = parameters(i);
//	bitset <geneLength> A = (long long unsigned int) ((value - obj_fun.lower_bound(i)) / (obj_fun.upper_bound(i) - obj_fun.lower_bound(i))*(pow(2.0, geneLength) - 1));
//	return A;
//}



vector<bool> DNA::dec2bin(const int i) {
    long double value = parameters(i);
    unsigned long long int x = (unsigned long long int)((value - obj_fun.lower_bound(i)) / (obj_fun.upper_bound(i) - obj_fun.lower_bound(i))*(pow(2.0, geneLength) - 1));
    std::vector<bool> A;
    for(unsigned int i = 0; i < sizeof(x) * CHAR_BIT; ++i, x >>= 1) {
        A.push_back(x & 1);
    }
    std::reverse(A.begin(), A.end());
    //Truncate
    return vector<bool>(A.end() - geneLength, A.end());


//    std::string chars( std::bitset< sizeof(long) * CHAR_BIT >( x )
//                               .to_string( char(0), char(1) ) );
//    vector<bool> A = std::vector< bool >( chars.begin(), chars.end() );
//    size_t  nIndex;
//    for (nIndex = 0; nIndex < A.size (); ++ nIndex) {
//        cout << A[nIndex];
//    }
//    cout<< " nIndex = " << nIndex << endl;
//    return std::vector< bool >( chars.begin(), chars.end() );


//    while(modified_value) {
//        A.push_back(modified_value & 1);
//        modified_value >>= 1;
//    }



//    vector<bool> A = ((value - obj_fun.lower_bound(i)) / (obj_fun.upper_bound(i) - obj_fun.lower_bound(i))*(pow(2.0, geneLength) - 1));

//    bitset <geneLength> A = (long long unsigned int) ((value - obj_fun.lower_bound(i)) / (obj_fun.upper_bound(i) - obj_fun.lower_bound(i))*(pow(2.0, geneLength) - 1));
    //return A;
}


void DNA::set_parameter(const int i, const long double p) {
	//Set one parameter with a double
	parameters(i) = p;
	chromosomes[i] = dec2bin(i);
}
				
void DNA::update_parameters() {
	for (int i = 0; i < nGenes; i++) {
		parameters(i) = bin2dec(i);
	}
}

void DNA::set_parameters(const Array<long double, Dynamic, 1> &p) {
	//Set all parameters at once with an ArrayXd
	parameters = p;
	for (int i = 0; i < nGenes; i++) {
		chromosomes[i] = dec2bin(i);
	}
}

void DNA::randomize_dna(){
    chromosomes.resize(nGenes);
    parameters.resize(nGenes);
    for (int i = 0; i < nGenes; i++) {
        parameters(i) = uniform_double(obj_fun.lower_bound(i), obj_fun.upper_bound(i));
        chromosomes[i] = dec2bin(i);
    }
}



DNA::DNA(objective_function &ref):obj_fun(ref) {
    chromosomes.resize(nGenes);
    chromosomes.resize(nGenes);
	parameters.resize(nGenes);
	for (int i = 0; i < nGenes; i++) {
		parameters(i) = uniform_double(obj_fun.lower_bound(i),obj_fun.upper_bound(i));
		chromosomes[i] = dec2bin(i);

	}

}

DNA::DNA(objective_function &ref,bool ):obj_fun(ref) {
    chromosomes.resize(nGenes);
    parameters.resize(nGenes);
}
