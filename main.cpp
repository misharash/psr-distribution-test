#include "params.hpp"
#include "pulsar.hpp"
#include <vector>
#include <random>
#include "create.hpp"
#include "dump.hpp"

int main(int argc, char** argv) {
	//random setup
	std::random_device rd;
	std::mt19937 e2(rd());
	std::uniform_real_distribution<> dist(0, 1);
    
	std::vector<Pulsar> p(Nstart); //main array
    
    dump_init();
    
    create_all(p, dist, e2); //create initial pulsars
    
    dump(p, 0, 0); //dump initial distribution
    
	return 0;
}
