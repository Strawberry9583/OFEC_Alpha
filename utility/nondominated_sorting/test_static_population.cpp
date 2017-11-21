#include "static_population.h"
#include <fstream>

int main() {
	NDS::uniform_population u3(2, 30, 0.5);
	auto data = u3.generate_new(3);
	std::ofstream outputfile("E://uniform_population_2_30_3.csv");
	for (auto& row : data)
		outputfile << row[0] << ", " << row[1] << std::endl;
	outputfile.close();
}