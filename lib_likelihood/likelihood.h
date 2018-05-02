#include<string>
#include<sstream>

namespace likelihood {
double calculate_reconstruction_likelihood(std::string atoms_filename, std::string align_dir,
  std::stringstream& reconstruction_stream, int improvements = 0);
}
