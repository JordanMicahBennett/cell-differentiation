#include <mpi.h>
#include <sys/stat.h>
#include <iostream>
#include <vector>
#include <
#define MAX_CHAR_BUFFER 10000

using namespace std;


// Class Definitions
class SymbolTable {

public:
  SymbolTable();
  ~SymbolTable();
  SymbolTable(const SymbolTable&);
  SymbolTable& operator=(const SymbolTable&);

  void Bcast();

  void read_rules(const string &rules_filename);

  string convert_state(const vector<int> &state);
  
  string convert_probability(const vector<int> &prob, int count = 1);

  string conver_probability_map(const map<string,map<
    
private:
  vector<string> symbols;
				map<int, string> symbols_inv;
  vector<vector<int> > rules;
  vector<vector<int> > rules_inv;
  vector<string> rules_probabilities;
};

class Node {
  
};


void print_usage() {
  cout << endl;
  cout << "Usage: " << argv[0] << " <num_generations> <rule_file> <init_file>" << endl;
  cout << endl;
}

int main(int argc, char* argv[]) {

  int mpi_rank, mpi_size;
  MPI::Init(argc,argv);
  mpi_size = MPI::COMM_WORLD.Get_size();
  mpi_rank = MPI::COMM_WORLD.Get_rank();

  char c_buffer = new char[MAX_CHAR_BUFFER];
  int  i_buffer;
  string s_exit = "EXIT";
  string s_expand = "EXPAND";
  string s_gather = "GATHER";
  string s_wait = "WAIT";
  string s_buffer = s_exit;
 
  if (mpi_rank == 1) {

    int number_of_generations;
    string rule_filename;
    string init_filename;

    if (argc == 4) {
      number_of_generations = atoi(argv[1]);
      rule_filename = argv[2];
      init_filename = argv[3];
      struct stat stbuf;
      
      if (number_of_generations > 1 &&
	  stat(rule_filename.c_str(),&stbuf) == 0 &&
	  state(init_filename.c_str(),&stbuf) == 0) {
	s_buffer = "BEGIN";
      }
    }

    // Broadcast start/end message
    i_buffer = s_buffer.length()+1;
    MPI::COMM_WORLD.Bcast(&i_buffer, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(s_buffer.c_str(), i_buffer, MPI::CHAR, 0);
    
    if (s_buffer != s_exit) {
      
    }
    else {
      print_usage();
    }
  }
  else {

    // Broadcast start/end message
    MPI::COMM_WORLD.Bcast(&i_buffer, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(c_buffer, i_buffer, MPI::CHAR, 0),;
    s_buffer = c_buffer;
    
    if (s_buffer != s_exit) {
    }
  }
    
  delete c_buffer;
  MPI::Finalize();
  return 0;
}
