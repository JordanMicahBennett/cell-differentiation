#include <mpi.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <map>

#define MAX_CHAR_BUFFER 10000

using namespace std;

struct SymbolTable {
  vector<string> symbols;
  map<string,int> symbols_inv;
  vector<vector<int> > rules;
  vector<vector<int > > rules_inv;
  vector<string> rules_probabilities;

  void Bcast_Root(int* i_buffer, char* c_buffer) {
    // Symbols
    *i_buffer = symbols.size();
    MPI::COMM_WORLD.Bcast(&i_buffer, 1, MPI::INT, 0);
    for (int x = 0; x < symbols.size(); x++) {
      *i_buffer = symbols[x].length() + 1;
      MPI::COMM_WORLD.Bcast(&i_buffer, 1, MPI::INT, 0);
      strcpy(c_buffer,symbols[x].c_str());
      MPI::COMM_WORLD.Bcast(c_buffer, *i_buffer, MPI::CHAR, 0);
    }  
    // Symbols Inv.
    *i_buffer = symbols_inv.size();
    MPI::COMM_WORLD.Bcast(&i_buffer, 1, MPI::INT, 0);
    for (map<string,int>::iterator x = symbols_inv.begin();
	 x != symbols_inv.end(); x++) {
      *i_buffer = x->first.length() + 1;
      MPI::COMM_WORLD.Bcast(&i_buffer, 1, MPI::INT, 0);
      strcpy(c_buffer,x->first.c_str());
      MPI::COMM_WORLD.Bcast(c_buffer, *i_buffer, MPI::CHAR, 0);
      *i_buffer = x->second;
      MPI::COMM_WORLD.Bcast(&i_buffer, 1, MPI::INT, 0);
    } 
  }
};

struct Node {
  vector<int> state;
  vector<int> expandable;
  vector<int> selected;
  SymbolTable* symbol_table;
};

int main(int argc, char* argv[]) {

  int mpi_rank, mpi_size;
  MPI::Init(argc,argv);
  mpi_size = MPI::COMM_WORLD.Get_size();
  mpi_rank = MPI::COMM_WORLD.Get_rank();

  char *c_buffer = new char[MAX_CHAR_BUFFER];
  int  i_buffer;
  string s_exit = "EXIT";
  string s_expand = "EXPAND";
  string s_gather = "GATHER";
  string s_wait = "WAIT";
  string s_buffer = s_exit;
 
  if (mpi_rank == 0) {

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
	  stat(init_filename.c_str(),&stbuf) == 0) {
	s_buffer = "BEGIN";
      }
    }

    // Broadcast start/end message
    i_buffer = s_buffer.length()+1;
    strcpy(c_buffer,s_buffer.c_str());
    MPI::COMM_WORLD.Bcast(&i_buffer, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(c_buffer, i_buffer, MPI::CHAR, 0);
    
    if (s_buffer != s_exit) {
      cout << "Process " << mpi_rank << " alive" << endl;      
    }
    else {
      cout << endl;
      cout << "Usage: " << argv[0] << " <num_generations> <rule_file> <init_file>" << endl;
      cout << endl;
    }
  }
  else {

    // Broadcast start/end message
    MPI::COMM_WORLD.Bcast(&i_buffer, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(c_buffer, i_buffer, MPI::CHAR, 0);
    s_buffer = c_buffer;
    
    if (s_buffer != s_exit) {
      cout << "Process " << mpi_rank << " alive" << endl;
    }
  }
    
  delete [] c_buffer;
  MPI::Finalize();
  return 0;
}
