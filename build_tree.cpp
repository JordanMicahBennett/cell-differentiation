#include <mpi.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <map>

#define MAX_BUFFER 10000

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
    MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
    for (int x = 0; x < symbols.size(); x++) {
      *i_buffer = symbols[x].length() + 1;
      MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
      strcpy(c_buffer,symbols[x].c_str());
      MPI::COMM_WORLD.Bcast(c_buffer, *i_buffer, MPI::CHAR, 0);
    }  
    // Symbols Inv.
    *i_buffer = symbols_inv.size();
    MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
    for (map<string,int>::iterator x = symbols_inv.begin();
	 x != symbols_inv.end(); x++) {
      *i_buffer = x->first.length() + 1;
      MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
      strcpy(c_buffer,x->first.c_str());
      MPI::COMM_WORLD.Bcast(c_buffer, *i_buffer, MPI::CHAR, 0);
      *i_buffer = x->second;
      MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
    } 
    // Rules
    *i_buffer = rules.size();
    MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
    for (int x = 0; x < rules.size(); x++) {
      *i_buffer = rules[x].size();
      MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
      memcpy(i_buffer,&(rules[x][0]),sizeof(int)*rules[x].size());
      MPI::COMM_WORLD.Bcast(i_buffer, rules[x].size(), MPI::INT, 0);
    }
    // Rules_Inv
    *i_buffer = rules_inv.size();
    MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
    for (int x = 0; x < rules_inv.size(); x++) {
      *i_buffer = rules_inv[x].size();
      MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
      memcpy(i_buffer,&(rules_inv[x][0]),sizeof(int)*rules_inv[x].size());
      MPI::COMM_WORLD.Bcast(i_buffer, rules_inv[x].size(), MPI::INT, 0);
    }
    // Rules Probabilities
    *i_buffer = rules_probabilities.size();
    MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
    for (int x = 0; x < rules_probabilities.size(); x++) {
      *i_buffer = rules_probabilities[x].length() + 1;
      MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
      strcpy(c_buffer, rules_probabilities[x].c_str());
      MPI::COMM_WORLD.Bcast(c_buffer, *i_buffer, MPI::CHAR, 0);
    }
    return;
  }

  void Bcast(int* i_buffer, char* c_buffer) {
    int i,j;
    // Symbols
    MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
    i = *i_buffer;
    for (int x = 0; x < i; x++) {
      MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
      MPI::COMM_WORLD.Bcast(c_buffer, *i_buffer, MPI::CHAR, 0);
      symbols.push_back(string(c_buffer));
    }  
    // Symbols Inv.
    MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
    i = *i_buffer;
    for (int x = 0; x < i; x++) {
      string symbol;
      MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
      MPI::COMM_WORLD.Bcast(c_buffer, *i_buffer, MPI::CHAR, 0);
      symbol = c_buffer;
      MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
      symbols_inv.insert(pair<string,int>(symbol,*i_buffer));
    } 
    // Rules
    MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
    i = *i_buffer;
    for (int x = 0; x < i; x++) {
      MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
      j = *i_buffer;
      MPI::COMM_WORLD.Bcast(i_buffer, j, MPI::INT, 0);
      rules.push_back(vector<int>());
      for (int y = 0; y < j; y++) {
	rules[x].push_back(i_buffer[y]);
      }
    }
    // Rules_Inv
    MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
    i = *i_buffer;
    for (int x = 0; x < i; x++) {
      MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
      j = *i_buffer;
      MPI::COMM_WORLD.Bcast(i_buffer, j, MPI::INT, 0);
      rules_inv.push_back(vector<int>());
      for (int y = 0; y < rules_inv[x].size(); y++) {
	rules_inv[x].push_back(i_buffer[y]);
      }
    }
    // Rules Probabilities
    MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
    i = *i_buffer;
    for (int x = 0; x < rules_probabilities.size(); x++) {
      MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
      MPI::COMM_WORLD.Bcast(c_buffer, *i_buffer, MPI::CHAR, 0);
      rules_probabilities.push_back(string(c_buffer));
    }
    return;
  }

};

struct Node {
  vector<int> state;
  vector<int> expandable;
  vector<int> selected;
  SymbolTable* symbol_table;

  void Send(int dest, int* i_buffer, char* c_buffer) {
    // State
    *i_buffer = state.size();
    MPI::COMM_WORLD.Send(i_buffer, 1, MPI::INT, dest, 0);
    memcpy(i_buffer, &(state[0]), sizeof(int)*state.size());
    MPI::COMM_WORLD.Send(i_buffer, state.size(), MPI::INT, dest, 0);
    // Expandable
    *i_buffer = expandable.size();
    MPI::COMM_WORLD.Send(i_buffer, 1, MPI::INT, dest, 0);
    memcpy(i_buffer, &(expandable[0]), sizeof(int)*expandable.size());
    MPI::COMM_WORLD.Send(i_buffer, expandable.size(), MPI::INT, dest, 0);
    // Selected
    *i_buffer = selected.size();
    MPI::COMM_WORLD.Send(i_buffer, 1, MPI::INT, dest, 0);
    memcpy(i_buffer, &(selected[0]), sizeof(int)*selected.size());
    MPI::COMM_WORLD.Send(i_buffer, selected.size(), MPI::INT, dest, 0);
  }

  void Recv(int source, int* i_buffer, char* c_buffer) {
    int i;
    // State
    MPI::COMM_WORLD.Recv(i_buffer, 1, MPI::INT, source, 0);
    i = *i_buffer;
    MPI::COMM_WORLD.Recv(i_buffer, i, MPI::INT, source, 0);
    for (int x = 0; x < i; x ++)
      state.push_back(i_buffer[x]);
    // Expandable
    MPI::COMM_WORLD.Recv(i_buffer, 1, MPI::INT, source, 0);
    i = *i_buffer;
    MPI::COMM_WORLD.Recv(i_buffer, i, MPI::INT, source, 0);
    for (int x = 0; x < i; x ++)
      expandable.push_back(i_buffer[x]);
    // Selected
    MPI::COMM_WORLD.Recv(i_buffer, 1, MPI::INT, source, 0);
    i = *i_buffer;
    MPI::COMM_WORLD.Recv(i_buffer, i, MPI::INT, source, 0);
    for (int x = 0; x < i; x ++)
      selected.push_back(i_buffer[x]);
  }

};

int main(int argc, char* argv[]) {

  int mpi_rank, mpi_size;
  MPI::Init(argc,argv);
  mpi_size = MPI::COMM_WORLD.Get_size();
  mpi_rank = MPI::COMM_WORLD.Get_rank();

  char* c_buffer = new char[MAX_BUFFER];
  int*  i_buffer = new int[MAX_BUFFER];
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
    *i_buffer = s_buffer.length()+1;
    strcpy(c_buffer,s_buffer.c_str());
    MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(c_buffer, *i_buffer, MPI::CHAR, 0);
    
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
    MPI::COMM_WORLD.Bcast(i_buffer, 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(c_buffer, *i_buffer, MPI::CHAR, 0);
    s_buffer = c_buffer;
    
    if (s_buffer != s_exit) {
      cout << "Process " << mpi_rank << " alive" << endl;
    }
  }
    
  delete [] c_buffer;
  delete [] i_buffer;
  MPI::Finalize();
  return 0;
}
