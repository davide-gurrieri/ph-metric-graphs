#include "GetPot"
#include "network_analysis.h"

int main(int argc, char *argv[]) {
  GetPot command_line(argc, argv);
  const string file_name = command_line("file", "../options1.txt");
  GetPot datafile(file_name.c_str());

  string network_name = datafile("network_name", "");
  string vertices_file_name = datafile("vertices_file_name", "../");
  string edges_file_name = datafile("edges_file_name", "../");
  char vertex_sep = datafile("vertex_sep", ',');
  char edge_sep = datafile("edge_sep", ',');
  bool indexes_from_0 = datafile("indexes_from_0", 0);

  NetworkData data;
  data.parse(vertices_file_name, vertex_sep, edges_file_name, edge_sep,
             indexes_from_0);

  data.print_network_info();

  NetworkAnalysis net(&data);

  net.analyse();

  net.print_global_analysis();

  net.save_persistence(network_name);

  return 0;
}