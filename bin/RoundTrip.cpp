#include <map>
#include "TSystem.h"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/TFileIO.h"

namespace po = boost::program_options;

using namespace std;

int main(int argc, char* argv[]){
  string datacard   = "";
  string mass       = "";
  string output     = "";
  string output_file = "";
  bool print        = false;

  gSystem->Load("libHiggsAnalysisCombinedLimit.dylib");

  po::options_description help_config("Help");
  help_config.add_options()
    ("help,h", "produce help message");
  po::options_description config("Configuration");
  config.add_options()
    ("input,i",      po::value<string>(&datacard)->required(),
        "The datacard .txt file [REQUIRED]")
    ("output,o",       po::value<string>(&output)->required(),
        "Name of the output txt file to create [REQUIRED]")
    ("root,r ",       po::value<string>(&output_file)->required(),
        "Name of the output root file to create [REQUIRED]")
    ("mass,m",          po::value<string>(&mass)->default_value(""),
        "Signal mass point of the input datacard")
    ("print",    po::value<bool>(&print)->default_value(print)->implicit_value(true),
        "Print CombineHarvester instance");

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv)
                .options(help_config)
                .allow_unregistered()
                .run(),
            vm);
  po::notify(vm);
  if (vm.count("help")) {
    cout << config << "\n";
    cout << "Example usage: " << endl;
    cout << "RoundTrip -i htt_mt_125.txt -o htt_mt_125_test.txt -r "
            "htt_mt_125_test.root -m 125 --print \n";
    return 1;
  }
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  ch::CombineHarvester cmb;
  // cmb.SetVerbosity(2);

  cmb.ParseDatacard(datacard, "", "", "", 0, mass);

  if (print) cmb.PrintAll();

  TH1::AddDirectory(false);
  cmb.WriteDatacard(output.c_str(), output_file.c_str());

  return 0;
}

