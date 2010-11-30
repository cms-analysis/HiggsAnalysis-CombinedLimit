#include <iostream>
#include <exception>
#include <RooStats/HLFactory.h>

int main(int argc, char **argv) {
    using namespace std;

    if (argc != 2) { std::cerr << "Usage: testhlf file.hlf " << std::endl; return 1; }

    try {
        RooStats::HLFactory hlf("factory", argv[1]);
        RooWorkspace *w = hlf.GetWs();

        if (w == 0) {
            std::cerr << "Could not read HLF from file " <<  argv[1] << std::endl;
            return 2;
        }

        w->Print("V");

        return 0;
    } catch (const std::exception &ex) {
        std::cerr << "Error " << ex.what() << " when reading " << argv[1] << std::endl;
        return 2;
    } catch (...) {
        std::cerr << "Error when reading " << argv[1] << std::endl;
        return 2;
    }
}


