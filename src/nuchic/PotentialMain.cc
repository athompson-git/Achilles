#include "nuchic/Version.hh"
#include "nuchic/Logging.hh"
#include "nuchic/RunPotential.hh"

#include "docopt.h"
#include "fmt/color.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

void Splash() {

    fmt::print(R"splash(
+=================================================================+
|                                                                 | 
|       d8b   db db    db  .o88b. db   db d888888b  .o88b.        |
|       888o  88 88    88 d8P  Y8 88   88   `88'   d8P  Y8        |
|       88V8o 88 88    88 8P      88ooo88    88    8P             |
|       88 V8o88 88    88 8b      88~~~88    88    8b             |
|       88  V888 88b  d88 Y8b  d8 88   88   .88.   Y8b  d8        |
|       VP   V8P ~Y8888P'  `Y88P' YP   YP Y888888P  `Y88P'        |
|                                                                 |
+-----------------------------------------------------------------+
|                                                                 |
|    Version: {:52}|
|    Authors: Joshua Isaacson, William Jay, Alessandro Lovato,    | 
|             Pedro A. Machado, Noemi Rocco                       | 
|                                                                 |
+=================================================================+

)splash", NUCHIC_VERSION);

    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::red),
               "NOTE: You are running in potential test mode. Hard interactions are not included\n");
}

static const std::string USAGE =
R"(
    Usage:
      nuchic-potential [<input>] [-v | -vv]
      nuchic-potential (-h | --help)
      nuchic-potential --version

    Options:
      -v[v]            Increase verbosity level.
      -h --help        Show this screen.
      --version        Show version.
)";

int main(int argc, char *argv[]) {

    Splash();
    std::map<std::string, docopt::value> args = docopt::docopt(USAGE,
                                                    { argv + 1, argv + argc },
                                                    true, // show help if requested
                                                    fmt::format("Nuchic {}", NUCHIC_VERSION)); //version string

    std::string runcard = "potential.yml";
    if(args["<input>"].isString()) runcard = args["<input>"].asString();

    auto verbosity = static_cast<int>(2 - args["-v"].asLong());
    CreateLogger(verbosity, 5);

    nuchic::RunPotential(runcard);

    return 0;
}