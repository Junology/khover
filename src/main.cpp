#include <iostream>
#include <iomanip>
#include <exception>
#include <functional>
#include <vector>
#include <map>
#include <cctype>
#include <string>
#include <string_view>

#include "config.hpp"
#include "linkdiagram.hpp"
#include "khovanov.hpp"

#include "debug/debug.hpp"

using namespace khover;
/*************************************
 *** Command line option processor ***
 *************************************/
/*!
 * A simple class to manipulate commandline arguments.
 */
class CommandOpts {
public:
    enum ProcessResult : unsigned int {
        Pass = 0b00,
        Error = 0b01,
        Terminate = 0b10,
        TerminateSuccess = Terminate | Pass,
        TerminateError = Terminate | Error,
    };

    struct Command
    {
        std::string description;
        bool take_parameter;
        std::function<ProcessResult(std::string_view const&)> func;
    };

private:
    std::map<char,Command> cmds{};
    std::map<std::string,char> lsmapper{};

public:
    CommandOpts() {
        cmds.emplace(
            'h',
            Command{
                "Show help.", false,
                [this](auto) -> ProcessResult {
                    this->showHelp();
                    return TerminateSuccess;
                }
            });
        lsmapper.emplace("help", 'h');
    }

    void addOpt(
        char sname,
        std::string lname,
        std::string const &description,
        bool take_parameter,
        decltype(Command::func) &&func
        )
    {
        cmds.emplace(
            sname,
            Command{
                description,
                take_parameter,
                std::move(func)
            });
        lsmapper.emplace(lname, sname);
    }

    void showHelp()
    {
        std::cout << "Usage: "
                  << appconf::name << " [OPTION] GAUSSCODE"
                  << std::endl;
        std::cout << appconf::brief << std::endl;
        std::cout << "GAUSSCODE must be a valid Gauss code of knots or links enclosed by the brackets '[' and ']'. For links, codes for different components are separated by '0'." << std::endl;
        std::cout << std::endl;
        std::cout << "List of options" << std::endl;
        for(auto& opt : lsmapper) {
            std::cout << "  -" << opt.second
                      << (cmds[opt.second].take_parameter ? "<ARG>" : "")
                      << ", --" << opt.first
                      << (cmds[opt.second].take_parameter ? " <ARG>" : "")
                      << std::endl
                      << "\t" << cmds[opt.second].description << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Project home page: " << appconf::url << std::endl;
    }

    std::pair<ProcessResult,std::vector<std::string_view>>
    process(int argc, char** argv) {
        std::vector<std::string_view> args(argv, argv+argc);
        ProcessResult terminate_flag = ProcessResult::Pass;

        for(auto itr = std::next(std::begin(args)); itr != std::end(args);) {
            if(itr->size() <= 1 || (*itr)[0] != '-') {
                ++itr;
                continue;
            }

            if(std::isalnum((*itr)[1])) {
                if(auto cmd = cmds.find((*itr)[1]); cmd != std::end(cmds)) {
                    terminate_flag = cmd->second.func(itr->substr(2));
                    itr = args.erase(itr);

                    if(terminate_flag & ProcessResult::Terminate)
                        return std::make_pair(terminate_flag, std::move(args));

                    continue;
                }
            }
            else if((*itr)[1] == '-') {
                if(auto name = lsmapper.find(std::string(itr->substr(2)));
                   name != std::end(lsmapper))
                {
                    if(cmds[name->second].take_parameter
                       && std::next(itr) != std::end(args))
                    {
                        terminate_flag = cmds[name->second].func(*std::next(itr));
                        itr = args.erase(itr, std::next(itr,2));
                    }
                    else {
                        terminate_flag = cmds[name->second].func("");
                        itr = args.erase(itr);
                    }

                    if(terminate_flag & ProcessResult::Terminate)
                        return std::make_pair(terminate_flag, std::move(args));
                    continue;
                }
            }

            ++itr;
        }

        return std::make_pair(ProcessResult::Pass,std::move(args));
    }
};


/*************************
 *** Utility functions ***
 *************************/
std::vector<int> read_vector(
    std::string_view const& strview, char delim)
    noexcept(false)
{
    using size_type = std::string_view::size_type;

    std::vector<int> result{};
    size_type beg = 0, len;

    do {
        len = strview.substr(beg).find(delim);
        result.push_back(std::stoi(std::string(strview.substr(beg,len))));
        beg += len+1;
    } while(len != std::string_view::npos);

    return result;
}


/*********************
 *** Main function ***
 *********************/
int main(int argc, char* argv[])
{
    enum class AppMode {
        Khovanov,
        Crux,
        Derivative,
        CruxImage,
    };

    AppMode mode = AppMode::Khovanov;
    int target_crossing = -1;
    std::vector<std::pair<std::size_t,bool>> is_positive_specifiers;

    CommandOpts opts;

    opts.addOpt(
        'v', "version", "Show version.", false,
        [](auto) {
            std::cout << appconf::version << std::endl;
            return CommandOpts::TerminateSuccess;
        });
    opts.addOpt(
        'p', "positive",
        "The crossing <ARG> becomes positive. If several -p and -n options are supplied, the latter one has higher priority.",
        true,
        [&is_positive_specifiers](std::string_view const& str) {
            try {
                std::size_t c = static_cast<std::size_t>(
                    std::abs(std::stoi(std::string(str))));
                is_positive_specifiers.emplace_back(c, true);
                return CommandOpts::Pass;
            }
            catch(std::invalid_argument const& err) {
                std::cerr << "ERROR: Invalid argument." << std::endl;
                return CommandOpts::TerminateError;
            }
            catch(std::out_of_range const& err) {
                std::cerr << "ERROR: Out of range." << std::endl;
                return CommandOpts::TerminateError;
            }
        });
    opts.addOpt(
        'n', "negative",
        "The crossing <ARG> becomes negative. If several -p and -n options are supplied, the latter one has higher priority.",
        true,
        [&is_positive_specifiers](std::string_view const& str) {
            try {
                std::size_t c = static_cast<std::size_t>(
                    std::abs(std::stoi(std::string(str))));
                is_positive_specifiers.emplace_back(c, false);
                return CommandOpts::Pass;
            }
            catch(std::invalid_argument const& err) {
                std::cerr << "ERROR: Invalid argument." << std::endl;
                return CommandOpts::TerminateError;
            }
            catch(std::out_of_range const& err) {
                std::cerr << "ERROR: Out of range." << std::endl;
                return CommandOpts::TerminateError;
            }
        });
    opts.addOpt(
        'c', "crux",
        "Compute the crux complex at the crossing <ARG>.",
        true,
        [&mode, &target_crossing](std::string_view const& str) {
            try {
                mode = AppMode::Crux;
                target_crossing = std::stoi(std::string(str));
                return CommandOpts::Pass;
            }
            catch(std::invalid_argument const& err) {
                std::cerr << "ERROR: Invalid argument." << std::endl;
                return CommandOpts::TerminateError;
            }
            catch(std::out_of_range const& err) {
                std::cerr << "ERROR: Out of range." << std::endl;
                return CommandOpts::TerminateError;
            }
        } );
    opts.addOpt(
        'd', "derivative",
        "Compute the first Vassiliev derivative at the crossing <ARG>.",
        true,
        [&mode, &target_crossing](std::string_view const& str) {
            try {
                mode = AppMode::Derivative;
                target_crossing = std::stoi(std::string(str));
                return CommandOpts::Pass;
            }
            catch(std::invalid_argument const& err) {
                std::cerr << "ERROR: Invalid argument." << std::endl;
                return CommandOpts::TerminateError;
            }
            catch(std::out_of_range const& err) {
                std::cerr << "ERROR: Out of range." << std::endl;
                return CommandOpts::TerminateError;
            }
        } );
    opts.addOpt(
        'i', "image",
        "Compute the crux image at the crossing <ARG>.",
        true,
        [&mode, &target_crossing](std::string_view const& str) {
            try {
                mode = AppMode::CruxImage;
                target_crossing = std::stoi(std::string(str));
                return CommandOpts::Pass;
            }
            catch(std::invalid_argument const& err) {
                std::cerr << "ERROR: Invalid argument." << std::endl;
                return CommandOpts::TerminateError;
            }
            catch(std::out_of_range const& err) {
                std::cerr << "ERROR: Out of range." << std::endl;
                return CommandOpts::TerminateError;
            }
        } );

    auto [procres,args] = opts.process(argc, argv);

    // Check if command line argument processor was terminated.
    if(procres & CommandOpts::Terminate) {
        return (procres & CommandOpts::Error) ? EXIT_FAILURE : EXIT_SUCCESS;
    }

    // Find a Gauss code in command line arguments.
    auto code_itr = std::find_if(
        std::next(std::begin(args)), std::end(args),
        [](std::string_view const& arg) {
            return arg.size() >= 2
            && arg.front() == '[' && arg.back() == ']';
        } );

    // No Gauss code given.
    if (code_itr == std::end(args)) {
        std::cerr << "ERROR: No Gauss code given." << std::endl;
        std::cerr << std::endl;
        opts.showHelp();
        return -1;
    }

    // Remove brackets '[' and ']' from the Gauss code.
    code_itr->remove_prefix(1);
    code_itr->remove_suffix(1);

    // Read Gauss code
    std::optional<LinkDiagram> diagram;
    try {
        auto gcode = read_vector(*code_itr, ',');
        diagram = read_gauss_code(gcode, is_positive_specifiers);
        if (!diagram)
            throw std::invalid_argument("Invalid Gauss code");
    }
    catch(std::invalid_argument const& err) {
        std::cerr << "ERROR: Invalid Gauss code." << std::endl;
        return EXIT_FAILURE;
    }
    catch(std::out_of_range const& err) {
        std::cerr << "ERROR: Index out of range." << std::endl;
        return EXIT_FAILURE;
    }

    DBG_MSG("Diagram successfully loaded.");

    switch(mode) {
    case AppMode::Khovanov: {
        int qmin =
            - static_cast<int>(diagram->smoothing(0u).first)
            - static_cast<int>(diagram->nnegative())
            + diagram->writhe();
        int qmax =
            static_cast<int>(diagram->smoothing(~state_t(0u)).first)
            + static_cast<int>(diagram->npositive())
            + diagram->writhe();
        for(int q = qmin; q <= qmax; q+=2) {
            auto ch = khChain(
                *diagram,
                SmoothCube::fromDiagram(*diagram),
                q);
            if(!ch)
                continue;
            std::cout << "q-degree: " << q << std::endl;
            int i = ch->mindeg();
            for(auto h : ch->compute()) {
                std::cout << std::setw(4) << std::right << (-i) << ": "
                          << std::resetiosflags(std::ios_base::adjustfield | std::ios_base::basefield)
                          << h.compute().pretty()
                          << std::endl;
                ++i;
            }
        }
    }
        break;

    case AppMode::Crux: {
        if (target_crossing <= 0 || target_crossing > static_cast<int>(diagram->ncrosses())) {
            std::cerr << "Invalid double point index." << std::endl;
            return EXIT_FAILURE;
        }

        int qmin =
            - static_cast<int>(diagram->smoothing(0u).first)
            - static_cast<int>(diagram->nnegative())
            + diagram->writhe();
        int qmax =
            static_cast<int>(diagram->smoothing(~state_t(0u)).first)
            + static_cast<int>(diagram->npositive())
            + diagram->writhe();

        for(int q = qmin; q <= qmax; q+=2) {
            auto ch = cruxChain(
                *diagram,
                target_crossing-1,
                CruxCube::fromDiagram(*diagram, target_crossing-1),
                q);
            if(!ch)
                continue;
            std::cout << "q-degree: " << q << std::endl;
            int i = ch->mindeg();
            for(auto h : ch->compute()) {
                std::cout << std::setw(4) << std::right << (-i) << ": "
                          << std::resetiosflags(std::ios_base::adjustfield | std::ios_base::basefield)
                          << h.compute().pretty()
                          << std::endl;
                ++i;
            }
        }
    }
        break;

    case AppMode::Derivative:
        std::cout << "Computation of the first Vassiliev derivative is WIP." << std::endl;
        std::cout << "Target crossing: " << target_crossing << std::endl;
        break;

    case AppMode::CruxImage:
        std::cout << "Computation of the crux image is WIP." << std::endl;
        std::cout << "Target crossing: " << target_crossing << std::endl;
        break;

    default:
        break;
    }

    return 0;
}
