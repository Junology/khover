#include <iostream>
#include <exception>
#include <functional>
#include <vector>
#include <map>
#include <cctype>
#include <string>
#include <string_view>

#include "config.hpp"

class CommandOpts {
public:
    struct Command
    {
        std::string description;
        bool take_parameter;
        std::function<void(std::string_view const&)> func;
    };

private:
    std::map<char,Command> cmds{};
    std::map<std::string,char> lsmapper{};

public:
    CommandOpts() {
        cmds.emplace(
            'h',
            Command{"Show help.", false, [this](auto) { this->showHelp(); }}
            );
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
                std::move(func)});
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

    std::pair<bool,std::vector<std::string_view>>
    process(int argc, char** argv) {
        std::vector<std::string_view> args(argv, argv+argc);
        bool done_something = false;

        for(auto itr = std::next(std::begin(args)); itr != std::end(args);) {
            if(itr->size() <= 1 || (*itr)[0] != '-') {
                ++itr;
                continue;
            }

            if(std::isalnum((*itr)[1])) {
                if(auto cmd = cmds.find((*itr)[1]); cmd != std::end(cmds)) {
                    cmd->second.func(itr->substr(2));
                    itr = args.erase(itr);
                    done_something = true;
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
                        cmds[name->second].func(*std::next(itr));
                        itr = args.erase(itr, std::next(itr,2));
                    }
                    else {
                        cmds[name->second].func("");
                        itr = args.erase(itr);
                    }
                    done_something = true;
                    continue;
                }
            }

            ++itr;
        }

        return std::make_pair(done_something,std::move(args));
    }
};

int main(int argc, char* argv[])
{
    enum class AppMode {
        Khovanov,
        Crux,
        Derivative,
        CruxImage,
        Misc
    };

    AppMode mode = AppMode::Misc;
    int target_crossing = -1;

    CommandOpts opts;

    opts.addOpt(
        'v', "version", "Show version.", false,
        [](auto) { std::cout << appconf::version << std::endl; }
        );
    opts.addOpt(
        'c', "crux",
        "Compute the crux complex at the crossing <ARG>.",
        true,
        [&mode, &target_crossing](std::string_view const& str) {
            try {
                mode = AppMode::Crux;
                target_crossing = std::stoi(std::string(str));
            }
            catch(std::invalid_argument const& err) {
                std::cerr << "Invalid argument." << std::endl;
            }
            catch(std::out_of_range const& err) {
                std::cerr << "Out of range." << std::endl;
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
            }
            catch(std::invalid_argument const& err) {
                std::cerr << "Invalid argument." << std::endl;
            }
            catch(std::out_of_range const& err) {
                std::cerr << "Out of range." << std::endl;
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
            }
            catch(std::invalid_argument const& err) {
                std::cerr << "Invalid argument." << std::endl;
            }
            catch(std::out_of_range const& err) {
                std::cerr << "Out of range." << std::endl;
            }
        } );
    auto [done,args] = opts.process(argc, argv);

    // help/version are shown.
    if(done && mode == AppMode::Misc) {
        return 0;
    }

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

    // If no option is given, enter the Khovanov homology mode.
    if(!done)
        mode = AppMode::Khovanov;

    switch(mode) {
    case AppMode::Khovanov:
        std::cout << "Computation of Khovanov homology is WIP." << std::endl;
        break;

    case AppMode::Crux:
        std::cout << "Computation of Crux homology is WIP." << std::endl;
        std::cout << "Target crossing: " << target_crossing << std::endl;
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
