#include "rlm/rlm.hpp"

#include <cctype>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

int rlm_check();

namespace {

class UsageError : public std::runtime_error {
public:
    explicit UsageError(const std::string& message) : std::runtime_error(message) {}
};

void print_usage(std::ostream& out) {
    out << "usage:\n"
           "  rlm run continuum [flags]\n"
           "  rlm run lattice [flags]\n"
           "  rlm check\n"
           "  rlm --help\n"
           "\n"
           "shared flags:\n"
           "  --mesh N          theta mesh (default 16)\n"
           "  --samples N       disorder samples (default 1)\n"
           "  --seed N          RNG seed (default 1)\n"
           "  --solver mkl|gpu  eigensolver (default mkl)\n"
           "  --threads N       OpenMP threads over theta2; <=0 uses the OpenMP default\n"
           "  --out PATH        output file (default energy_Chern_No.txt)\n"
           "\n"
           "continuum flags:\n"
           "  --nphi N\n"
           "  --impurity X      impurity concentration\n"
           "  --quanta X        quanta concentration\n"
           "  --off-head N      <0 selects the historical default\n"
           "\n"
           "lattice flags:\n"
           "  --lx N\n"
           "  --ly N\n"
           "  --q N\n"
           "  --band N\n";
}

int parse_int(const std::string& flag, const std::string& text) {
    std::size_t used = 0;
    int value = 0;
    try {
        value = std::stoi(text, &used, 10);
    } catch (const std::exception&) {
        throw UsageError(flag + " expects an integer, got '" + text + "'");
    }
    if (used != text.size()) {
        throw UsageError(flag + " expects an integer, got '" + text + "'");
    }
    return value;
}

unsigned parse_unsigned(const std::string& flag, const std::string& text) {
    std::size_t start = 0;
    while (start < text.size() && std::isspace(static_cast<unsigned char>(text[start]))) {
        ++start;
    }
    if (start >= text.size() || text[start] == '-' || text[start] == '+') {
        throw UsageError(flag + " expects an unsigned integer, got '" + text + "'");
    }
    std::size_t used = 0;
    unsigned long value = 0;
    try {
        value = std::stoul(text, &used, 10);
    } catch (const std::exception&) {
        throw UsageError(flag + " expects an unsigned integer, got '" + text + "'");
    }
    if (used != text.size() || value > static_cast<unsigned long>(std::numeric_limits<unsigned>::max())) {
        throw UsageError(flag + " expects an unsigned integer, got '" + text + "'");
    }
    return static_cast<unsigned>(value);
}

double parse_double(const std::string& flag, const std::string& text) {
    std::size_t used = 0;
    double value = 0.0;
    try {
        value = std::stod(text, &used);
    } catch (const std::exception&) {
        throw UsageError(flag + " expects a number, got '" + text + "'");
    }
    if (used != text.size() || !std::isfinite(value)) {
        throw UsageError(flag + " expects a number, got '" + text + "'");
    }
    return value;
}

struct CommonOpts {
    int samples = 1;
    int threads = 0;
    rlm::Backend backend = rlm::Backend::Mkl;
    std::string out_path = "energy_Chern_No.txt";
};

bool parse_shared_flag(const std::string& name, const std::string& value, CommonOpts& opt, int& mesh, unsigned& seed) {
    if (name == "mesh") {
        mesh = parse_int("--mesh", value);
    } else if (name == "samples") {
        opt.samples = parse_int("--samples", value);
        if (opt.samples < 0) {
            throw UsageError("--samples must be >= 0");
        }
    } else if (name == "seed") {
        seed = parse_unsigned("--seed", value);
    } else if (name == "solver") {
        if (value == "mkl") {
            opt.backend = rlm::Backend::Mkl;
        } else if (value == "gpu") {
            opt.backend = rlm::Backend::Gpu;
        } else {
            throw UsageError("--solver must be mkl or gpu");
        }
    } else if (name == "threads") {
        opt.threads = parse_int("--threads", value);
    } else if (name == "out") {
        opt.out_path = value;
    } else {
        return false;
    }
    return true;
}

void read_flag(int argc, char** argv, int& i, std::string& name, std::string& value) {
    const std::string arg = argv[i];
    if (arg.size() < 3 || arg.compare(0, 2, "--") != 0) {
        throw UsageError("unexpected argument '" + arg + "'");
    }
    const std::size_t eq = arg.find('=');
    if (eq != std::string::npos) {
        name = arg.substr(2, eq - 2);
        value = arg.substr(eq + 1);
        if (name.empty() || value.empty()) {
            throw UsageError(arg + " missing value");
        }
        return;
    }
    name = arg.substr(2);
    if (name.empty()) {
        throw UsageError("missing flag name");
    }
    if (i + 1 >= argc) {
        throw UsageError("--" + name + " missing value");
    }
    value = argv[++i];
    if (value.compare(0, 2, "--") == 0) {
        throw UsageError("--" + name + " missing value");
    }
}

template <typename Params, typename MakeModel>
int run_samples(const Params& params_in, const CommonOpts& opt, MakeModel make) {
    std::ofstream out(opt.out_path, std::ios::out | std::ios::trunc);
    if (!out) {
        std::cerr << "cannot open " << opt.out_path << " for writing\n";
        return 1;
    }
    out << std::scientific << std::setprecision(16);
    for (int sample = 0; sample < opt.samples; ++sample) {
        Params params = params_in;
        params.sample = sample;
        auto model = make(params);
        auto solver = rlm::make_solver(opt.backend);
        const rlm::SampleResult result = rlm::run_sample(*model, *solver, opt.threads);
        const int ns = model->n_state();
        if (static_cast<int>(result.energy.size()) != ns || static_cast<int>(result.chern.size()) != ns) {
            throw std::runtime_error(model->name() + ": run_sample size does not match n_state");
        }
        double chern_sum = 0.0;
        for (double value : result.chern) {
            chern_sum += value;
        }
        std::cerr << model->name() << " sample " << sample << " chern_sum " << std::scientific << std::setprecision(16)
                  << chern_sum << '\n';
        for (int state = 0; state < ns; ++state) {
            out << sample << ' ' << result.energy[static_cast<std::size_t>(state)] << ' '
                << result.chern[static_cast<std::size_t>(state)] << '\n';
        }
        out.flush();
    }
    if (!out) {
        std::cerr << "failed writing " << opt.out_path << '\n';
        return 1;
    }
    return 0;
}

int run_continuum(int argc, char** argv) {
    rlm::ContinuumParams params;
    CommonOpts opt;
    for (int i = 3; i < argc; ++i) {
        std::string name;
        std::string value;
        read_flag(argc, argv, i, name, value);
        if (parse_shared_flag(name, value, opt, params.n_mesh, params.seed)) {
            continue;
        }
        if (name == "nphi") {
            params.n_phi = parse_int("--nphi", value);
        } else if (name == "impurity") {
            params.impurity_concentration = parse_double("--impurity", value);
        } else if (name == "quanta") {
            params.quanta_concentration = parse_double("--quanta", value);
        } else if (name == "off-head") {
            params.off_head = parse_int("--off-head", value);
        } else {
            throw UsageError("unknown flag --" + name);
        }
    }
    return run_samples(params, opt, [](const rlm::ContinuumParams& p) { return rlm::make_continuum(p); });
}

int run_lattice(int argc, char** argv) {
    rlm::LatticeParams params;
    CommonOpts opt;
    for (int i = 3; i < argc; ++i) {
        std::string name;
        std::string value;
        read_flag(argc, argv, i, name, value);
        if (parse_shared_flag(name, value, opt, params.n_mesh, params.seed)) {
            continue;
        }
        if (name == "lx") {
            params.lx = parse_int("--lx", value);
        } else if (name == "ly") {
            params.ly = parse_int("--ly", value);
        } else if (name == "q") {
            params.Q = parse_int("--q", value);
        } else if (name == "band") {
            params.band = parse_int("--band", value);
        } else {
            throw UsageError("unknown flag --" + name);
        }
    }
    return run_samples(params, opt, [](const rlm::LatticeParams& p) { return rlm::make_lattice(p); });
}

}  // namespace

int main(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(std::cout);
            return 0;
        }
    }
    if (argc < 2) {
        print_usage(std::cout);
        return 0;
    }
    try {
        const std::string command = argv[1];
        if (command == "check") {
            if (argc != 2) {
                throw UsageError("check takes no arguments");
            }
            return rlm_check();
        }
        if (command != "run") {
            throw UsageError("unknown command '" + command + "'");
        }
        if (argc < 3) {
            throw UsageError("run requires continuum or lattice");
        }
        const std::string model = argv[2];
        if (model == "continuum") {
            return run_continuum(argc, argv);
        }
        if (model == "lattice") {
            return run_lattice(argc, argv);
        }
        throw UsageError("unknown model '" + model + "'");
    } catch (const UsageError& ex) {
        std::cerr << ex.what() << '\n';
        print_usage(std::cerr);
        return 2;
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << '\n';
        return 1;
    }
}
