// InitializeParticles.cpp
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <cctype>
#include <cstdlib>
#include <cstdint>

// ----------------------------------------------------------------------------------
// Distribution struct & diameter generator
// ----------------------------------------------------------------------------------
struct Distribution {
    std::string type = "mono";
    double d       = 1.0;
    double mu      = 1.0, sigma   = 0.2;
    double mu1     = 1.0, sigma1  = 0.2;
    double mu2     = 2.0, sigma2  = 0.3;
    double p       = 0.5;
    double d1      = 1.0, d2     = 2.0;
    double d_min   = 0.5, d_max  = 1.5;
    double exponent= -2.5;
    double scale   = 1.0, shape  = 2.0;
    std::vector<double> custom;
};

std::vector<double> generate_diameter_distribution(
    size_t N,
    const Distribution& dist,
    std::mt19937& rng)
{
    std::vector<double> D(N);
    std::uniform_real_distribution<> ur(0.0,1.0);

    if (dist.type == "mono") {
        std::fill(D.begin(), D.end(), dist.d);
    }
    else if (dist.type == "gaussian") {
        std::normal_distribution<> nd(dist.mu, dist.sigma);
        for (auto &v : D) v = std::abs(nd(rng));
    }
    else if (dist.type == "bigaussian") {
        size_t N1 = static_cast<size_t>(std::round(dist.p * N));
        std::normal_distribution<> g1(dist.mu1, dist.sigma1);
        std::normal_distribution<> g2(dist.mu2, dist.sigma2);
        for (size_t i = 0; i < N1; ++i) D[i] = std::abs(g1(rng));
        for (size_t i = N1; i < N; ++i) D[i] = std::abs(g2(rng));
        std::shuffle(D.begin(), D.end(), rng);
    }
    else if (dist.type == "bidisperse") {
        size_t N1 = static_cast<size_t>(std::round(dist.p * N));
        for (size_t i = 0; i < N1; ++i) D[i] = dist.d1;
        for (size_t i = N1; i < N; ++i) D[i] = dist.d2;
        std::shuffle(D.begin(), D.end(), rng);
    }
    else if (dist.type == "lognormal") {
        std::lognormal_distribution<> ld(dist.mu, dist.sigma);
        for (auto &v : D) {
            double tmp = ld(rng);
            v = std::clamp(tmp, 0.0025, 30.0);
        }
    }
    else if (dist.type == "flat") {
        std::uniform_real_distribution<> ud(dist.d_min, dist.d_max);
        for (auto &v : D) v = ud(rng);
        // --- Rescale to ensure min(D) = d_min and max(D) = d_max ---
        double D_min = *std::min_element(D.begin(), D.end());
        double D_max = *std::max_element(D.begin(), D.end());

        if (D_max > D_min + 1e-12) {  // avoid divide-by-zero
            for (size_t i = 0; i < N; ++i) {
                D[i] = dist.d_min + (D[i] - D_min) * (dist.d_max - dist.d_min) / (D_max - D_min);
            }
        } else {
            std::fill(D.begin(), D.end(), dist.d_min);  // fallback: all equal
        }          
    }
    else if (dist.type == "powerlaw") {
        std::uniform_real_distribution<> ud(0.0,1.0);
        double a = dist.exponent + 1.0;
        for (size_t i = 0; i < N; ++i) {
            double u = ud(rng);
            if (std::abs(a) < 1e-12) {
                D[i] = dist.d_min * std::exp(u * std::log(dist.d_max / dist.d_min));
            } else {
                double top = std::pow(dist.d_max, a) - std::pow(dist.d_min, a);
                D[i] = std::pow(top * u + std::pow(dist.d_min, a), 1.0 / a);
            }
        }

        // --- Rescale to ensure min(D) = d_min and max(D) = d_max ---
        double D_min = *std::min_element(D.begin(), D.end());
        double D_max = *std::max_element(D.begin(), D.end());

        if (D_max > D_min + 1e-12) {  // avoid divide-by-zero
            for (size_t i = 0; i < N; ++i) {
                D[i] = dist.d_min + (D[i] - D_min) * (dist.d_max - dist.d_min) / (D_max - D_min);
            }
        } else {
            std::fill(D.begin(), D.end(), dist.d_min);  // fallback: all equal
        }        

    }
    else if (dist.type == "exponential") {
        double lambda = -std::log(0.5) / (dist.d_max - dist.d_min);
        std::uniform_real_distribution<> ud(0.0,1.0);
        for (size_t i = 0; i < N; ++i) {
            double u = ud(rng);
            D[i] = dist.d_min - (1.0/lambda) * std::log(1.0 - u);
            if (D[i] > dist.d_max) D[i] = dist.d_max;
        }

        // --- Rescale to ensure min(D) = d_min and max(D) = d_max ---
        double D_min = *std::min_element(D.begin(), D.end());
        double D_max = *std::max_element(D.begin(), D.end());

        if (D_max > D_min + 1e-12) {  // avoid divide-by-zero
            for (size_t i = 0; i < N; ++i) {
                D[i] = dist.d_min + (D[i] - D_min) * (dist.d_max - dist.d_min) / (D_max - D_min);
            }
        } else {
            std::fill(D.begin(), D.end(), dist.d_min);  // fallback: all equal
        }           
    }
    else if (dist.type == "weibull") {
        std::weibull_distribution<> wd(dist.shape, dist.scale);
        for (auto &v : D) v = wd(rng);
    }
    else if (dist.type == "custom") {
        if (dist.custom.size() != N)
            throw std::runtime_error("custom distribution needs N entries");
        D = dist.custom;
    }
    else {
        std::cerr << "Unsupported distribution: " << dist.type << std::endl;
        std::exit(1);
    }
    return D;
}

std::pair<std::vector<double>, double> scale_diametersND(
    const std::vector<double>& D,
    double phi_target,
    const std::vector<double>& box,
    size_t Ndim)
{
    double V_box = 1.0;
    for (double L : box) V_box *= L;
    double sumVol = 0.0;
    for (double Di : D) {
        double r = Di * 0.5;
        double vol_unit = std::pow(M_PI, Ndim/2.0) / std::tgamma(Ndim/2.0 + 1.0);
        sumVol += vol_unit * std::pow(r, double(Ndim));
    }
    double factor = std::pow(phi_target * V_box / sumVol, 1.0 / double(Ndim));
    std::vector<double> scaled = D;
    for (auto &v : scaled) v *= factor;
    return {scaled, factor};
}

double sqDist(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        s += (a[i] - b[i]) * (a[i] - b[i]);
    return s;
}


// ----------------------------------------------------------------------------------
// Binned initialization to reduce collision checks
// ----------------------------------------------------------------------------------
void initialize_positions_binned(
    std::vector<std::vector<double>>& x,
    const std::vector<double>&        D,
    const std::vector<double>&        box,
    std::vector<int8_t> &walls,
    bool fix_height)
{
    size_t N    = x.size();
    size_t Ndim = box.size();
    size_t NB   = 1u << Ndim;  // 2^Ndim bins by midpoint cut per dim
    std::vector<std::vector<size_t>> bins(NB);

    std::vector<double> bin_size(Ndim);
    for (size_t d = 0; d < Ndim; ++d)
        bin_size[d] = box[d] * 0.5;


    auto getBin = [&](const std::vector<double>& pt){
        size_t id = 0;
        for (size_t d = 0; d < Ndim; ++d)
            if (pt[d] > box[d] * 0.5) id |= (1u << d);
        return id;
    };

    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<> ud(0.0,1.0);

    // place first particle in a random bin
    std::vector<double> xt0(Ndim);
    size_t bin_id0 = rng() % NB;
    for (size_t d = 0; d < Ndim; ++d) {
        size_t in_bin = (bin_id0 >> d) & 1;
        double origin = in_bin * bin_size[d];
        xt0[d] = origin + D[0]/2 + ud(rng) * (bin_size[d] - D[0]);
        x[0][d] = xt0[d];
    }
    bins[bin_id0].push_back(0);


    
    size_t count = 1;
    while (count < N) {
        // propose a random position
        std::vector<double> xt(Ndim);
        size_t bin_id = rng() % NB;
        for (size_t d = 0; d < Ndim; ++d) {
            size_t in_bin = (bin_id >> d) & 1;
            double origin = in_bin * bin_size[d];
            xt[d] = origin + D[count]/2 + ud(rng) * (bin_size[d] - D[count]);
        }

        size_t b = getBin(xt);
        bool ok = true;

        // 1) Curved‐boundary check
        if (walls[0] < 0) {
            size_t K = static_cast<size_t>(-walls[0]);  // number of curved dims
            double R = box[0] * 0.5;                    // hypersphere radius
            double sumsq = 0.0;
            for (size_t d = 0; d < K; ++d) {
                double diff = xt[d] - R;
                sumsq += diff * diff;
            }
            if (sumsq > R * R) {
                ok = false;  // outside the curved boundary
            }
        }

        if (ok)
        {
            // only check collisions against particles in the same bin
            for (size_t idx : bins[b]) {
                double r_ij = 0.5 * (D[idx] + D[count]);
                if (sqDist(x[idx], xt) <= r_ij * r_ij) {
                    ok = false;
                    break;
                }
            }

            if (ok) {
                x[count] = std::move(xt);
                bins[b].push_back(count);
                ++count;
            }
        }
    }
}

void initialize_positions_naive(
    std::vector<std::vector<double>>& x,
    const std::vector<double>& D,
    const std::vector<double>& box,
    std::vector<int8_t> &walls,
    bool fix_height)
{
    size_t N = x.size(), Ndim = box.size();
    std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<> ud(0.0,1.0);
    // place first
    for (size_t d = 0; d < Ndim; ++d) x[0][d] = ud(rng) * (box[d]-D[0]) + D[0]/2;
    size_t count = 1;
    while (count < N) {
        std::vector<double> xt(Ndim);
        for (size_t d = 0; d < Ndim; ++d) xt[d] = ud(rng) * (box[d]-D[count]) + D[count]/2;
        bool ok = true;

        // 1) Curved‐boundary check
        if (walls[0] < 0) {
            size_t K = static_cast<size_t>(-walls[0]);  // number of curved dims
            double R = box[0] * 0.5;                    // hypersphere radius
            double sumsq = 0.0;
            for (size_t d = 0; d < K; ++d) {
                double diff = xt[d] - R;
                sumsq += diff * diff;
            }
            if (sumsq > (R-D[count]/2) * (R-D[count]/2)) {
                ok = false;  // outside the curved boundary
            }
        }        

        if (ok)
        {
            for (size_t i = 0; i < count; ++i) {
                double r_ij = 0.5 * (D[i] + D[count]);
                if (sqDist(x[i], xt) <= r_ij*r_ij) { ok = false; break; }
            }
            if (ok) { x[count] = std::move(xt); ++count; }
        }
    }
}

void print_help(const char* prog) {
    std::cout << "Usage: " << prog
              << " --phi <0-1> --N <#particles> --Ndim <dim>"
              << " [--box x,y,...] --dist <type> [--param value ...] --fix_height <true/false>\n\n";
    std::cout << "Flags:\n";
    std::cout << "  --phi        target packing fraction (default=0.05)\n";
    std::cout << "  --N          number of particles\n";
    std::cout << "  --Ndim       dimension (>=2)\n";
    std::cout << "  --box        box dimensions comma-separated (default=1 repeated Ndim)\n";
    std::cout << "  --dist       distribution type: mono, bidisperse, gaussian, biGaussian,\n";
    std::cout << "               lognormal, flat, powerlaw, exponential, weibull, custom\n";
    std::cout << "  --fix_height fix particle heights? true/false (default=false)\n\n";
    std::cout << "Params for distributions:\n";
    std::cout << "  mono:           --d <value>\n";
    std::cout << "  bidisperse:     --d1 <value> --d2 <value> --p <fraction>\n";
    std::cout << "  gaussian:       --mu <value> --sigma <value>\n";
    std::cout << "  biGaussian:     --mu1 <value> --sigma1 <value> --mu2 <value> --sigma2 <value> --p <fraction>\n";
    std::cout << "  lognormal:      --mu <value> --sigma <value>\n";
    std::cout << "                  (mu and sigma of underlying normal, diameters = exp(N(mu, sigma^2)))\n";
    std::cout << "  flat:           --d_min <min> --d_max <max>\n";
    std::cout << "  powerlaw:       --d_min <min> --d_max <max> --exponent <exp>\n";
    std::cout << "  exponential:    --d_min <min> --d_max <max>\n";
    std::cout << "  weibull:        --scale <value> --shape <value>\n";
    std::cout << "  custom:         --custom <file> (file with N diameters)\n";
    std::cout << "\n";
    std::exit(0);
}

std::vector<double> parse_box(const std::string& s) {
    std::vector<double> out;
    size_t p = 0;
    while (p < s.size()) {
        size_t q = s.find(',', p);
        out.push_back(std::stod(s.substr(p, q-p)));
        if (q == std::string::npos) break;
        p = q + 1;
    }
    return out;
}

// trim whitespace from both ends
static inline std::string trim(const std::string &s) {
    auto wsfront = std::find_if_not(s.begin(), s.end(),
                                    [](int c){ return std::isspace(c); });
    auto wsback  = std::find_if_not(s.rbegin(), s.rend(),
                                    [](int c){ return std::isspace(c); }).base();
    return (wsback <= wsfront ? std::string() : std::string(wsfront, wsback));
}

std::vector<std::int8_t>
parse_walls(const std::string &s, std::size_t Ndim, std::vector<double> &box) {
    std::vector<std::int8_t> walls;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ',')) {
        token = trim(token);
        // lowercase for true/false
        std::string lower = token;
        std::transform(lower.begin(), lower.end(), lower.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        if (lower == "true") {
            walls.push_back(1);
        }
        else if (lower == "false") {
            walls.push_back(0);
        }
        else {
            // parse as integer
            try {
                int v = std::stoi(lower);
                if (v < -128 || v > 127) throw std::out_of_range("int8_t");
                walls.push_back(static_cast<std::int8_t>(v));
            }
            catch (...) {
                std::cerr << "Invalid --walls entry: \"" << token << "\"\n";
                std::exit(EXIT_FAILURE);
            }
        }
    }

    // Default fill to zeros if nothing was supplied
    if (walls.empty()) {
        walls.assign(Ndim, 0);
    }

    // Your special initialization when walls[0] < 0
    if (walls[0] < 0) {
        if (walls[0] == -1) {
            walls[0] = -2;
        }
        std::size_t limit = static_cast<std::size_t>(-walls[0]);
        for (std::size_t k = 1; k <= limit && k < Ndim; ++k) {
            walls[k] = -1;
            box[k]   = box[0];
        }
    }

    // If user-specified fewer entries than Ndim, pad the rest with 0
    if (walls.size() < Ndim) {
        walls.resize(Ndim, 0);
    }

    return walls;
}

inline double sphereVolume(double r, size_t Ndim) {
    return std::pow(M_PI, Ndim/2.0) * std::pow(r, Ndim) /
           std::tgamma(Ndim/2.0 + 1.0);
}

int main(int argc, char** argv) {
    double phi = 0.05;
    size_t N = 0, Ndim = 0;
    std::vector<double> box;
    std::vector<std::int8_t> walls;
    bool fix_height=false;
    Distribution dist;

    if (argc < 2) { print_help(argv[0]); return 1; }
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--help") { print_help(argv[0]); return 0; }
        else if (a == "--phi"  && i+1<argc) phi    = std::stod(argv[++i]);
        else if (a == "--N"    && i+1<argc) N      = std::stoul(argv[++i]);
        else if (a == "--Ndim" && i+1<argc) Ndim   = std::stoul(argv[++i]);
        else if (a == "--box"  && i+1<argc) box    = parse_box(argv[++i]);
        else if (a == "--walls"  && i+1<argc) walls = parse_walls(argv[++i], Ndim, box);
        else if (a == "--dist" && i+1<argc) dist.type = argv[++i];
        else if (a == "--fix_height") {fix_height = true;}
        else if (a == "--d"    && i+1<argc) dist.d    = std::stod(argv[++i]);
        else if (a == "--d1"   && i+1<argc) dist.d1   = std::stod(argv[++i]);
        else if (a == "--d2"   && i+1<argc) dist.d2   = std::stod(argv[++i]);
        else if (a == "--p"    && i+1<argc) dist.p    = std::stod(argv[++i]);
        else if (a == "--mu"   && i+1<argc) dist.mu   = std::stod(argv[++i]);
        else if (a == "--sigma"&& i+1<argc) dist.sigma= std::stod(argv[++i]);
        else if (a == "--mu1"  && i+1<argc) dist.mu1  = std::stod(argv[++i]);
        else if (a == "--sigma1"&&i+1<argc) dist.sigma1=std::stod(argv[++i]);
        else if (a == "--mu2"  && i+1<argc) dist.mu2  = std::stod(argv[++i]);
        else if (a == "--sigma2"&&i+1<argc) dist.sigma2=std::stod(argv[++i]);
        else if (a == "--d_min"&& i+1<argc) dist.d_min=std::stod(argv[++i]);
        else if (a == "--d_max"&& i+1<argc) dist.d_max=std::stod(argv[++i]);
        else if (a == "--exponent"&&i+1<argc) dist.exponent=std::stod(argv[++i]);
        else {
            std::cerr << "Unknown flag: " << a << "\n";
            print_help(argv[0]);
            return 1;
        }
    }

    if (Ndim < 2) { std::cerr << "--Ndim must be >=2\n"; return 1; }
    if (box.empty()) box.assign(Ndim, 1.0);
    if (phi <= 0 || phi >= 1 || N == 0) {
        std::cerr << "Missing or invalid phi/N\n";
        print_help(argv[0]);
        return 1;
    }

    double Llast = box[0];
    double h = 1;
    if (fix_height)
    {
        h = box[Ndim-1];
        box[Ndim-1] = box[0];
    }

    double phi_modifier = 1.0;
    if (walls[0] < 0);{phi_modifier = sphereVolume(1/2., -walls[0]);}

    std::mt19937 rng{std::random_device{}()};
    auto rawD = generate_diameter_distribution(N, dist, rng);
    auto [D_scaled, factor] = scale_diametersND(rawD, phi*phi_modifier, box, Ndim);

    // sort indices by descending diameter
    std::vector<size_t> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
        [&](size_t a, size_t b){ return D_scaled[a] > D_scaled[b]; });

    // reorder D_scaled into descending order for initialization
    std::vector<double> D_sorted(N);
    for (size_t i = 0; i < N; ++i)
        D_sorted[i] = D_scaled[idx[i]];
    D_scaled = std::move(D_sorted);

    //Rescale diameters and box size to match targeted fixed height
    if (fix_height)
    {

        // 1) scale for diameters
        double s = std::pow(h * D_scaled[0] / Llast, 1.0 / (Ndim-1));

        // 2) scale for the last box side
        double eps = std::pow(s, Ndim);

        // apply
        for (size_t i = 0; i < D_scaled.size(); ++i) {
            D_scaled[i] *= s; 
        }

        box[Ndim-1] *= eps;
    }

    // --- use the binned initializer here ---
    std::vector<std::vector<double>> x(N, std::vector<double>(Ndim));
    initialize_positions_binned(x, D_scaled, box, walls, fix_height);
    // initialize_positions_naive(x, D_scaled, box, walls, fix_height);

    // Reorder x and D_scaled back to original
    // std::vector<std::vector<double>> x_original(N, std::vector<double>(Ndim));
    // std::vector<double> D_original(N);
    // for (size_t i = 0; i < N; ++i) {
    //     x_original[idx[i]] = x[i];
    //     D_original[idx[i]] = D_scaled[i];
    // }
    // x = std::move(x_original);
    // D_scaled = std::move(D_original);

    for (size_t i = 0; i < N; ++i) {
        for (size_t d = 0; d < Ndim; ++d)
        {
            std::cout << x[i][d] << ' ';
        }
        std::cout << D_scaled[i] << '\n';
    }

    return 0;
}
