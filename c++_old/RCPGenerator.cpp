// CreatePacking.cpp — Exact MATLAB-to-C++ Translation using 2D arrays
// ----------------------------------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdint>
#include <iomanip>

//------------------------------------------------------------------------------
// System parameters (from MATLAB Initialization)
//------------------------------------------------------------------------------
static constexpr double ALPHA_MAX     = 0.0025;        // Adam learning rate
static constexpr double BETA1        = 0.9;           // Exponential decay rate for first moment
static constexpr double BETA2        = 0.999;         // Exponential decay rate for second moment
static constexpr double EPSILON      = 1e-8;          // Small value to prevent division by zero
static constexpr size_t N_STEPS      = 150000;         // Max optimization steps
static constexpr double DT           = 0.1;           // Time step for Verlet
static const std::string METHOD      = "ADAM";       // Optimization method
static const std::vector<uint32_t> MAX_NEIGHBORS = {300, 750, 5500, 5500};
static constexpr double DELTA_PHI0   = 1.5 * 1.5e-3;   // Initial phi step size

// Global box dimensions for periodic boundaries
static std::vector<double> g_box;

//------------------------------------------------------------------------------
// Function prototypes (match MATLAB signatures)
//------------------------------------------------------------------------------
void parseArgs(int argc, char** argv,
               std::string &filePath,
               std::string &outPath,
               bool &hasOutput,
               std::vector<double> &box,
               std::vector<int8_t> &walls,
               uint32_t &neighborMax,
               uint32_t &seed,
               bool &verbose,
               bool &fix_height,
               uint32_t &save_interval);

void loadData(const std::string &filePath,
              size_t &N,
              size_t &Ndim,
              std::vector<std::vector<double>> &x,
              std::vector<double> &D);

std::string stripTxt(std::string filename);

void writePositions(
    std::string outPath,
    const std::vector<std::vector<double>> &x,
    const std::vector<double> &D
); 

void printStatus(
    size_t step,
    double phi,
    double dphi,
    double U,
    double max_min_dist,
    double F_mag,
    size_t max_neighbors,
    size_t num_changes,
    double alpha,
    double max_delta_x
);

double delta_x(
    const std::vector<std::vector<double>> &x,
    const std::vector<std::vector<double>> &x_old,
    const std::vector<double>              &D);

double norm(const std::vector<double> &a,
    const std::vector<double> &b);

double computeMeanForce(
    const std::vector<std::vector<double>> &F,
    double                                  Lc,
    const std::vector<size_t>             &z,
    size_t                                  Ndim);

double mean(const std::vector<double> &v,
    size_t start,
    size_t end);


double maxElem(const std::vector<std::vector<double>> &M);
double computeMaxMinDist(const std::vector<std::vector<double>> &min_dist);
size_t computeMaxNeighbors(const std::vector<std::vector<uint32_t>> &pairs);
size_t computeNumChanges(const std::vector<uint32_t> &refresh);

void computeMinMax(const std::vector<std::vector<double>> &x);

std::vector<size_t> sortIndicesByColumn(
    const std::vector<std::vector<double>> &x,
    size_t col);

inline double sphereVolume(double r, size_t Ndim);

std::pair<std::vector<double>, double> scale_diametersND(
    const std::vector<double> &D,
    double phi_target,
    const std::vector<double> &box,
    size_t Ndim,
    bool fix_height);


void GetPairsND_3(
    size_t N,
    size_t Ndim,
    const std::vector<std::vector<double>> &x,
    const std::vector<double> &D,
    const std::vector<double> &box,
    const std::vector<int8_t> &walls,
    const std::vector<uint32_t> &refresh,
    size_t neighborMax,                           // <-- new parameter
    std::vector<std::vector<uint32_t>> &pairs);


void GetForcesND_3(
    const std::vector<std::vector<uint32_t>> &pairs,
    const std::vector<std::vector<double>>   &x,
    const std::vector<double>                &D,
    const std::vector<double>                &box,
    const std::vector<int8_t>                &walls,
    std::vector<std::vector<double>>         &F,
    double                                   &U,
    std::vector<std::vector<double>>         &min_dist,
    double                                   &max_min_dist, 
    double                                   &Lc,
    std::vector<size_t>                      &z);

void AdamUpdate(
    const std::string &method,
    size_t N,
    size_t Ndim,
    const std::vector<std::vector<double>> &F,
    double beta1,
    double beta2,
    size_t t,
    double dt,
    double verlet_drag,
    std::vector<std::vector<double>> &m,
    std::vector<std::vector<double>> &v,
    std::vector<std::vector<double>> &v_max,
    std::vector<std::vector<double>> &v_update,
    std::vector<std::vector<double>> &m_hat,
    std::vector<std::vector<double>> &v_hat,
    std::vector<std::vector<double>> &a,
    std::vector<std::vector<double>> &v_verlet,
    std::vector<std::vector<double>> &a_old);

//------------------------------------------------------------------------------
int main(int argc, char** argv) {
    
    // ###################################################
    // ---- Define Parameters and Configure System -----
    // ###################################################

    // 1) Parse command-line arguments
    std::string filePath, outPath;
    std::vector<double> box;
    std::vector<int8_t> walls;
    uint32_t neighborMax = 0, seed = 0;
    bool verbose = false;
    bool hasOutput = false;
    bool fix_height = false;
    uint32_t save_interval = 0;
    parseArgs(argc, argv, filePath, outPath, hasOutput, box, walls, neighborMax, seed, verbose, fix_height, save_interval);
    if (hasOutput == false){save_interval=0;}

    // 2) Load positions x[N][Ndim] and diameters D[N]
    size_t N = 0, Ndim = 0;
    std::vector<std::vector<double>> x;
    std::vector<double> D;
    loadData(filePath, N, Ndim, x, D);
    std::vector<std::vector<double>> x_last = x;

    // 3) Initialize box
    if (box.empty()) box.assign(Ndim, 1.0);
    else if (box.size() != Ndim) {
        std::cerr << "Error: --box length ("<<box.size()<<") != Ndim("<<Ndim<<")"<<std::endl;
        return EXIT_FAILURE;
    }
    if (walls.empty()) {
        walls.assign(Ndim, 0);
    } 
    else if (walls.size() != Ndim) {
        std::cerr << "Error: --walls expects " << Ndim << " entries, got " << walls.size() << "\n";
        return EXIT_FAILURE;
    }

    g_box = box;

    // 4) Initial phi and scale diameters
    double phi0 = 0.025;
    double phi = phi0;
    double phi_modifier = 1;
    if (walls[0] < 0)
    {
        if (walls[0] == -1){walls[0] = -2;}
        for (size_t k = 1; k <= -walls[0]; ++k) {walls[k] = -1;box[k]=box[0];}

        phi_modifier = sphereVolume(1/2., -walls[0]);// provide correction to container volume now its a hypersphere and not a box
    }

    if (fix_height)
    {
        box[Ndim-1] = box[Ndim-1]*D[0];
    }

    if      (Ndim == 2) phi = 0.575;
    else if (Ndim == 3) phi = 0.33;
    else if (Ndim == 4) phi = 0.15;
    else if (Ndim == 5) phi = 0.10;
    else if (Ndim > 5)  phi = 0.10 / std::sqrt(double(Ndim));
    double delta_phi0 = DELTA_PHI0;
    double delta_phi = DELTA_PHI0;
    double dphi = 0.0;
    size_t count = 0;
    bool update_flag = false;
    int direction_flag = 0;
    double U_threshold = (2.5E-4)*(2.5E-4); //% Energy threshold for minimal overlap
    double F_tol = std::sqrt(U_threshold)/50; //% Force threshold for concluding phi too large and needs to be lowered    
    // D = scale_diametersND(D, phi, box, Ndim);
    auto [D_scaled, factor] = scale_diametersND(D, phi*phi_modifier, box, Ndim, fix_height);
    if (fix_height){
        box[Ndim-1] = box[Ndim-1]*factor;
        for (size_t i = 0; i<N; ++i)
        {
            x[i][Ndim-1] = x[i][Ndim-1]*factor;
        }
    }
    D = D_scaled;
    double Dmin=0;
    auto it = std::min_element(D.begin(), D.end());
    Dmin = *std::min_element(D.begin(), D.end());
    double Dmin_last = Dmin;
    double max_min_dist = 0;
    if (verbose) std::cout<<"Initial phi0="<<phi0<<std::endl;

    // 5) Default NeighborMax if unset, based on Ndim
    if (neighborMax == 0) {
        size_t idx = (Ndim >= 2 ? std::min<size_t>(Ndim-2, MAX_NEIGHBORS.size()-1) : 0);
        neighborMax = MAX_NEIGHBORS[idx];
        if (verbose) std::cout << "[info] --NeighborMax default=" << neighborMax << std::endl;
    }
    
    // 6) Prepare refresh flags and pairs matrix
    std::vector<uint32_t> refresh(N, 1);
    std::vector<std::vector<uint32_t>> pairs(
        N, std::vector<uint32_t>(neighborMax+1, 0)
    );
    
    // 7) Allocate force & contact outputs
    std::vector<std::vector<double>> F(N, std::vector<double>(Ndim, 0.0));
    std::vector<std::vector<double>> min_dist(
        N, std::vector<double>(neighborMax, 0.0)
    );
    std::vector<size_t> z(N, 0);
    double U = 0.0, Lc = 0.0;

    // 8) Initialize Adam/Verlet state arrays
    std::string method=METHOD;
    double beta1=BETA1, beta2=BETA2, t = 0, R = 1.00005, alpha = ALPHA_MAX, alpha_max = ALPHA_MAX, N_steps = N_STEPS, epsilon = EPSILON;
    double dt = DT; 
    double verlet_drag = 2;
    std::vector<std::vector<double>>
        m(N, std::vector<double>(Ndim,0.0)),
        v(N, std::vector<double>(Ndim,0.0)),
        v_max(N, std::vector<double>(Ndim,0.0)),
        v_update(N, std::vector<double>(Ndim,0.0)),
        m_hat(N, std::vector<double>(Ndim,0.0)),
        v_hat(N, std::vector<double>(Ndim,0.0)),
        a(N, std::vector<double>(Ndim,0.0)),
        v_verlet(N, std::vector<double>(Ndim,0.0)),
        a_old(N, std::vector<double>(Ndim,0.0));

    size_t LastPhiUpdate = 0;
    size_t LastAlphaUpdate = 0;
    // Allocate and zero‐initialize U_history to length N_steps
    std::vector<double> U_history(N_steps, 0.0);

    // Similarly, if you need phi_history and F_history:
    std::vector<double> phi_history(N_steps, 0.0);
    std::vector<double> F_history(N_steps,   0.0);

    double phi_min = 0.8;
    if (Ndim == 2){ phi_min = 0.80; }
    if (Ndim == 3){ phi_min = 0.61; }
    if (Ndim == 4){ phi_min = 0.43; }
    if (Ndim == 5){ phi_min = 0.28; }
    if (Ndim > 5){ phi_min = 0.28/(Ndim-4); }
    
    // ----

    // ###################################################
    // ---- Optimize Packing ----
    // ###################################################

    // Build and set parameters
    GetPairsND_3(N, Ndim, x, D, box, walls, refresh, neighborMax, pairs);

    GetForcesND_3(pairs, x, D, box, walls, F, U, min_dist, max_min_dist, Lc, z); 
    
    // AdamStep(method, N, Ndim, F, beta1, beta2, t, dt, verlet_drag, m, v, v_max, v_update, m_hat, v_hat,a, v_verlet, a_old);

    size_t steps = 0;
    for (size_t step = 1; step <= N_steps; ++step) {
        // --- Set Optimizer Given Current Conditions ---
        steps = step;
        if (dphi != 0) {
            if (method != "ADAM") {
                method = "ADAM";
                LastPhiUpdate = step;
            }
            if (dphi > 0) {
                alpha = std::min(alpha * 1.25, alpha_max);
            }
        }
        if (step - LastPhiUpdate > 2500 && method == "ADAM") {
            method = "AMSGrad";
        }
        if (step - LastPhiUpdate > 4000 && method == "AMSGrad") {
            method = "Verlet";
            for (size_t k = 0; k < N; ++k) {
                std::fill(v_verlet[k].begin(), v_verlet[k].end(), 0.0);
                std::fill(a[k].begin(), a[k].end(), 0.0);
                std::fill(a_old[k].begin(), a_old[k].end(), 0.0);
            }
        }
    
        // --- Learning Rate Management ---
        t += 1;
        alpha *= R;
        if (dphi > 0) {
            LastAlphaUpdate = step;
        }
        if (step - LastAlphaUpdate > 500 && method == "ADAM") {
            double trend1 = mean(U_history, step-450, step-250);
            double trend2 = mean(U_history, step-75,  step-1);
            if (trend1 / trend2 < 0.85) {
                LastAlphaUpdate = step;
                alpha *= 0.5;
            }
            else if (step - LastAlphaUpdate > 4000) {
                double trend3 = mean(U_history, step-200, step-1);
                double trend4 = mean(U_history, step-3000, step-2500);
                if (trend3 / trend4 > 0.5) {
                    LastAlphaUpdate = step;
                    alpha *= 1.5;
                }
            }
        }
        if (step % 15000 == 0) {
            alpha = std::min(alpha_max, 100*alpha);
        }
        if (alpha < 0.5e-5) {
            alpha = alpha_max/10;
            for (size_t k = 0; k < N; ++k)
                std::fill(m[k].begin(), m[k].end(), 0.0),
                std::fill(v[k].begin(), v[k].end(), 0.0);
            t = 1;
        }
    
        // --- Updating Phi ---
        if (step == 150) {
            delta_phi = delta_phi0;
        }
        phi = std::max(phi + dphi, phi_min);
        // D_new = scale_diametersND(D, phi, box);
        auto [D_new, scale] = scale_diametersND(D, phi*phi_modifier, box, Ndim, fix_height);
        if (fix_height){
            box[Ndim-1] = box[Ndim-1]*scale;
            for (size_t i = 0; i<N; ++i)
            {
                x[i][Ndim-1] = x[i][Ndim-1]*scale;
            }
        }
        D = D_new;
        Dmin *= scale;
        phi_history[step-1] = phi;
    
        // --- Managing Pairs List ---
        if (step % 750 == 0 || Dmin_last / Dmin > 1.05) {
            std::fill(refresh.begin(), refresh.end(), 1u);
            GetPairsND_3(N, Ndim, x, D, box, walls, refresh, neighborMax, pairs);
            x_last = x;
        } else {
            std::fill(refresh.begin(), refresh.end(), 0u);
            for (size_t k = 0; k < N; ++k) {
                if (norm(x[k], x_last[k]) > Dmin/4) {
                    refresh[k] = 1;
                    x_last[k] = x[k];
                    update_flag = true;
                }
            }
            if (update_flag) {
                GetPairsND_3(N, Ndim, x, D, box, walls, refresh, neighborMax, pairs);
            }
        }
        Dmin_last = Dmin;
        update_flag = false;
    
        // --- Forces ---
        GetForcesND_3(pairs, x, D, box, walls, F, U, min_dist, max_min_dist, Lc, z);
        double F_magnitude = computeMeanForce(F, Lc, z, Ndim);
        F_history[step-1] = F_magnitude;
 

        // --- Optimizer updates ---
        AdamUpdate(method, N, Ndim,
                 F, beta1, beta2, t, dt, verlet_drag,
                 m, v, v_max, v_update,
                 m_hat, v_hat,
                 a, v_verlet, a_old);
    
        // --- Termination and phi updating conditions ---
        if (U < U_threshold || max_min_dist < std::sqrt(U_threshold)*10) {
            dphi = delta_phi; //* (1 + (randUniform()-0.5)/10);
            if ( verbose){
                std::cout << "Increasing phi:" << std::endl;
                printStatus(step, phi, delta_phi, U, max_min_dist, F_magnitude, computeMaxNeighbors(pairs), computeNumChanges(refresh), alpha, 0);
            }
            direction_flag = 1;
            ++count;
        } else if (U > U_threshold && step > 125 && phi > phi_min) {
            if (F_magnitude < F_tol) {
                if (direction_flag == 1)
                    delta_phi *= 0.5;
                dphi = -delta_phi/2; // * (1 + (randUniform()-0.5)/10);
                direction_flag = -1;
            } else {
                dphi = 0;
            }
            count = 0;
        }
        if (count > 10 && Ndim == 2) {
            delta_phi = std::min(delta_phi * 1.5, delta_phi0);
            count = 0;
        }
    
        // --- Print / Terminate if converged ---
        if (delta_phi < 5e-6 &&
           (U < U_threshold || max_min_dist < std::sqrt(U_threshold)*10)) {
            break;
        }
    
        // --- Positions updates ---
        auto x_old = x;
        if (method != "Verlet") {
            for (size_t k = 0; k < N; ++k)
                for (size_t d = 0; d < Ndim; ++d)
                    x[k][d] -= alpha * m_hat[k][d] / (std::sqrt(v_hat[k][d]) + epsilon);
        } else {
            for (size_t k = 0; k < N; ++k)
                for (size_t d = 0; d < Ndim; ++d)
                    x[k][d] += v_verlet[k][d]*dt + 0.5*a_old[k][d]*dt*dt;
        }
        for (size_t k = 0; k < N; ++k)
            for (size_t d = 0; d < Ndim; ++d)
                x[k][d] = std::fmod(x[k][d] + box[d], box[d]);
    
        U_history[step-1] = U;
    
        if (dphi == 0 && delta_x(x, x_old, D) < std::sqrt(U_threshold)/2500) {
            dphi = -delta_phi/2; //* (1 + (randUniform()-0.5)/10);
            if (delta_phi < 5e-6 &&
               (U < U_threshold || max_min_dist < std::sqrt(U_threshold)*10)) {
                break;
            }
        }
    
        if ((step < 100 || step % 1000 == 0) && verbose) {
            printStatus(step, phi, delta_phi, U, max_min_dist, F_magnitude, computeMaxNeighbors(pairs), computeNumChanges(refresh), alpha, delta_x(x, x_old, D));
        }

        // every save_interval steps, dump out a file
        if (save_interval > 0){
          if ((step-1) % save_interval == 0) {
              std::ostringstream oss;
              oss << outPath
                  << "_"
                  << std::setw(6)        // pad to width=6
                  << std::setfill('0')
                  << step-1
                  << ".txt";            // add extension back
  
              std::string filename = oss.str(); 
              writePositions(filename, x, D);
          }
        }
        
    }

    if (hasOutput) {
        
        if (steps == N_steps){
            std::cout << "Could not relax packing fully ... Reached max " << steps << " steps" << std::endl;
        }else{
            std::cout << "Done! Took " << steps << " steps" << std::endl;
        }
        std::cout << "Final Packing Fraction " << phi << std::endl;
        std::cout << "Max overlap: " << max_min_dist << std::endl;

        writePositions(outPath, x, D);
        std::cout << "Output written to " << outPath << ".txt" << std::endl;
    }
    
    else {
        // no --output: just dump positions & radii to stdout
        for (size_t i = 0; i < N; ++i) {
            for (size_t d = 0; d < Ndim; ++d)
                std::cout << x[i][d] << ' ';
            std::cout << D[i] << '\n';
        }
    }

    return 0;
}

//------------------------------------------------------------------------------
// parseArgs implementation
//------------------------------------------------------------------------------
void parseArgs(int argc, char** argv,
               std::string &filePath,
               std::string &outPath,
               bool &hasOutput,
               std::vector<double> &box,
               std::vector<int8_t> &walls,
               uint32_t &neighborMax,
               uint32_t &seed,
               bool &verbose,
               bool &fix_height,
               uint32_t &save_interval) {
    outPath = "packing_out.txt";
    if (argc < 2) {
        std::cerr<<"Usage: CreatePacking.exe --file input.txt [--output out.txt] [--box v1,..] [--walls true/false,..] [--fix-height] [--save-interval] [--NeighborMax N] [--seed S] [--verbose]"<<std::endl;
        std::exit(EXIT_FAILURE);
    }
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--file" && i+1 < argc) filePath = argv[++i];
        else if (arg == "--output" && i+1 < argc) {
            outPath = argv[++i];
            outPath = stripTxt(outPath);
            hasOutput  = true;
        } else if (arg == "--box" && i+1 < argc) {
            box.clear(); double val; std::stringstream ss(argv[++i]);
            while (ss>>val) { box.push_back(val); if (ss.peek()==',') ss.ignore(); }
        } else if (arg == "--NeighborMax" && i+1 < argc) {
            neighborMax = std::stoul(argv[++i]);
        } else if (arg == "--seed" && i+1 < argc) {
            seed = std::stoul(argv[++i]);
        } else if (arg == "--verbose") {
            verbose = true;
        } else if (arg == "--fix-height") {
            fix_height = true;
        }else if (arg == "--save-interval") {
            save_interval = std::stoul(argv[++i]);
        }else if (arg == "--walls" && i+1 < argc) {
            // parse comma‐separated 0/1 (or true/false) flags
            walls.clear();
            std::string token;
            std::stringstream ss(argv[++i]);
            while (std::getline(ss, token, ',')) {
                if (token == "1" || token == "true"){
                    walls.push_back(1);}
                else if (token == "0" || token == "false"){
                    walls.push_back(0);}
                else if (token == "0" || token == "false"){
                    walls.push_back(0);}
                else if (token == "-1"){
                    walls.push_back(-1);}
                else if (token == "-2"){
                    walls.push_back(-2);}
                else if (token == "-3"){
                    walls.push_back(-3);}
                else if (token == "-4"){
                    walls.push_back(-4);}
                else if (token == "-5"){
                    walls.push_back(-5);}
                else if (token == "-6"){
                    walls.push_back(-6);}
                else if (token == "-7"){
                    walls.push_back(-7);}
                else if (token == "-8"){
                    walls.push_back(-8);}
                else if (token == "-9"){
                    walls.push_back(-9);}
                else {
                    std::cerr << "Invalid --walls entry: " << token << "\n";
                    std::exit(EXIT_FAILURE);
                }
            }
        } 
        else {
            std::cerr<<"Unknown flag: "<<arg<<std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
   
}

//------------------------------------------------------------------------------
// loadData implementation
//------------------------------------------------------------------------------
void loadData(const std::string &filePath,
    size_t &N,
    size_t &Ndim,
    std::vector<std::vector<double>> &x,
    std::vector<double> &D)
{
    // 1) Read every non‐empty line into a buffer
    std::vector<std::string> lines;
    {
        std::istream* in = &std::cin;
        std::ifstream fin;
        if (!filePath.empty()) {
            fin.open(filePath);
            if (!fin) {
                std::cerr << "Cannot open " << filePath << "\n";
                std::exit(EXIT_FAILURE);
            }
            in = &fin;
        }
        std::string line;
        while (std::getline(*in, line)) {
            if (line.find_first_not_of(" \t\r\n") != std::string::npos)
                lines.push_back(line);
        }
    }

    // 2) Figure out N and Ndim from the buffered lines
    N = lines.size();
    if (N == 0) {
        std::cerr << "No data lines found\n";
        std::exit(EXIT_FAILURE);
    }
    {
        std::stringstream ss(lines[0]);
        size_t cols = 0;
        double tmp;
        while (ss >> tmp) ++cols;
        if (cols < 2) {
            std::cerr << "Expect at least 2 columns per line\n";
            std::exit(EXIT_FAILURE);
        }
        Ndim = cols - 1;
    }

    // 3) Allocate storage
    x.assign(N, std::vector<double>(Ndim));
    D.assign(N, 0.0);

    // 4) Parse again from the buffered lines
    for (size_t i = 0; i < N; ++i) {
        std::stringstream ss(lines[i]);
        for (size_t d = 0; d < Ndim; ++d)
        ss >> x[i][d];
        ss >> D[i];
    }
}

//------------------------------------------------------------------------------
// Strip .txt
//------------------------------------------------------------------------------

std::string stripTxt(std::string filename) {
    const std::string ext = ".txt";
    if (filename.size() >= ext.size() &&
        filename.substr(filename.size() - ext.size()) == ext)
    {
        filename.erase(filename.size() - ext.size());
    }
    return filename;
}

//------------------------------------------------------------------------------
// Save Data
//------------------------------------------------------------------------------
void writePositions(
    std::string outPath,
    const std::vector<std::vector<double>> &x,
    const std::vector<double> &D
) {

    if (outPath.size() < 4 || outPath.substr(outPath.size()-4) != ".txt") {
        outPath += ".txt";
    }
    
    const size_t N = x.size();
    if (D.size() != N) {
        std::cerr << "writePositions error: x.size()=" << N
                  << " but D.size()=" << D.size() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    const size_t Ndim = x.empty() ? 0 : x[0].size();

    std::ofstream fout(outPath);
    if (!fout) {
        std::cerr << "Cannot open output file: " << outPath << std::endl;
        std::exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < N; ++i) {
        if (x[i].size() != Ndim) {
            std::cerr << "writePositions error: inconsistent dimension at row "
                      << i << std::endl;
            std::exit(EXIT_FAILURE);
        }
        // write all coordinates
        for (size_t d = 0; d < Ndim; ++d) {
            fout << x[i][d] << ' ';
        }
        // write diameter and newline
        fout << D[i] << '\n';
    }
    fout.close();

    //std::cout << "Output written to " << outPath << std::endl;
}

//------------------------------------------------------------------------------
// Print Status implementation
//------------------------------------------------------------------------------
void printStatus(
    size_t step,
    double phi,
    double dphi,
    double U,
    double max_min_dist,
    double F_mag,
    size_t max_neighbors,
    size_t num_changes,
    double alpha,
    double max_delta_x
) {
    std::printf(
        "Step %zu: Phi = %.4f, dPhi = %0.3e, U = %.2e, min_dist = %.2e, "
        "|F| = %.2e, Neighbors = %zu, N changes = %zu, alpha = %.2e, max d = %.2e\n",
        step,
        phi,
        dphi,
        U,
        max_min_dist,
        F_mag,
        max_neighbors,
        num_changes,
        alpha,
        max_delta_x
    );
}

double delta_x(
    const std::vector<std::vector<double>> &x,
    const std::vector<std::vector<double>> &x_old,
    const std::vector<double>              &D)
{
    size_t N = x.size();
    if (x_old.size() != N || D.size() != N) {
        throw std::invalid_argument("delta_x: mismatched sizes");
    }

    double sum = 0.0;
    for (size_t i = 0; i < N; ++i) {
        // compute Euclidean distance between x_old[i] and x[i]
        double d2 = 0.0;
        for (size_t d = 0; d < x[i].size(); ++d) {
            double diff = x_old[i][d] - x[i][d];
            d2 += diff * diff;
        }
        double dist = std::sqrt(d2);
        sum += dist / D[i];
    }

    return sum / double(N);
}

// Return ||a – b||₂ for two N‑dim vectors
double norm(const std::vector<double> &a,
    const std::vector<double> &b)
{
if (a.size() != b.size()) {
throw std::invalid_argument("norm: vector sizes differ");
}
double sum2 = 0.0;
for (size_t i = 0; i < a.size(); ++i) {
double d = a[i] - b[i];
sum2 += d*d;
}
return std::sqrt(sum2);
}

double mean(const std::vector<double> &v,
    size_t start,
    size_t end) {
double sum = 0.0;
for (size_t k = start; k <= end; ++k) {
sum += v[k];
}
return sum / double(end - start + 1);
}

double computeMeanForce(
    const std::vector<std::vector<double>> &F,
    double                                  Lc,
    const std::vector<size_t>             &z,
    size_t                                  Ndim)
{
    size_t N = F.size();
    double sum_mag = 0.0;
    double sum_z   = 0.0;

    for (size_t i = 0; i < N; ++i) {
        // compute ||F[i]||₂
        double d2 = 0.0;
        for (size_t d = 0; d < Ndim; ++d) {
            d2 += F[i][d] * F[i][d];
        }
        sum_mag += std::sqrt(d2);
        sum_z   += double(z[i]);
    }

    double mean_mag = sum_mag / double(N);
    double mean_z   = sum_z   / double(N);

    return mean_mag / Lc / mean_z * double(Ndim) / std::sqrt(double(Ndim));
}

double maxElem(const std::vector<std::vector<double>> &M) {
    if (M.empty()) {
        throw std::invalid_argument("maxElem: empty matrix");
    }
    double m = M[0][0];
    for (const auto &row : M) {
        if (row.empty()) continue;
        double row_max = *std::max_element(row.begin(), row.end());
        if (row_max > m) m = row_max;
    }
    return m;
}

double computeMaxMinDist(const std::vector<std::vector<double>> &min_dist) {
    double m = 0.0;
    for (auto &row : min_dist)
        for (double v : row)
            if (v > m) m = v;
    return m;
}

size_t computeMaxNeighbors(const std::vector<std::vector<uint32_t>> &pairs) {
    size_t m = 0;
    for (auto &row : pairs)
        if (row[0] > m) m = row[0];
    return m;
}

size_t computeNumChanges(const std::vector<uint32_t> &refresh) {
    size_t sum = 0;
    for (auto v : refresh)
        sum += v;
    return sum;
}

//------------------------------------------------------------------------------
// computeMinMax implementation
//------------------------------------------------------------------------------
void computeMinMax(const std::vector<std::vector<double>> &x) {
    size_t N = x.size();
    size_t Ndim = x.empty()?0:x[0].size();
    for (size_t d=0; d<Ndim; ++d) {
        double mn=x[0][d], mx=x[0][d];
        for (size_t i=1; i<N; ++i) { mn=std::min(mn,x[i][d]); mx=std::max(mx,x[i][d]); }
        std::cout<<"Column "<<d<<": ["<<mn<<", "<<mx<<"]"<<std::endl;
    }
}

//------------------------------------------------------------------------------
// sortIndicesByColumn implementation
//------------------------------------------------------------------------------
std::vector<size_t> sortIndicesByColumn(
    const std::vector<std::vector<double>> &x,
    size_t col)
{
    size_t N = x.size();
    std::vector<size_t> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&](size_t a, size_t b){ return x[a][col] < x[b][col]; });
    return idx;
}


//------------------------------------------------------------------------------
// sphereVolume implementation
//------------------------------------------------------------------------------
inline double sphereVolume(double r, size_t Ndim) {
    return std::pow(M_PI, Ndim/2.0) * std::pow(r, Ndim) /
           std::tgamma(Ndim/2.0 + 1.0);
}

//------------------------------------------------------------------------------
// scale_diametersND implementation
//------------------------------------------------------------------------------
#include <vector>
#include <cmath>
#include <utility>

// Assumes you have a sphereVolume(radius, Ndim) function defined elsewhere.

/// Scales the diameters D so that the resulting packing fraction becomes phi_target.
/// Returns a pair: (scaled_D_vector, scale_factor).
std::pair<std::vector<double>, double> scale_diametersND(
    const std::vector<double> &D,
    double phi_target,
    const std::vector<double> &box,
    size_t Ndim,
    bool fix_height)
{
    // Compute box volume
    double V_box = 1.0;
    for (double L : box) {
        V_box *= L;
    }

    // Compute total volume of spheres given current diameters
    double sumVol = 0.0;
    for (double Di : D) {
        double radius = Di * 0.5;
        sumVol += sphereVolume(radius, Ndim);
    }

    // Determine scale factor so that phi_target = sumVol_scaled / V_box
    //
    // sumVol_scaled = sumVol * factor^Ndim
    // => factor = (phi_target * V_box / sumVol)^(1/Ndim)
    double factor = 1.0;
    if (fix_height)
    {
        factor = std::pow(phi_target * V_box / sumVol, 1.0 / double(Ndim-1));
    }
    else
    {
        factor = std::pow(phi_target * V_box / sumVol, 1.0 / double(Ndim));
    }

    // Apply factor to each diameter
    std::vector<double> D_scaled = D;
    for (double &Di : D_scaled) {
        Di *= factor;
    }

    return { D_scaled, factor };
}

void GetPairsND_3(
    size_t N,
    size_t Ndim,
    const std::vector<std::vector<double>> &x,
    const std::vector<double> &D,
    const std::vector<double> &box,
    const std::vector<int8_t> &walls,
    const std::vector<uint32_t> &refresh,
    size_t neighborMax,
    std::vector<std::vector<uint32_t>> &pairs)
{
    double Dmax = *std::max_element(D.begin(), D.end());
    double Dmin = *std::min_element(D.begin(), D.end());
    std::vector<double> Dsorted = D;
    std::sort(Dsorted.begin(), Dsorted.end());
    double median = Dsorted[N/2];
    double t = std::min(std::max(median*1.05, Dmin*3.0), Dmax);

    auto sort_idx = sortIndicesByColumn(x, 0);
    std::vector<size_t> sort_loc(N);
    for (size_t k = 0; k < N; ++k)
        sort_loc[sort_idx[k]] = k;

    bool warningIssued = false;

    for (size_t i = 0; i < N; ++i) {
        if (refresh[i] != 1) continue;
        pairs[i][0] = 0;

        double r_c_max = (D[i] + 1.01 * Dmax) / 2 + t;
        for (int dir : {-1, 1}) {
            size_t jdx = 0;
            bool go = true;
            while (go) {
                ++jdx;
                if (jdx > N/2) { go = false; break; }

                size_t pos = (sort_loc[i] + dir * (int)jdx + N) % N;
                size_t j = sort_idx[pos];

                double dx = x[j][0] - x[i][0];
                if (walls[0] < 10){dx -= std::round(dx / box[0]) * box[0];}
                if (std::abs(dx) > r_c_max) { go = false; break; }

                double r_c = (D[i] + D[j]) / 2 + t;
                if (std::abs(dx) < r_c) {
                    bool proceed = true;
                    double d2 = dx * dx;
                    for (size_t d = 1; d < Ndim; ++d) {
                        double dz = x[j][d] - x[i][d];
                        if (walls[d] < 10){dz -= std::round(dz / box[d]) * box[d];}
                        if (std::abs(dz) > r_c) { proceed = false; break; }
                        d2 += dz * dz;
                    }
                    if (!proceed) continue;

                    uint32_t &cnt = pairs[i][0];
                    if (cnt < neighborMax) {
                        ++cnt;
                        pairs[i][cnt] = j;
                    } else if (!warningIssued) {
                        std::cerr << "Error: found more neighbors than allotted ("
                                  << neighborMax
                                  << "). Restart with a larger --NeighborMax.\n";
                        std::exit(EXIT_FAILURE);
                    }

                    // mirror to j's list
                    if (cnt <= neighborMax) {
                        bool duplicate = false;
                        for (size_t k = 1; k <= pairs[j][0]; ++k)
                            if (pairs[j][k] == i) { duplicate = true; break; }
                        if (!duplicate) {
                            uint32_t &cntj = pairs[j][0];
                            if (cntj < neighborMax) {
                                ++cntj;
                                pairs[j][cntj] = i;
                            } else if (!warningIssued) {
                                std::cerr << "Error: found more neighbors than allotted ("
                                          << neighborMax
                                          << "). Restart with a larger --NeighborMax.\n";
                                std::exit(EXIT_FAILURE);
                            }
                        }
                    }
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
// GetForcesND_3 implementation
//------------------------------------------------------------------------------
void GetForcesND_3(
    const std::vector<std::vector<uint32_t>> &pairs,
    const std::vector<std::vector<double>>   &x,
    const std::vector<double>                &D,
    const std::vector<double>                &box,
    const std::vector<int8_t>                  &walls,
    std::vector<std::vector<double>>         &F,
    double                                   &U,
    std::vector<std::vector<double>>         &min_dist,
    double                                   &max_min_dist, 
    double                                   &Lc,
    std::vector<size_t>                      &z)
{
    size_t N      = pairs.size();
    size_t Ndim   = x[0].size();
    size_t Mplus1 = pairs[0].size();

    // init outputs
    F.assign(N, std::vector<double>(Ndim, 0.0));
    min_dist.assign(N, std::vector<double>(Mplus1, 0.0));
    z.assign(N, 0);
    U = 0.0;
    Lc = 0.0;
    size_t count = 1;

    std::vector<double> dx(Ndim);

    max_min_dist = 0;
    bool circle_flag = (walls[0] < 0);

    for (size_t i = 0; i < N; ++i) {
        uint32_t numNbr = pairs[i][0];
        for (uint32_t jdx = 1; jdx <= numNbr; ++jdx) {
            size_t j = pairs[i][jdx];
            if (j <= i) continue;

            double r_ij = 0.5 * (D[i] + D[j]);
            bool   flag = true;
            double d2   = 0.0;

            // PBC wrap using floor(m + 0.5) to match MATLAB round
            for (size_t d = 0; d < Ndim; ++d) {
                double delta = x[j][d] - x[i][d];
                if (walls[d] < 10){
                    double m     = delta / box[d];
                    delta       -= std::floor(m + 0.5) * box[d];
                }
                dx[d]        = delta;
                if (flag) {
                    if (std::abs(delta) > r_ij) {
                        flag = false;
                    }
                    d2 += delta * delta;
                }
            }

            if (flag) {
                double dist = std::sqrt(d2);
                if (dist < r_ij) {
                    // update contacts & energy
                    z[i]++; 
                    z[j]++;
                    double F_mag = - (r_ij / dist - 1.0);
                    U          += (1.0 - dist/r_ij);
                    ++count;
                    min_dist[i][jdx] = (1.0 - dist/r_ij);
                    Lc             += r_ij;
                    if (min_dist[i][jdx] > max_min_dist){max_min_dist=min_dist[i][jdx];};

                    // update forces
                    for (size_t d = 0; d < Ndim; ++d) {
                        double fcomp = F_mag * dx[d];
                        F[i][d]     += fcomp;
                        F[j][d]     -= fcomp;
                    }
               }
            }
        }

        //Boundary
        //Loop dimensions
        for (size_t d = 0; d < Ndim; ++d) {
            if (walls[d] == 1)
            {
                double r_ij = D[i]/2;
                //Loop two ends of box
                for (size_t wall = 0; wall < 2; ++wall) {
                    dx[d] = box[d]*wall - x[i][d];
                    double dist = std::abs(dx[d]);
                    if (dist < r_ij) {
                        // update contacts & energy
                        double F_mag = - 2.0*(r_ij / dist - 1.0); 
                        U          += 2.0*(1.0 - dist/r_ij);
                        ++count;
                        // update forces
                        double fcomp = F_mag * dx[d]; //since dx is half the size spring constant K = 2.0
                        F[i][d]     += fcomp;
                        }
                }
            }
        }
        
        if (circle_flag)
        {
            double r_ij = D[i]/2;
            double d_ij = 0;
            double R = box[0]/2;
            for (size_t M = 0; M < -walls[0]; ++M)
            {
                d_ij = d_ij + (x[i][M]-R)*(x[i][M]-R);
            }
            d_ij = std::sqrt(d_ij);
            double delta = R - d_ij;
            if (delta < r_ij)
            {
                delta = d_ij;
                d_ij = 0;
                for (size_t M = 0; M < -walls[0]; ++M)
                {
                    dx[M] = std::abs(R/delta - 1)*(x[i][M]-R);
                    d_ij = d_ij + dx[M]*dx[M];
                }
                d_ij = sqrt(d_ij);

                // Compute pairwise force magnitude
                double F_mag = - 2.0*std::abs(r_ij / d_ij - 1.0);  //abs ensures force always two center of container even if particle were outside container

                // Compute potential energy contribution
                U = U + (1 - d_ij/r_ij);
                ++count;

                // Update forces in all dimensions
                for (size_t M = 0; M < -walls[0]; ++M)
                {
                    F[i][M] = F[i][M] + F_mag * dx[M]; // Force on particle i
                }
        }
        }


    }

    U  = std::pow(U / double(count), 2);
    Lc = Lc / double(count);
}


//------------------------------------------------------------------------------
// AdamStep implementation
//------------------------------------------------------------------------------
void AdamUpdate(
    const std::string &method,
    size_t N,
    size_t Ndim,
    const std::vector<std::vector<double>> &F,
    double beta1,
    double beta2,
    size_t t,
    double dt,
    double verlet_drag,
    std::vector<std::vector<double>> &m,
    std::vector<std::vector<double>> &v,
    std::vector<std::vector<double>> &v_max,
    std::vector<std::vector<double>> &v_update,
    std::vector<std::vector<double>> &m_hat,
    std::vector<std::vector<double>> &v_hat,
    std::vector<std::vector<double>> &a,
    std::vector<std::vector<double>> &v_verlet,
    std::vector<std::vector<double>> &a_old) {
    for(size_t k=0;k<N;++k){
        for(size_t kk=0;kk<Ndim;++kk){ size_t idx=k*Ndim+kk; double Fkd=F[k][kk];
            m[k][kk]=beta1*m[k][kk] - (1.0-beta1)*Fkd;
            v[k][kk]=beta2*v[k][kk] + (1.0-beta2)*Fkd*Fkd;
            v_max[k][kk]=std::max(v_max[k][kk],v[k][kk]);
            if(method=="AMSGrad") v_update[k][kk]=v_max[k][kk]; else v_update[k][kk]=v[k][kk];
            m_hat[k][kk]=m[k][kk]/(1.0-std::pow(beta1,double(t)));
            v_hat[k][kk]=v_update[k][kk]/(1.0-std::pow(beta2,double(t)));
            a[k][kk] = Fkd - verlet_drag*v_verlet[k][kk];
            v_verlet[k][kk] += 0.5*(a_old[k][kk]+a[k][kk])*dt;
            a_old[k][kk] = a[k][kk];
        }
    }
}
