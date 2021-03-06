#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <any>
#include <chrono>
#include <list>
#include <set>
#include <random>
#include <array>
#include <functional>

std::mt19937 random_generator;

auto arg = [](int argc, char **argv, std::string name, auto default_value) -> decltype(default_value)
{
    using namespace std;
    string paramname = "";
    any ret = default_value;
    for (auto argument : vector<string>(argv, argv + argc))
    {
        if ((argument.size() > 0) && (argument[0] == '-'))
        {
            if (paramname != "")
            {
                if (name == argument.substr(1))
                    ret = true;
            }
            paramname = argument.substr(1);
        }
        else if (name == paramname)
        {
            if (std::is_same_v<decltype(default_value), int>)
                ret = stoi(argument);
            else if (std::is_same_v<decltype(default_value), double>)
                ret = stod(argument);
            else if (std::is_same_v<decltype(default_value), char>)
                ret = argument.at(0);
            else if (std::is_same_v<decltype(default_value), bool>)
                ret = (argument == "true") || (argument == "1") || (argument == "yes");
            else
                ret = argument;
            paramname = "";
        }
    }
    return std::any_cast<decltype(default_value)>(ret);
};

void show_help()
{
#ifndef __OPTIMIZE__
#define OPTIMIZED "optimized "
#else
#define OPTIMIZED ""
#endif
    std::cout << "# " << OPTIMIZED << "binary with date: " << __TIME__ << " " << __DATE__ << std::endl;
    std::cout << "-fname filename" << std::endl;
    std::cout << "-show_progress true/false" << std::endl;
}

using coordinates_t = std::array<double, 2>;
using solution_t = std::vector<int>;
using problem_t = std::vector<coordinates_t>;
/**
 * @brief helper operator for printing the solution
 *
 * @param o stream
 * @param sol solution to print on stream
 * @return std::ostream& the same stream as o
 */
inline std::ostream &operator<<(std::ostream &o, solution_t sol)
{
    for (int i = 0; i < sol.size(); i++)
    {
        o << (i == 0 ? "" : ", ") << sol[i];
    }
    return o;
}

/**
 * @brief substracting two vectors
 *
 * @param a
 * @param b
 * @return coordinates_t
 */
inline coordinates_t operator-(coordinates_t a, coordinates_t b)
{
    return {a[0] - b[0], a[1] - b[1]};
}
/**
 * @brief vector length
 *
 * @param p
 * @return double
 */
inline double length(coordinates_t p)
{
    return std::sqrt(p[0] * p[0] + p[1] * p[1]);
}

/**
 * @brief the cost function for the Travelling Salesman Problem. It calculates the cost of following the selected order in the list of cities represented by the problem.
 *
 * @param problem list of the cities
 * @param route the visitation order of cities
 * @return double the cost of following the route and comming back to the first city
 */
inline double tsp_problem_cost_function(const problem_t problem, const solution_t route)
{
    double distance = 0;
    for (int i = 0; i < route.size(); i++)
    {
        auto p1 = problem.at(route[i]);
        auto p2 = problem.at(route[(i + 1) % route.size()]);
        distance += length(p2 - p1);
    }
    return distance;
};

/**
 * @brief load the TSP problem from file or standard input. Each line represents one city. The city is represented by its coordinates.
 *
 * @param fname input filename. If this is empty, then we load from standard input
 * @return problem_t loaded array representing cities coordinates.
 */
inline problem_t load_problem(std::string fname)
{
    using namespace std;
    problem_t problem;
    istream *fst = &cin;
    ifstream f;
    if (fname != "")
    {
        f = ifstream(fname);
        fst = &f;
    }
    string line;
    while (getline(*fst, line))
    {
        stringstream sline(line);
        double x, y;
        sline >> x >> y;
        problem.push_back({x, y});
    }
    return problem;
}

inline std::string vis(solution_t sol, problem_t problem)
{
    using namespace std;
    /*
    graph G {
    A [ pos = "0,0!"]
    B [ pos = "10,10!"]
    A -- B
    C -- A
    }
    */
    stringstream out;
    out << "graph G {" << endl;
    int num = 0;
    for (auto [x, y] : problem)
    {
        out << num << " [ pos = \"" << x << "," << y << "!\" ] " << endl;
        num++;
    }
    for (int i = 0; i < sol.size(); i++)
    {
        out << sol[i] << " -- " << sol[(i + 1) % sol.size()] << endl;
    }

    out << "}" << endl;
    return out.str();
}

/**
 * @brief full review method that will use goal as a measure how bad is the solution (minimization).
 *
 * @tparam T the type of the solution
 * @param start_solution first solution to check
 * @param cost the goal function. Here it is cost function
 * @param next_solution method that will generate the next solution based on the current one
 * @param term_condition terminate iterating when this value is false
 * @param progress show progress on the console per iteration in the format suitable for gnuplot
 * @return T the best solution found so far
 */
template <typename T>
T calculate_full_review(T start_solution, std::function<double(T)> cost,
                        std::function<T(T)> next_solution, std::function<bool(T)> term_condition, bool progress = false)
{
    using namespace std;
    auto solution = start_solution;
    auto best_solution = solution;
    int iteration_counter = 0;
    do
    {
        if ((cost(solution) < cost(best_solution)) || (iteration_counter == 0))
        {
            if (progress)
                cout << (iteration_counter - 1) << " " << cost(best_solution) << " # " << solution << endl;
            best_solution = solution;
            if (progress)
                cout << iteration_counter << " " << cost(best_solution) << " # " << solution << endl;
        }
        iteration_counter++;
        solution = next_solution(solution);
    } while (term_condition(solution));
    if (progress)
        cout << iteration_counter << " " << cost(best_solution) << " # " << solution << endl;
    return best_solution;
}

std::pair<solution_t, std::chrono::duration<double>> experiment_full_revision(problem_t problem, bool show_progress)
{
    using namespace std;
    /// generate initial solution
    solution_t solution;
    for (int i = 0; i < problem.size(); i++)
    {
        solution.push_back(i);
    }

    /// prepare important functions:
    /// goal function for solution - how good or bad it is. We minimize this function
    function<double(solution_t)> goal = [problem](auto s)
    {
        return tsp_problem_cost_function(problem, s);
    };
    /// generate next solution for the method
    function<solution_t(solution_t)> next_solution =
        [](auto s)
    {
        next_permutation(s.begin(), s.end());
        return s;
    };
    /// what is the termination condition
    function<bool(solution_t)> term_condition = [solution](auto s)
    {
        return !(s == solution);
    };
    /// run the full review method
    auto calculation_start = chrono::steady_clock::now();
    solution_t best_solution = calculate_full_review<solution_t>(solution, goal, next_solution, term_condition, show_progress);
    auto calculation_end = chrono::steady_clock::now();
    /// show the result
    chrono::duration<double> calculation_duration = calculation_end - calculation_start;
    return {best_solution, calculation_duration};
}

solution_t generate_random_solution(problem_t problem)
{
    solution_t solution;
    for (int i = 0; i < problem.size(); i++)
    {
        solution.push_back(i);
    }

    std::shuffle(solution.begin(), solution.end(), random_generator);

    return solution;
}

std::pair<solution_t, std::chrono::duration<double>> experiment_random_hillclimb(problem_t problem,
                                                                                 int iterations, bool show_progress)
{
    using namespace std;
    /// generate initial solution
    solution_t solution = generate_random_solution(problem);

    /// prepare important functions:
    /// goal function for solution - how good or bad it is. We minimize this function
    function<double(solution_t)> goal = [problem](auto s)
    {
        return tsp_problem_cost_function(problem, s);
    };
    /// generate next solution for the method
    function<solution_t(solution_t)> next_solution = [problem](auto s)
    {
        uniform_int_distribution<int> distr(0, problem.size() - 1);
        int a = distr(random_generator);
        int b = distr(random_generator);
        swap(s[a], s[b]);
        return s;
    };
    /// what is the termination condition
    int i = 0;
    function<bool(solution_t)> term_condition = [&](auto s)
    {
        i++;
        return i < iterations;
    };
    /// run the full review method
    auto calculation_start = chrono::steady_clock::now();
    solution_t best_solution = calculate_full_review<solution_t>(solution,
                                                                 goal, next_solution,
                                                                 term_condition, show_progress);
    auto calculation_end = chrono::steady_clock::now();
    /// show the result
    chrono::duration<double> calculation_duration = calculation_end - calculation_start;
    return {best_solution, calculation_duration};
}

template <typename T>
T calculate_simulated_annealing(T start_solution,
                                std::function<double(T)> cost,
                                std::function<T(T)> next_solution,
                                std::function<double(int)> temperature,
                                std::function<bool(T)> term_condition,
                                bool progress = false)
{
    using namespace std;
    auto solution = start_solution;
    auto best_solution = solution;
    int iteration_counter = 1;
    do
    {
        auto neighbor = next_solution(solution);
        if (cost(neighbor) < cost(solution))
        {
            solution = neighbor;
        }
        else
        {
            uniform_real_distribution<double> u(0, 1);
            if (u(random_generator) < exp(-abs(cost(solution) - cost(neighbor)) / temperature(iteration_counter)))
            {
                solution = neighbor;
            }
        }
        if (cost(solution) < cost(best_solution))
            best_solution = solution;
        iteration_counter++;
        if (progress)
            cout << iteration_counter << " " << cost(solution) << " " << cost(best_solution) << " # SOL: " << solution << endl;

    } while (term_condition(solution));
    if (progress)
        cout << iteration_counter << " " << cost(best_solution) << " # " << solution << endl;
    return best_solution;
}

std::pair<solution_t, std::chrono::duration<double>> experiment_simulated_annealing(problem_t problem,
                                                                                    int iterations,
                                                                                    bool show_progress)
{
    using namespace std;
    /// generate initial solution
    solution_t solution = generate_random_solution(problem);

    /// prepare important functions:
    /// goal function for solution - how good or bad it is. We minimize this function
    function<double(solution_t)> goal = [problem](auto s)
    {
        return tsp_problem_cost_function(problem, s);
    };
    /// generate next solution for the method
    function<solution_t(solution_t)> next_solution = [problem](auto s)
    {
        uniform_int_distribution<int> distr(0, problem.size() - 1);
        int a = distr(random_generator);
        int b = distr(random_generator);
        swap(s[a], s[b]);
        return s;
    };
    /// what is the termination condition
    int i = 0;
    function<bool(solution_t)> term_condition = [&](auto s)
    {
        i++;
        return i < iterations;
    };

    function<double(int)> temperature = [](int t)
    {
        return 1000 / t;
    };
    /// run the full review method
    auto calculation_start = chrono::steady_clock::now();
    solution_t best_solution = calculate_simulated_annealing<solution_t>(solution,
                                                                         goal, next_solution,
                                                                         temperature,
                                                                         term_condition,
                                                                         show_progress);
    auto calculation_end = chrono::steady_clock::now();
    /// show the result
    chrono::duration<double> calculation_duration = calculation_end - calculation_start;
    return {best_solution, calculation_duration};
}

std::pair<solution_t, std::chrono::duration<double>> experiment_random_probe(problem_t problem,
                                                                             int iterations, bool show_progress)
{
    using namespace std;
    /// generate initial solution
    solution_t solution = generate_random_solution(problem);

    /// prepare important functions:
    /// goal function for solution - how good or bad it is. We minimize this function
    function<double(solution_t)> goal = [problem](auto s)
    {
        return tsp_problem_cost_function(problem, s);
    };
    /// generate next solution for the method
    function<solution_t(solution_t)> next_solution = [problem](auto s)
    {
        return generate_random_solution(problem);
    };
    /// what is the termination condition
    int i = 0;
    function<bool(solution_t)> term_condition = [&](auto s)
    {
        i++;
        return i < iterations;
    };
    /// run the full review method
    auto calculation_start = chrono::steady_clock::now();
    solution_t best_solution = calculate_full_review<solution_t>(solution,
                                                                 goal, next_solution,
                                                                 term_condition, show_progress);
    auto calculation_end = chrono::steady_clock::now();
    /// show the result
    chrono::duration<double> calculation_duration = calculation_end - calculation_start;
    return {best_solution, calculation_duration};
}

template <typename T>
T calculate_tabu(T start_solution,
                 std::function<double(T)> cost,
                 std::function<std::vector<T>(T)> get_neighbors,
                 std::function<bool(T)> term_condition,
                 int max_tabu_size = 100,
                 bool progress = false)
{
    using namespace std;
    auto solution = start_solution;
    auto best_solution = solution;
    std::list<T> tabu_list;
    std::set<T> tabu_set;
    int iteration_counter = 0;
    auto is_not_in_tabu = [&](T s) -> bool
    {
        return tabu_set.count(s) == 0;
    };

    do
    {
        auto neighbors = get_neighbors(solution);
        vector<T> good_neighbors;
        copy_if(neighbors.begin(), neighbors.end(), std::back_inserter(good_neighbors), is_not_in_tabu);
        if (good_neighbors.size() == 0)
        {
            cerr << "empty list..." << endl;
            break;
        }
        solution = *min_element(good_neighbors.begin(), good_neighbors.end(), [&](auto a, auto b) -> int
                                { return (cost(a) < cost(b)); });
        tabu_list.push_back(solution);
        tabu_set.insert(solution);
        if (tabu_list.size() > max_tabu_size)
        {
            tabu_set.erase(tabu_list.front());
            tabu_list.pop_front();
        }
        if ((cost(solution) < cost(best_solution)) || (iteration_counter == 0))
        {
            if (progress)
                cout << (iteration_counter - 1) << " " << cost(best_solution) << " # " << solution << endl;
            best_solution = solution;
            if (progress)
                cout << iteration_counter << " " << cost(best_solution) << " # " << solution << endl;
        }
        iteration_counter++;
    } while (term_condition(solution));
    if (progress)
        cout << iteration_counter << " " << cost(best_solution) << " # " << solution << endl;
    return best_solution;
}

std::pair<solution_t, std::chrono::duration<double>> experiment_tabu(problem_t problem,
                                                                     int iterations, int max_tabu_size, bool show_progress)
{
    using namespace std;
    /// generate initial solution
    solution_t solution = generate_random_solution(problem);

    /// prepare important functions:
    /// goal function for solution - how good or bad it is. We minimize this function
    function<double(solution_t)> goal = [problem](auto s)
    {
        return tsp_problem_cost_function(problem, s);
    };
    /// generate next solution for the method
    function<vector<solution_t>(solution_t)> get_neighbors = [problem](auto s) -> vector<solution_t>
    {
        vector<solution_t> n;
        for (int i = 0; i < s.size(); i++)
        {
            auto tmp = s;
            swap(tmp[i], tmp[(i + 1) % s.size()]);
            n.push_back(tmp);
        }
        return n;
    };
    /// what is the termination condition
    int i = 0;
    function<bool(solution_t)> term_condition = [&](auto s)
    {
        i++;
        return i < iterations;
    };
    /// run the full review method
    auto calculation_start = chrono::steady_clock::now();
    solution_t best_solution = calculate_tabu<solution_t>(solution,
                                                          goal, get_neighbors,
                                                          term_condition, max_tabu_size, show_progress);
    auto calculation_end = chrono::steady_clock::now();
    /// show the result
    chrono::duration<double> calculation_duration = calculation_end - calculation_start;
    return {best_solution, calculation_duration};
}

int selection_tournament(std::vector<double> fitnesses)
{
    std::uniform_int_distribution<int> distr(0, fitnesses.size() - 1);
    int idx1 = distr(random_generator);
    int idx2 = distr(random_generator);
    if (fitnesses[idx1] > fitnesses[idx2])
        return idx1;
    return idx2;
}

template <typename T>
T calculate_genetic_algorithm(std::vector<T> initial_population,
                              std::function<double(T)> fitness,
                              std::function<int(std::vector<double>)> selection,
                              std::function<std::pair<T, T>(T, T)> crossover,
                              std::function<T(T)> mutation,
                              std::function<bool(std::vector<double>, std::vector<T>)> term_condition,
                              bool progress = false)
{
    using namespace std;
    auto population = initial_population;
    auto calculate_fitness = [&](auto pop)
    {
        std::vector<double> ret;
        for (auto e : pop)
            ret.push_back(fitness(e));
        return ret;
    };
    auto population_fit = calculate_fitness(population);

    while (term_condition(population_fit, population))
    {
        vector<T> parents;
        vector<T> offspring;

        for (int i = 0; i < population.size(); i++)
        {
            parents.push_back(population[selection(population_fit)]);
        }
        for (int i = 0; i < parents.size(); i += 2)
        {
            auto [a, b] = crossover(parents[i], parents[i + 1]);
            offspring.push_back(a);
            offspring.push_back(b);
        }
        for (int i = 0; i < parents.size(); i++)
        {
            offspring[i] = mutation(offspring[i]);
        }

        population = offspring;
        population_fit = calculate_fitness(population);
    };
    //if (progress) cout << iteration_counter << " " << cost(best_solution) << " # " << solution << endl;
    return population[max_element(population_fit.begin(), population_fit.end()) - population_fit.begin()];
}

std::pair<solution_t, std::chrono::duration<double>> experiment_genetic_algorithm(problem_t problem,
                                                                                  int iterations,
                                                                                  int population_size,
                                                                                  double crossover_probability,
                                                                                  double mutation_probability,
                                                                                  bool show_progress)
{
    using namespace std;
    /// generate initial solution
    vector<solution_t> init_population;
    for (int i = 0; i < population_size; i++)
        init_population.push_back(generate_random_solution(problem));

    /// prepare important functions:
    /// goal function for solution - how good or bad it is. We minimize this function
    function<double(solution_t)> fitness_function = [problem](auto s)
    {
        map<int, int> stats;
        double error = 0;
        for (auto i : s)
        {
            stats[i] = stats[i] + 1;
        }
        for (auto [k, v] : stats)
        {
            if (v > 1)
                error += v;
        }
        double fit = 10000.0 + 1000.0 / (1 + tsp_problem_cost_function(problem, s)) - error * 100;
        return (fit > 0) ? fit : 0.0;
    };
    /// mutation method
    function<solution_t(solution_t)> mutation = [=](auto s)
    {
        uniform_real_distribution<double> u(0.0, 1.0);
        if (u(random_generator) < mutation_probability)
        {
            uniform_int_distribution<int> distr(0, problem.size() - 1);
            int a = distr(random_generator);
            int b = distr(random_generator);
            swap(s[a], s[b]);
        }
        return s;
    };
    auto crossover = [=](solution_t a, solution_t b) -> pair<solution_t, solution_t>
    {
        uniform_real_distribution<double> u(0.0, 1.0);
        if (u(random_generator) < crossover_probability)
        {
            uniform_int_distribution<int> cross_distr(0, a.size() - 1);
            int crossover_pt = cross_distr(random_generator);
            for (int i = crossover_pt; i < a.size(); i++)
                swap(a[i], b[i]);
        }
        return {a, b};
    };
    /// what is the termination condition
    int i = 0;
    std::function<bool(std::vector<double>, std::vector<solution_t>)>
        term_condition = [&](auto pop_fit, auto pop)
    {
        i++;
        return i < iterations;
    };

    /// run the full review method
    auto calculation_start = chrono::steady_clock::now();
    solution_t best_solution = calculate_genetic_algorithm<solution_t>(init_population,
                                                                       fitness_function,
                                                                       selection_tournament,
                                                                       crossover,
                                                                       mutation,
                                                                       term_condition,
                                                                       show_progress);
    auto calculation_end = chrono::steady_clock::now();
    /// show the result
    chrono::duration<double> calculation_duration = calculation_end - calculation_start;
    return {best_solution, calculation_duration};
}

int main(int argc, char **argv)
{
    using namespace std;

    /// get arguments
    auto fname = arg(argc, argv, "fname", std::string(""));
    auto show_progress = arg(argc, argv, "show_progress", false);
    auto iterations = arg(argc, argv, "iterations", 1000);
    auto method = arg(argc, argv, "method", std::string("random_hillclimb"));
    auto tabu_size = arg(argc, argv, "tabu_size", 100);
    auto dot = arg(argc, argv, "dot", false);

    auto pop_size = arg(argc, argv, "pop_size", 50);
    auto crossover_p = arg(argc, argv, "crossover_p", 0.9);
    auto mutation_p = arg(argc, argv, "mutation_p", 0.1);
    auto help = arg(argc, argv, "help", true);

    /// do we need to show some help?
    if (arg(argc, argv, "h", false))
    {
        show_help();
        return 0;
    }
    cout << "# fname = " << fname << ";" << std::endl;

    /// load the problem at hand
    auto problem = load_problem(fname);
    solution_t best_solution;
    std::chrono::duration<double> calculation_duration;
    if (method == "full_revision")
    {
        auto [a, b] = experiment_full_revision(problem, show_progress);
        best_solution = a;
        calculation_duration = b;
    }
    else if (method == "random_probe")
    {
        auto [a, b] = experiment_random_probe(problem, iterations, show_progress);
        best_solution = a;
        calculation_duration = b;
    }
    else if (method == "random_hillclimb")
    {
        auto [a, b] = experiment_random_hillclimb(problem, iterations, show_progress);
        best_solution = a;
        calculation_duration = b;
    }
    else if (method == "tabu")
    {
        auto [a, b] = experiment_tabu(problem, iterations, tabu_size, show_progress);
        best_solution = a;
        calculation_duration = b;
    }
    else if (method == "simulated_annealing")
    {
        auto [a, b] = experiment_simulated_annealing(problem, iterations, show_progress);
        best_solution = a;
        calculation_duration = b;
    }
    else if (method == "genetic_algorithm")
    {
        auto [a, b] = experiment_genetic_algorithm(problem, iterations,
                                                   pop_size,
                                                   crossover_p,
                                                   mutation_p,
                                                   show_progress);
        best_solution = a;
        calculation_duration = b;
    }
    else
    {
        std::cerr << "unknown method" << std::endl;
    }
    cout << "# " << method << ": best_cost: " << tsp_problem_cost_function(problem, best_solution) << " calculation_time: " << calculation_duration.count() << endl;
    cout << "# solution: " << best_solution << endl;
    if (dot)
        cout << vis(best_solution, problem);
    return 0;
}
