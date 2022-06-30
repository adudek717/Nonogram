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

using x_axis = std::vector<int>;
using y_axis = std::vector<int>;
using nono_problem_t = std::vector<std::vector<int>>;
using nono_solution_t = std::vector<std::vector<int>>;

/**
 * @brief helper operator for printing the solution
 *
 * @param o stream
 * @param sol solution to print on stream
 * @return std::ostream& the same stream as o
 */
inline std::ostream &operator<<(std::ostream &o, nono_solution_t sol)
{
    for (int i = 0; i < sol.size(); i++)
        for (int j = 0; j < sol.at(i).size(); j++)
            o << (i == 0 ? "" : ", ") << sol[i][j];
    return o;
}

/**
 * @brief load the TSP problem from file or standard input. Each line represents one city. The city is represented by its coordinates.
 *
 * @param fname input filename. If this is empty, then we load from standard input
 * @return problem_t loaded array representing cities coordinates.
 */
inline nono_problem_t nono_load_problem(std::string fname)
{
    using namespace std;
    nono_problem_t problem;
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
        int x;
        int i = 0;
        std::vector<int> row;
        while (sline >> x)
        {
            row.push_back(x);
        }
        problem.push_back(row);
    }
    return problem;
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
T calculate_full_review_nono(T start_solution, std::function<double(T)> cost,
                             std::function<T(T)> next_solution, std::function<bool(T)> term_condition, bool progress = false)
{
    using namespace std;
    auto solution = start_solution;
    auto best_solution = solution;
    int iteration_counter = 0;
    do
    {
        cout << "Solution cost: " << cost(solution) << endl;
        cout << "Best cost: " << cost(best_solution) << endl;
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

nono_solution_t generate_random_solution_nono(nono_problem_t problem)
{
    nono_solution_t solution;
    std::uniform_int_distribution<int> distr(0, 1);
    for (int i = 0; i < problem.size(); i++)
    {
        std::vector<int> nums;
        for (int j = 0; j < problem.at(0).size(); j++)
        {
            nums.push_back(distr(random_generator));
            // solution.at(i).at(j) = distr(random_generator);
        }
        solution.push_back(nums);
    }

    return solution;
}

/**
 * @brief the cost function for the Travelling Salesman Problem. It calculates the cost of following the selected order in the list of cities represented by the problem.
 *
 * @param problem list of the cities
 * @param route the visitation order of cities
 * @return double the cost of following the route and comming back to the first city
 */
inline int nono_problem_cost_function(const nono_problem_t problem, const nono_solution_t solution)
{
    int bad_cells = 0;
    if (problem.size() != solution.size())
    {
        std::cout << "Problem size and solution size are different!" << std::endl;
        std::cout << "Problem size: " << problem.size() << std::endl;
        std::cout << "Solution size: " << solution.size() << std::endl;
    }

    for (int i = 0; i < problem.size(); i++)
    {
        for (int j = 0; j < problem.at(i).size(); j++)
        {
            if (problem.at(i).size() != solution.at(i).size())
            {
                std::cout << "Problem size and solution size are different!" << std::endl;
                std::cout << "Problem at i size: " << problem.at(i).size() << std::endl;
                std::cout << "Solution at i size: " << solution.at(i).size() << std::endl;
            }
            if (problem.at(i).at(j) != solution.at(i).at(j))
                bad_cells++;
        }
    }
    return bad_cells;
};

std::pair<nono_solution_t, std::chrono::duration<double>> experiment_random_hillclimb_nono(nono_problem_t problem,
                                                                                           int iterations, bool show_progress)
{
    using namespace std;

    /// generate initial solution
    nono_solution_t solution = generate_random_solution_nono(problem);

    /// prepare important functions:
    /// goal function for solution - how good or bad it is. We minimize this function
    function<double(nono_solution_t)> goal = [problem](auto s)
    {
        return nono_problem_cost_function(problem, s);
    };

    /// generate next solution for the method
    function<nono_solution_t(nono_solution_t)> next_solution = [problem](auto s)
    {
        uniform_int_distribution<int> distr(0, problem.size() - 1);
        int a = distr(random_generator);
        int b = distr(random_generator);
        swap(s[a], s[b]);
        return s;
    };

    /// what is the termination condition
    int i = 0;
    function<bool(nono_solution_t)> term_condition = [&](auto s)
    {
        i++;
        return i < iterations;
    };
    nono_solution_t best_solution;

    /// run the full review method
    auto calculation_start = chrono::steady_clock::now();
    best_solution = calculate_full_review_nono<nono_solution_t>(solution,
                                                                goal, next_solution,
                                                                term_condition, show_progress);
    auto calculation_end = chrono::steady_clock::now();
    //
    ///// show the result
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
    auto help = arg(argc, argv, "help", true);

    /// do we need to show some help?
    if (arg(argc, argv, "h", false))
    {
        show_help();
        return 0;
    }
    cout << "# fname = " << fname << ";" << std::endl;

    /// load the problem at hand
    auto nono_problem = nono_load_problem(fname);
    nono_solution_t best_nono_solution;
    std::chrono::duration<double> calculation_duration;

    // Testing
    cout << endl
         << "Problem: " << endl;
    for (int i = 0; i < nono_problem.size(); i++)
    {
        for (int j = 0; j < nono_problem.at(i).size(); j++)
        {
            cout << nono_problem.at(i).at(j) << " ";
        }
        cout << endl;
    }
    cout << endl;

    if (method == "random_hillclimb_nono")
    {

        auto [a, b] = experiment_random_hillclimb_nono(nono_problem, iterations, show_progress);
        best_nono_solution = a;
        calculation_duration = b;
    }
    else
    {
        std::cerr << "unknown method" << std::endl;
    }

    cout << "# " << method << ": best_cost: " << nono_problem_cost_function(nono_problem, best_nono_solution) << " calculation_time: " << calculation_duration.count() << endl;
    cout << "# solution: " << best_nono_solution << endl;

    return 0;
}