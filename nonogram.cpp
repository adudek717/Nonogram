#include <iostream>
#include <string>
#include <any>
#include <vector>
#include <functional>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <climits>
#include <stack>

using namespace std;

random_device rd;
mt19937 random_generator(rd());

auto arg = [](int argc, char **argv, std::string name, auto default_value) -> decltype(default_value)
{
    // using namespace std;
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

vector<vector<int>> load_problem(const string &fname)
{
    cout << "Loading problem..." << endl;
    vector<vector<int>> problem;
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
        string text;
        string delimiter = ",";
        size_t pos = 0;
        string token;
        while (sline >> text)
        {
            vector<int> nums;
            if ((pos = text.find(delimiter)) != string::npos)
            {
                token = text.substr(0, pos);
                // cout << token << endl; // Test output
                nums.push_back(stoi(token));
                text.erase(0, pos + delimiter.length());
            }
            // cout << "TEXT: " << text << endl; // Test output
            nums.push_back(stoi(text));
            problem.push_back(nums);
        }
    }
    return problem;
}

vector<vector<bool>> load_valid_solution(const string &fname)
{
    cout << "Loading valid solution..." << endl;
    vector<vector<bool>> valid_solution;
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
        std::vector<bool> row;
        while (sline >> x)
        {
            if (x == 0)
                row.push_back(false);
            if (x == 1)
                row.push_back(true);
        }
        valid_solution.push_back(row);
    }
    return valid_solution;
}

vector<vector<bool>> generate_random_solution(const vector<vector<int>> &problem, const int rows, const int cols)
{
    vector<vector<bool>> random_solution;

    // Initialize random_solution vector
    for (int i = 0; i < rows; i++)
    {
        vector<bool> col(cols, 0);
        random_solution.push_back(col);
    }

    // Mark ones for the rows
    // Looping valid solution rows
    // i = problem vector   AND     i = valid_solution vector
    for (int i = 0; i < rows; i++) // 0-10
    {
        int offset = 0;
        int sum = 0;
        // Looping through nono problem(inner vector) values (for sum)
        for (int k = 0; k < problem.at(i).size(); k++)
        {
            sum += problem.at(i).at(k);
        }
        // range is a size
        int range = cols - sum;

        // range is a index
        if (--range < 0)
            range = 0;

        uniform_int_distribution<int> distr(0, range);
        offset = distr(random_generator);
        // cout << "Offset: " << offset << endl;

        // Index idx
        int idx = 0;
        while (idx < static_cast<int>(problem.at(0).size()))
        {
            // Looping through nono problem(inner vector)
            for (int k = 0; k < problem.at(i).size(); k++)
            {
                // Looping through the value of the element
                for (int l = 0; l < problem.at(i).at(k); l++)
                {
                    random_solution.at(i).at(idx + offset) = 1; // later l + offset
                    idx++;
                }
                idx++;
            }
        }
    }

    return random_solution;
}

vector<vector<bool>> generate_completely_random_solution(const vector<vector<int>> &problem, const int rows, const int cols)
{
    vector<vector<bool>> random_solution;

    // Initialize random_solution vector
    for (int i = 0; i < rows; i++)
    {
        vector<bool> col(cols, 0);
        random_solution.push_back(col);
    }

    for (int i = 0; i < rows; i++) // 0-10
    {
        for (int j = 0; j < cols; j++)
        {
            uniform_int_distribution<int> distr(0, 1);
            random_solution.at(i).at(j) = distr(random_generator);
        }
    }

    return random_solution;
}

vector<vector<int>> get_problem_rows_from_solution(const vector<vector<bool>> &solution)
{
    vector<vector<int>> problem;

    for (int i = 0; i < solution.size(); i++)
    {
        int sum = 0;
        vector<int> vec;
        for (int j = 0; j < solution.at(i).size(); j++)
        {
            bool current_num = solution.at(i).at(j);
            if (current_num != 0)
                sum++;

            if (vec.empty() && j == solution.at(i).size() - 1)
            {
                vec.push_back(sum);
                sum = 0;
            }

            if (current_num == 0 && sum != 0)
            {
                vec.push_back(sum);
                sum = 0;
            }
            if (current_num != 0 && j == solution.at(i).size() - 1 && sum != 0)
            {
                vec.push_back(sum);
                sum = 0;
            }
        }
        problem.push_back(vec);
    }

    return problem;
}

vector<vector<int>> get_problem_cols_from_solution(const vector<vector<bool>> &solution)
{
    vector<vector<int>> problem;

    for (int i = 0; i < solution.at(0).size(); i++)
    {
        int sum = 0;
        vector<int> vec;
        for (int j = 0; j < solution.size(); j++)
        {
            bool current_num = solution.at(j).at(i);
            if (current_num != 0)
                sum++;

            if (vec.empty() && j == solution.size() - 1)
            {
                vec.push_back(sum);
                sum = 0;
            }

            if (current_num == 0 && sum != 0)
            {
                vec.push_back(sum);
                sum = 0;
            }
            if (current_num != 0 && j == solution.size() - 1 && sum != 0)
            {
                vec.push_back(sum);
                sum = 0;
            }
        }
        problem.push_back(vec);
    }
    return problem;
}

vector<vector<vector<bool>>> generate_neighbors(const vector<vector<bool>> &solution)
{
    vector<vector<vector<bool>>> neighbors;
    for (int i = 0; i < solution.size(); i++)
    {
        for (int j = 0; j < solution.at(i).size(); j++)
        {
            vector<vector<bool>> neighbor = solution;
            neighbor.at(i).at(j) = !neighbor.at(i).at(j);
            // swap(neighbor.at(i).at(j), neighbor.at(i).at(j + 1));
            neighbors.push_back(neighbor);
        }
    }
    return neighbors;
}

// Return the cost of a solution
int get_solution_cost(const vector<vector<bool>> &solution, const vector<vector<int>> &problem_rows, const vector<vector<int>> &problem_cols)
{
    int overall_cost = 0;
    int rows_cost = 0;
    int cols_cost = 0;

    vector<vector<int>> solution_to_problem_rows = get_problem_rows_from_solution(solution);
    vector<vector<int>> solution_to_problem_cols = get_problem_cols_from_solution(solution);

    for (int i = 0; i < solution_to_problem_rows.size(); i++)
    {
        for (int j = 0; j < solution_to_problem_rows.at(i).size(); j++)
        {
            int stp_rows_size = solution_to_problem_rows.at(i).size();
            int p_rows_size = problem_rows.at(i).size();
            if (stp_rows_size == p_rows_size)
            {
                if (solution_to_problem_rows.at(i).at(j) != problem_rows.at(i).at(j))
                    rows_cost++;
            }
            else
            {
                int size_difference = 0;
                if (stp_rows_size > p_rows_size)
                    size_difference = stp_rows_size - p_rows_size;
                else
                    size_difference = p_rows_size - stp_rows_size;
                rows_cost += size_difference;
                // rows_cost++;
            }
        }
    }

    for (int i = 0; i < solution_to_problem_cols.size(); i++)
    {
        for (int j = 0; j < solution_to_problem_cols.at(i).size(); j++)
        {
            int stp_cols_size = solution_to_problem_cols.at(i).size();
            int p_cols_size = problem_cols.at(i).size();
            if (stp_cols_size == p_cols_size)
            {
                if (solution_to_problem_cols.at(i).at(j) != problem_cols.at(i).at(j))
                    cols_cost++;
            }
            else
            {
                int size_difference = 0;
                if (stp_cols_size > p_cols_size)
                    size_difference = stp_cols_size - p_cols_size;
                else
                    size_difference = p_cols_size - stp_cols_size;
                cols_cost += size_difference;
                // cols_cost++;
            }
        }
    }

    overall_cost = rows_cost + cols_cost;

    return overall_cost;
}

void print_problem(const vector<vector<int>> &problem)
{
    for (auto vec : problem)
    {
        for (auto num : vec)
        {
            cout << num << " ";
        }
        cout << endl;
    }
}

void print_solution(const vector<vector<bool>> &solution)
{
    for (auto vec : solution)
    {
        for (auto num : vec)
        {
            if (num == 0)
                cout << "-";
            else
                cout << '#';
        }
        cout << endl;
    }
}

vector<vector<bool>> get_best_neighbor(const vector<vector<vector<bool>>> &neighbors, const vector<vector<int>> &problem_rows, const vector<vector<int>> &problem_cols)
{
    vector<vector<bool>> best_neighbor = neighbors.at(0);
    int best_cost = get_solution_cost(best_neighbor, problem_rows, problem_cols);
    int current_cost = best_cost;

    for (int i = 1; i < neighbors.size(); i++)
    {
        current_cost = get_solution_cost(neighbors.at(i), problem_rows, problem_cols);
        if (current_cost < best_cost)
        {
            best_cost = current_cost;
            best_neighbor = neighbors.at(i);
        }
    }

    return best_neighbor;
}

vector<vector<bool>> get_hillclimbing_solution(const vector<vector<bool>> &solution, const vector<vector<int>> &problem_rows, const vector<vector<int>> &problem_cols, int iterations)
{
    vector<vector<bool>> best_solution = solution;
    int current_iteration = 0;
    int best_cost = get_solution_cost(solution, problem_rows, problem_cols);

    while (best_cost > 0 && current_iteration <= iterations)
    {
        vector<vector<vector<bool>>> neighbors = generate_neighbors(best_solution);
        // cout << "Iteration: " << current_iteration << endl;
        // for (auto nei : neighbors)
        //{
        //     print_solution(nei);
        //     cout << get_solution_cost(nei, problem_rows, problem_cols) << endl;;
        // }
        //

        auto new_solution = get_best_neighbor(neighbors, problem_rows, problem_cols);
        auto new_cost = get_solution_cost(new_solution, problem_rows, problem_cols);
        auto current_cost = get_solution_cost(best_solution, problem_rows, problem_cols);
        if (new_cost <= current_cost)
        {
            best_solution = new_solution;
        }
        else
            break;

        // cout << "Best solution cost: " << best_cost << endl;
        current_iteration++;
    }
    return best_solution;
}

vector<vector<bool>> get_hillclimbing_random_solution(const vector<vector<bool>> &solution, const vector<vector<int>> &problem_rows, const vector<vector<int>> &problem_cols, int iterations)
{
    vector<vector<bool>> best_solution = solution;
    int current_iteration = 0;
    int best_cost = get_solution_cost(solution, problem_rows, problem_cols);

    while (best_cost > 0 && current_iteration <= iterations)
    {
        vector<vector<vector<bool>>> neighbors = generate_neighbors(best_solution);
        vector<int> costs;
        vector<vector<int>> costs_indices;
        int cost;
        int sum = 0;
        // cout << "Iteration: " << current_iteration << endl;
        // for (auto nei : neighbors)
        //{
        //     print_solution(nei);
        //     cout << get_solution_cost(nei, problem_rows, problem_cols) << endl;;
        // }

        for (auto nei : neighbors)
        {
            cost = get_solution_cost(nei, problem_rows, problem_cols);
            if (cost == 0)
                return nei;
            sum += cost;
            costs.push_back(cost);
        }

        for (int i = 0; i < costs.size(); i++)
        {
            costs_indices.push_back({sum - costs.at(i), i});
        }

        sort(costs_indices.begin(), costs_indices.end());
        reverse(costs_indices.begin(), costs_indices.end());

        uniform_int_distribution<int> distr(0, sum);
        int random = distr(random_generator);

        int pick_idx = 0;
        for (auto ci : costs_indices)
        {
            random -= ci.at(0);
            if (random <= 0)
            {
                pick_idx = ci.at(1);
                break;
            }
        }

        best_solution = neighbors.at(pick_idx);
        best_cost = get_solution_cost(neighbors.at(pick_idx), problem_rows, problem_cols);

        // for (auto ci : costs_indices)
        //{
        //     cout << ci.at(0) << "  " << ci.at(1) << endl;
        // }

        // uniform_int_distribution<int> distr(0, neighbors.size() - 1);
        // int idx1 = distr(random_generator);
        // int idx2 = distr(random_generator);
        //
        // vector<vector<bool>> random_neighbor1 = neighbors.at(idx1);
        // vector<vector<bool>> random_neighbor2 = neighbors.at(idx2);
        //
        // int random_neighbor1_cost = get_solution_cost(random_neighbor1, problem_rows, problem_cols);
        // int random_neighbor2_cost = get_solution_cost(random_neighbor2, problem_rows, problem_cols);
        //
        //
        // if (random_neighbor1_cost < random_neighbor2_cost)
        //{
        //     best_solution = random_neighbor1;
        // }
        // else
        //{
        //     best_solution = random_neighbor2;
        // }
        //
        // auto new_solution = get_best_neighbor(neighbors, problem_rows, problem_cols);
        // auto new_cost = get_solution_cost(new_solution, problem_rows, problem_cols);
        // auto current_cost = get_solution_cost(best_solution, problem_rows, problem_cols);
        // if (new_cost <= current_cost)
        //{
        //     best_solution = new_solution;
        // }
        // else
        //     break;

        // cout << "Best solution cost: " << best_cost << endl;
        current_iteration++;
    }
    return best_solution;
}

vector<vector<bool>> get_tabu_solution(const vector<vector<bool>> &solution, const vector<vector<int>> &problem_rows, const vector<vector<int>> &problem_cols, int iterations)
{
    vector<vector<bool>> best_solution = solution;
    int best_cost = get_solution_cost(solution, problem_rows, problem_cols);
    vector<vector<vector<bool>>> tabu;
    tabu.push_back(solution);
    stack<vector<vector<bool>>> steps;
    steps.push(solution);
    int current_iteration = 0;

    while (best_cost > 0 && current_iteration <= iterations && steps.size() > 0)
    {
        vector<vector<vector<bool>>> neighbors = generate_neighbors(best_solution);
        vector<vector<vector<bool>>> valid_neighbors;
        for (auto nei : neighbors)
        {
            // steps.push(nei);
            // tabu.push_back(nei);
            if (find(tabu.begin(), tabu.end(), nei) == tabu.end())
            {
                valid_neighbors.push_back(nei);
            }
        }

        vector<vector<bool>> best_neighbor;
        int best_neighbor_cost = INT_MAX;
        for (auto vnei : valid_neighbors)
        {
            int cost = get_solution_cost(vnei, problem_rows, problem_cols);
            if (cost < best_neighbor_cost)
            {
                best_neighbor_cost = cost;
                best_neighbor = vnei;
            }
        }

        if (best_neighbor_cost < best_cost)
        {
            best_solution = best_neighbor;
            best_cost = best_neighbor_cost;
            steps.push(best_solution);
            tabu.push_back(best_solution);
        }
        else
        {
            best_solution = steps.top();
            steps.pop();
        }

        // cout << "Best solution cost: " << best_cost << endl;
        current_iteration++;
    }
    return best_solution;
}

// vector<vector<bool>> get_hillclimbing_random_solution(const vector<vector<bool>>& solution, const vector<vector<int>>& problem_rows, const vector<vector<int>>& problem_cols, int iterations)
//{
//     vector<vector<bool>> best_solution = solution;
//     int current_iteration = 0;
//     int best_cost = get_solution_cost(solution, problem_rows, problem_cols);
//     int current_cost = best_cost;
//
//     while (best_cost > 0 && current_iteration <= iterations)
//     {
//         vector<vector<vector<bool>>> neighbors = generate_neighbors(best_solution);
//         //for (auto nei : neighbors)
//         //{
//         //    print_solution(nei);
//         //    cout << get_solution_cost(nei, problem_rows, problem_cols) << endl;;
//         //}
//         for (int i = 0; i < neighbors.size(); i++)
//         {
//             uniform_int_distribution<int> distr(0, neighbors.size() - 1);
//             int idx = distr(random_generator);
//             current_cost = get_solution_cost(neighbors.at(idx), problem_rows, problem_cols);
//             if (current_cost < best_cost)
//             {
//                 best_cost = current_cost;
//                 best_solution = neighbors.at(idx);
//                 break;
//             }
//         }
//         for (int i = 0; i < neighbors.size(); i++)
//         {
//             current_cost = get_solution_cost(neighbors.at(i), problem_rows, problem_cols);
//             if (current_cost < best_cost)
//             {
//                 best_cost = current_cost;
//                 best_solution = neighbors.at(i);
//                 break;
//             }
//         }
//         current_iteration++;
//     }
//     return best_solution;
// }

// vector<vector<bool>> get_tabu_solution(const vector<vector<bool>>& solution, const vector<vector<int>>& problem_rows, const vector<vector<int>>& problem_cols, int iterations)
//{
//     vector<vector<bool>> best_solution = solution;
//     vector<vector<int>> cost_index_vec;
//     vector<vector<vector<bool>>> tabu;
//     int current_iteration = 0;
//     int best_cost = get_solution_cost(solution, problem_rows, problem_cols);
//     int current_cost = best_cost;
//     int lowest_cost_idx = 0;
//
//     while (best_cost > 0 && current_iteration <= iterations)
//     {
//         vector<vector<vector<bool>>> neighbors = generate_neighbors(best_solution);
//
//         for (int i = 0; i < neighbors.size(); i++)
//         {
//             if (!tabu.empty())
//             {
//                 if (find(tabu.begin(), tabu.end(), neighbors.at(i)) != tabu.end())
//                 {
//                     remove(neighbors.begin(), neighbors.end(), neighbors.at(i));
//                 }
//             }
//
//         }
//         //for (auto nei : neighbors)
//         //{
//         //    print_solution(nei);
//         //    cout << get_solution_cost(nei, problem_rows, problem_cols) << endl;;
//         //}
//         // Create vector of cost and solution indexes
//         for (int i = 0; i < neighbors.size(); i++)
//         {
//             current_cost = get_solution_cost(neighbors.at(i), problem_rows, problem_cols);
//             if (current_cost < best_cost)
//             {
//                 cost_index_vec.push_back({ current_cost, i });
//                 best_cost = current_cost;
//             }
//         }
//         // Find best solution
//         if (!cost_index_vec.empty())
//         {
//             int min_cost = INT_MAX;
//             int min_idx = INT_MAX;
//             for (int i = 0; i < cost_index_vec.size(); i++)
//             {
//                 if (cost_index_vec.at(i).at(0) < min_cost)
//                 {
//                     min_idx = cost_index_vec.at(i).at(1);
//                     min_cost = cost_index_vec.at(i).at(0);
//                 }
//             }
//             best_solution = neighbors.at(min_idx);
//         }
//         tabu.push_back(best_solution);
//         current_iteration++;
//     }
//     return best_solution;
// }

int main(int argc, char **argv)
{
    // get arguments
    auto fname_rows = arg(argc, argv, "fname_rows", string(""));
    auto fname_cols = arg(argc, argv, "fname_cols", string(""));
    auto iterations = arg(argc, argv, "iterations", 1000);
    // auto perfect_solution = arg(argc, argv, "solution", string(""));
    cout << "# fname_rows = " << fname_rows << ";" << endl;
    cout << "# fname_cols = " << fname_cols << ";" << endl;

    // Load the nonogram problem rows
    vector<vector<int>> problem_rows = load_problem(fname_rows);

    // Load the nonogram problem cols
    vector<vector<int>> problem_cols = load_problem(fname_cols);

    // vector<vector<bool>> valid_solution = load_valid_solution(perfect_solution);
    // cout << "Cost perfect: " << get_solution_cost(valid_solution, problem_rows, problem_cols) << endl;

    // get rozmiar of nonogram
    const int rows = problem_rows.size();
    const int cols = problem_cols.size();

    // Generate random solution test...
    // vector<vector<bool>> random_solution = generate_random_solution(problem_rows, rows, cols);
    vector<vector<bool>> random_solution = generate_completely_random_solution(problem_rows, rows, cols);

    // Generate problem from solution rows test...
    vector<vector<int>> problem_from_random_solution_rows = get_problem_rows_from_solution(random_solution);
    // vector<vector<int>> problem_from_random_solution_rows = get_problem_rows_from_solution(valid_solution);

    // Generate problem from solution rows test...
    vector<vector<int>> problem_from_random_solution_cols = get_problem_cols_from_solution(random_solution);
    // vector<vector<int>> problem_from_random_solution_cols = get_problem_cols_from_solution(valid_solution);

    // Testing problem contents...
    // cout << "Initial rows " << endl;;
    // print_problem(problem_rows);
    //
    // cout << "Initial cols " << endl;;
    // print_problem(problem_cols);
    //
    // cout << "From solution rows " << endl;
    // print_problem(problem_from_random_solution_rows);
    //
    // cout << "From solution cols " << endl;
    // print_problem(problem_from_random_solution_cols);

    // Testing random solution contents...
    // cout << "Random cost " << get_solution_cost(random_solution, problem_rows, problem_cols) << endl;
    print_solution(random_solution);
    // cout << "Problem rows: " << endl;
    // print_problem(get_problem_rows_from_solution(random_solution));
    // cout << "Problem cols: " << endl;
    // print_problem(get_problem_cols_from_solution(random_solution));

    // Testing the cost function
    cout << "Random solution cost: " << get_solution_cost(random_solution, problem_rows, problem_cols) << endl;

    // Testing the generate neighbors function
    // vector<vector<vector<bool>>> neighbors = generate_neighbors(random_solution);
    // for (int i = 0; i < neighbors.size(); i++)
    //{
    //    cout << "Neighbor " << i << " cost " << get_solution_cost(neighbors.at(i), problem_rows, problem_cols) << endl;
    //    print_solution(neighbors.at(i));
    //    cout << "Problem rows: " << endl;
    //    print_problem(get_problem_rows_from_solution(neighbors.at(i)));
    //    cout << "Problem cols: " << endl;
    //    print_problem(get_problem_cols_from_solution(neighbors.at(i)));
    //}

    // cout << "Random neighbor at 0 ";
    // print_solution(neighbors.at(0));

    // cout << "Random neighbor at 44 ";
    // print_solution(neighbors.at(44));
    //
    // cout << "Random neighbor at 87 ";
    // print_solution(neighbors.at(87));

    // Hillclimbing deterministic
    cout << "Hillclimbing Deterministic solution: " << endl;
    auto start = chrono::steady_clock::now();
    vector<vector<bool>> hillclimbing_solution = get_hillclimbing_solution(random_solution, problem_rows, problem_cols, iterations);
    auto end = chrono::steady_clock::now();
    print_solution(hillclimbing_solution);
    cout << "Cost: " << get_solution_cost(hillclimbing_solution, problem_rows, problem_cols) << endl;
    cout << "Execution time in miliseconds: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << endl;
    cout << "Iterations: " << iterations << endl;

    // Hillclimbing random
    cout << "Hillclimbing Random solution: " << endl;
    auto start2 = chrono::steady_clock::now();
    vector<vector<bool>> hillclimbing_random_solution = get_hillclimbing_random_solution(random_solution, problem_rows, problem_cols, iterations);
    auto end2 = chrono::steady_clock::now();
    print_solution(hillclimbing_random_solution);
    cout << "Cost: " << get_solution_cost(hillclimbing_random_solution, problem_rows, problem_cols) << endl;
    cout << "Execution time in miliseconds: " << chrono::duration_cast<chrono::milliseconds>(end2 - start2).count() << endl;
    cout << "Iterations: " << iterations << endl;

    // Tabu solution
    cout << "Tabu solution: " << endl;
    auto start3 = chrono::steady_clock::now();
    vector<vector<bool>> tabu_solution = get_tabu_solution(random_solution, problem_rows, problem_cols, iterations);
    auto end3 = chrono::steady_clock::now();
    print_solution(tabu_solution);
    cout << "Cost: " << get_solution_cost(tabu_solution, problem_rows, problem_cols) << endl;
    cout << "Execution time in miliseconds: " << chrono::duration_cast<chrono::milliseconds>(end3 - start3).count() << endl;
    cout << "Iterations: " << iterations << endl;

    //// Tabu search
    // cout << "Tabu solution: " << endl;
    // auto start3 = chrono::steady_clock::now();
    // vector<vector<bool>> tabu_solution = get_tabu_solution(random_solution, problem_rows, problem_cols, iterations);
    // auto end3 = chrono::steady_clock::now();
    // print_solution(tabu_solution);
    // cout << "Cost: " << get_solution_cost(tabu_solution, problem_rows, problem_cols) << endl;
    // cout << "Execution time in miliseconds: " << chrono::duration_cast<chrono::milliseconds>(end3 - start3).count() << endl;
    // cout << "Iterations: " << iterations << endl;
}
