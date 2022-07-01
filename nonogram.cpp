#include <iostream>
#include <string>
#include <any>
#include <vector>
#include <functional>
#include <fstream>
#include <sstream>
#include <random>

using namespace std;

random_device rd;
mt19937 random_generator(rd());

auto arg = [](int argc, char** argv, std::string name, auto default_value) -> decltype(default_value)
{
    //using namespace std;
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

vector<vector<int>> load_problem(string fname)
{
    cout << "Loading problem..." << endl;
    vector<vector<int>> problem;
    istream* fst = &cin;
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
                //cout << token << endl; // Test output
                nums.push_back(stoi(token));
                text.erase(0, pos + delimiter.length());
            }
            //cout << "TEXT: " << text << endl; // Test output
            nums.push_back(stoi(text));
            problem.push_back(nums);
        }
    }
    return problem;
}

vector<vector<int>> load_valid_solution(string fname)
{
    cout << "Loading valid solution..." << endl;
    vector<vector<int>> valid_solution;
    istream* fst = &cin;
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
        valid_solution.push_back(row);
    }
    return valid_solution;
}

vector<int> mark_ones(vector<int> vec)
{
    return {};
}

vector<vector<int>> generate_random_solution(vector<vector<int>> problem, vector<vector<int>> valid_solution)
{
    vector<vector<int>> random_solution;
    int rows = static_cast<int>(valid_solution.size());
    int cols = static_cast<int>(valid_solution.at(0).size());
   

    // Initialize random_solution vector
    for (int i = 0; i < valid_solution.size(); i++)
    {
        vector<int> col(valid_solution.at(i).size(), 0);
        random_solution.push_back(col);
    }

    // Mark ones for the rows
    // Looping valid solution rows
    // i = problem vector   AND     i = valid_solution vector
    for (int i = 0; i < rows; i++)  // 0-10
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
        if (--range < 0) range = 0;

        uniform_int_distribution<int> distr(0, range);
        offset = distr(random_generator);
        cout << "Offset: " << offset << endl;

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
                    random_solution.at(i).at(idx+offset) = 1; // later l + offset
                    idx++;
                }
                idx++;
            }
            
        }
    }

    return random_solution;
}

void print_problem(vector<vector<int>> problem)
{
    cout << "Problem : " << endl;
    for (auto vec : problem)
    {
        cout << "Vector: ";
        for (auto num : vec)
        {
            cout << num << " ";
        }
        cout << endl;
    }
}

void print_solution(vector<vector<int>> solution)
{
    cout << "Solution: " << endl;
    for (auto vec : solution)
    {
        for (auto num : vec)
        {
            cout << num << " ";
        }
        cout << endl;
    }
}

int main(int argc, char **argv)
{
    // get arguments
    auto fname = arg(argc, argv, "fname", string(""));
    auto vsol = arg(argc, argv, "vsol", string(""));
    cout << "# fname = " << fname << ";" << endl;
    cout << "# vsol = " << vsol << ";" << endl;

    // Load the nonogram problem
    vector<vector<int>> problem = load_problem(fname);

    // Load the valid solution
    vector<vector<int>> valid_solution = load_valid_solution(vsol);

    // Generate random solution test...
    vector<vector<int>> random_solution = generate_random_solution(problem, valid_solution);

    // Testing problem contents...
    print_problem(problem);

    // Testing valid solution contents...
    print_solution(valid_solution);

    // Testing random solution contents...
    cout << "Random ";
    print_solution(random_solution);
}
 