#include "TimeIntegrator.h"

TimeIntegrator::TimeIntegrator(Parameters *param_, Ode *ode_)
    : param(param_)
{
    ode = ode_;
    n_output = 1;
}

TimeIntegrator::~TimeIntegrator()
{
}

Real TimeIntegrator::get_cpu_time()
{
    return (Real)clock() / CLOCKS_PER_SEC;
}

void TimeIntegrator::output_solution(Real t, Vector *y)
{
    // First we write the integrating variable in a matlab readable file
    unsigned int N = ode->get_system_size();
    ofstream outfile;

    if (param->matlab_output)
    {
        if (n_output == 1)
            outfile.open(param->output_path + string("_evolution.m"), ofstream::out);
        else
            outfile.open(param->output_path + string("_evolution.m"), ofstream::out | ofstream::app);
        outfile << setprecision(16) << "t(" << n_output << ")=" << t << ";" << endl;
        outfile << "y(" << n_output << ",:)=[";
        for (int i = 0; i < N - 1; i++)
            outfile << (*y)(i) << ",";
        outfile << (*y)(N - 1) << "];" << endl;
        outfile.close();
    }

    // Then we write it in a binary file
    if (param->bin_output)
    {
        if (n_output == 1)
            outfile.open(param->output_path + string("_evolution.bin"), ios::out | ios::binary);
        else
            outfile.open(param->output_path + string("_evolution.bin"), ios::out | ios::binary | ios::app);
        outfile.write((char *)&t, sizeof(double));
        outfile.write((char *)&(*y)(0), N * sizeof(double));
        outfile.close();
    }

    n_output++;
}

void TimeIntegrator::output_final_solution(Vector *y)
{
    unsigned int N = ode->get_system_size();
    ofstream outfile;

    outfile.open(param->output_path + string(".csv"), ofstream::out);
    outfile << "y" << endl;
    outfile << setprecision(16);
    for (int i = 0; i < N - 1; i++)
        outfile << (*y)(i) << endl;
    outfile << (*y)(N - 1);
    outfile.close();

    outfile.open(param->output_path + string(".bin"), ios::out | ios::binary);
    outfile.write((char *)&(*y)(0), N * sizeof(double));
    outfile.close();

    outfile.open(param->output_path + string("_statistics.csv"), ofstream::out);
    outfile << "cpu_time" << endl
            << elapsed_time;
    outfile.close();
}
