#ifndef _CELLHODGKIN_HUXLEY_SQUID_AXON_MODEL_1952_MODIFIEDFROMCELLML_HPP__
#define _CELLHODGKIN_HUXLEY_SQUID_AXON_MODEL_1952_MODIFIEDFROMCELLML_HPP__

// Model: hodgkin_huxley_squid_axon_model_1952_modified
// Processed by pycml - CellML Tools in Python
//     (translators: , pycml: , optimize: )
// on Fri Mar 11 12:48:24 2016
#include <iostream>
#include <cmath>
#include <utility>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/iterator/constant_iterator.h>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>
using namespace std;
using namespace boost::numeric::odeint;

typedef double value_type;

typedef thrust::host_vector<float> state_type;

struct output_observer {
  value_type* cumulative_error_;
  output_observer(value_type *error) : cumulative_error_(error) { }
  void operator()(const state_type xs, const value_type t) {
    cout << t << ",";
    for ( int i = 0 ; i < xs.size(); i++) {
        cout << xs[i] << ",";
    }
    cout << endl;
  }
};

struct ode_system
{
    struct ode_functor
    {
        template< class T >
        __host__ __device__
        void operator() ( T t ) const
        {
            // Time units: millisecond
            value_type  var_environment__time = thrust::get< 0 > ( thrust::get< 0 > ( t ));
            // Inputs:
            // Units: millivolt; Initial value: -75
            const value_type var_membrane__V = thrust::get< 1 > ( thrust::get< 0 > ( t ));
            // Units: dimensionless; Initial value: 0.05
            const value_type var_sodium_channel_m_gate__m = thrust::get< 2 > ( thrust::get< 0 > ( t ));
            // Units: dimensionless; Initial value: 0.6
            const value_type var_sodium_channel_h_gate__h = thrust::get< 3 > ( thrust::get< 0 > ( t ));
            // Units: dimensionless; Initial value: 0.325
            const value_type var_potassium_channel_n_gate__n = thrust::get< 4 > ( thrust::get< 0 > ( t ));
            // Mathematics
            const value_type var_membrane__E_R =  -75.0; // millivolt
            const value_type var_membrane__Cm = 1.0; // microF_per_cm2
            const value_type var_membrane__time = var_environment__time; // millisecond
            const value_type var_sodium_channel__h = var_sodium_channel_h_gate__h; // dimensionless
            const value_type var_sodium_channel__m = var_sodium_channel_m_gate__m; // dimensionless
            const value_type var_sodium_channel__g_Na = 120.0; // milliS_per_cm2
            const value_type var_sodium_channel__E_R = var_membrane__E_R; // millivolt
            const value_type var_sodium_channel__E_Na = var_sodium_channel__E_R + 115.0; // millivolt
            const value_type var_sodium_channel__V = var_membrane__V; // millivolt
            const value_type var_sodium_channel__i_Na = var_sodium_channel__g_Na * pow(var_sodium_channel__m, 3.0) * var_sodium_channel__h * (var_sodium_channel__V - var_sodium_channel__E_Na); // microA_per_cm2
            const value_type var_membrane__i_Na = var_sodium_channel__i_Na; // microA_per_cm2
            const value_type var_potassium_channel__n = var_potassium_channel_n_gate__n; // dimensionless
            const value_type var_potassium_channel__E_R = var_membrane__E_R; // millivolt
            const value_type var_potassium_channel__E_K = var_potassium_channel__E_R - 12.0; // millivolt
            const value_type var_potassium_channel__g_K = 36.0; // milliS_per_cm2
            const value_type var_potassium_channel__V = var_membrane__V; // millivolt
            const value_type var_potassium_channel__i_K = var_potassium_channel__g_K * pow(var_potassium_channel__n, 4.0) * (var_potassium_channel__V - var_potassium_channel__E_K); // microA_per_cm2
            const value_type var_membrane__i_K = var_potassium_channel__i_K; // microA_per_cm2
            const value_type var_leakage_current__g_L = 0.29999999999999999; // milliS_per_cm2
            const value_type var_leakage_current__E_R = var_membrane__E_R; // millivolt
            const value_type var_leakage_current__E_L = var_leakage_current__E_R + 10.613; // millivolt
            const value_type var_leakage_current__V = var_membrane__V; // millivolt
            const value_type var_leakage_current__i_L = var_leakage_current__g_L * (var_leakage_current__V - var_leakage_current__E_L); // microA_per_cm2
            const value_type var_membrane__i_L = var_leakage_current__i_L; // microA_per_cm2
            const value_type var_membrane__stim_end = 10000.0; // millisecond
            const value_type var_membrane__stim_amplitude =  -20.0; // microA_per_cm2
            const value_type var_membrane__stim_duration = 0.5; // millisecond
            const value_type var_membrane__stim_period = 1000.0; // millisecond
            const value_type var_membrane__stim_start = 10.0; // millisecond
            const value_type var_membrane__i_Stim = ((var_membrane__time >= var_membrane__stim_start) && (var_membrane__time <= var_membrane__stim_end) && (((var_membrane__time - var_membrane__stim_start) - (floor((var_membrane__time - var_membrane__stim_start) / var_membrane__stim_period) * var_membrane__stim_period)) <= var_membrane__stim_duration)) ? var_membrane__stim_amplitude : 0.0; // microA_per_cm2
            const value_type var_sodium_channel__time = var_environment__time; // millisecond
            const value_type var_sodium_channel_m_gate__V = var_sodium_channel__V; // millivolt
            const value_type var_sodium_channel_m_gate__alpha_m = ((-0.10000000000000001) * (var_sodium_channel_m_gate__V + 50.0)) / (exp((-(var_sodium_channel_m_gate__V + 50.0)) / 10.0) - 1.0); // per_millisecond
            const value_type var_sodium_channel_m_gate__beta_m = 4.0 * exp((-(var_sodium_channel_m_gate__V + 75.0)) / 18.0); // per_millisecond
            const value_type var_sodium_channel_m_gate__time = var_sodium_channel__time; // millisecond
            const value_type var_sodium_channel_h_gate__V = var_sodium_channel__V; // millivolt
            const value_type var_sodium_channel_h_gate__alpha_h = 0.070000000000000007 * exp((-(var_sodium_channel_h_gate__V + 75.0)) / 20.0); // per_millisecond
            const value_type var_sodium_channel_h_gate__beta_h = 1.0 / (exp((-(var_sodium_channel_h_gate__V + 45.0)) / 10.0) + 1.0); // per_millisecond
            const value_type var_sodium_channel_h_gate__time = var_sodium_channel__time; // millisecond
            const value_type var_potassium_channel__time = var_environment__time; // millisecond
            const value_type var_potassium_channel_n_gate__V = var_potassium_channel__V; // millivolt
            const value_type var_potassium_channel_n_gate__alpha_n = ((-0.01) * (var_potassium_channel_n_gate__V + 65.0)) / (exp((-(var_potassium_channel_n_gate__V + 65.0)) / 10.0) - 1.0); // per_millisecond
            const value_type var_potassium_channel_n_gate__beta_n = 0.125 * exp((var_potassium_channel_n_gate__V + 75.0) / 80.0); // per_millisecond
            const value_type var_potassium_channel_n_gate__time = var_potassium_channel__time; // millisecond
            const value_type d_dt_membrane__V = (-(var_membrane__i_Stim + var_membrane__i_Na + var_membrane__i_K + var_membrane__i_L)) / var_membrane__Cm; // 'millivolt per millisecond'
            const value_type d_dt_sodium_channel_m_gate__m = (var_sodium_channel_m_gate__alpha_m * (1.0 - var_sodium_channel_m_gate__m)) - (var_sodium_channel_m_gate__beta_m * var_sodium_channel_m_gate__m); // per_millisecond
            const value_type d_dt_sodium_channel_h_gate__h = (var_sodium_channel_h_gate__alpha_h * (1.0 - var_sodium_channel_h_gate__h)) - (var_sodium_channel_h_gate__beta_h * var_sodium_channel_h_gate__h); // per_millisecond
            const value_type d_dt_potassium_channel_n_gate__n = (var_potassium_channel_n_gate__alpha_n * (1.0 - var_potassium_channel_n_gate__n)) - (var_potassium_channel_n_gate__beta_n * var_potassium_channel_n_gate__n); // per_millisecond
            thrust::get< 1 > ( thrust::get< 1 > ( t )) = d_dt_membrane__V;
            thrust::get< 2 > ( thrust::get< 1 > ( t )) = d_dt_sodium_channel_m_gate__m;
            thrust::get< 3 > ( thrust::get< 1 > ( t )) = d_dt_sodium_channel_h_gate__h;
            thrust::get< 4 > ( thrust::get< 1 > ( t )) = d_dt_potassium_channel_n_gate__n;
        }
    };
    
    size_t m_N;
    ode_system(size_t  N) : m_N( N ) { }
    template< class State, class Deriv >
    void operator() ( const State &ys, Deriv &dydts, value_type t) const
    {
        thrust::constant_iterator< value_type > t_iterator(t);
        thrust::for_each(
            thrust::make_zip_iterator(thrust::make_tuple(
                thrust::make_zip_iterator(thrust::make_tuple(
                    t_iterator + 0 * m_N,
                    boost::begin(ys) + 0 * m_N,
                    boost::begin(ys) + 1 * m_N,
                    boost::begin(ys) + 2 * m_N,
                    boost::begin(ys) + 3 * m_N))
                ,
                thrust::make_zip_iterator(thrust::make_tuple(
                    t_iterator + 0 * m_N,
                    boost::begin(dydts) + 0 * m_N,
                    boost::begin(dydts) + 1 * m_N,
                    boost::begin(dydts) + 2 * m_N,
                    boost::begin(dydts) + 3 * m_N))
            )),
            thrust::make_zip_iterator(thrust::make_tuple(
                thrust::make_zip_iterator(thrust::make_tuple(
                    t_iterator + 1 * m_N,
                    boost::begin(ys) + 1 * m_N,
                    boost::begin(ys) + 2 * m_N,
                    boost::begin(ys) + 3 * m_N,
                    boost::begin(ys) + 4 * m_N))
                ,
                thrust::make_zip_iterator(thrust::make_tuple(
                    t_iterator + 1 * m_N,
                    boost::begin(dydts) + 1 * m_N,
                    boost::begin(dydts) + 2 * m_N,
                    boost::begin(dydts) + 3 * m_N,
                    boost::begin(dydts) + 4 * m_N))
            )),
            ode_functor());
        }
    };
    
    int main()
    {
        const size_t N = 1;
        cout << "t" << ",";
        cout << "membrane_voltage0"  << ",";
        cout << "m0"  << ",";
        cout << "h0"  << ",";
        cout << "n0"  << ",";
        cout << endl;
        state_type ys(4);
        // Variable membrane_voltage
        thrust::fill(ys.begin() + 0 * N, ys.begin() + 1 * N, -75);
        // Variable m
        thrust::fill(ys.begin() + 1 * N, ys.begin() + 2 * N, 0.05);
        // Variable h
        thrust::fill(ys.begin() + 2 * N, ys.begin() + 3 * N, 0.6);
        // Variable n
        thrust::fill(ys.begin() + 3 * N, ys.begin() + 4 * N, 0.325);
        
        value_type error = 0.0;
        output_observer observer = output_observer(&error);
        
        ode_system system = ode_system(N);
        
        typedef runge_kutta_dopri5< state_type , value_type , state_type , value_type > stepper_type;
        integrate_adaptive(make_controlled(1E-12, 1E-12, stepper_type()),
                        system,
                        ys, 0.0, 1000.0, 0.01,
                        observer);
    }
    
    

    #endif
