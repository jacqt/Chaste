#ifndef _CELLHODGKIN_HUXLEY_SQUID_AXON_MODEL_1952_MODIFIEDFROMCELLML_HPP__
#define _CELLHODGKIN_HUXLEY_SQUID_AXON_MODEL_1952_MODIFIEDFROMCELLML_HPP__

// Model: hodgkin_huxley_squid_axon_model_1952_modified
// Processed by pycml - CellML Tools in Python
//     (translators: , pycml: , optimize: )
// on Fri Mar 11 10:03:02 2016
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
            value_type  var_environment__time = thrust::get< 0 > ( t );
            // Inputs:
            // Units: millivolt; Initial value: -75
            const value_type var_membrane__V = thrust::get< 1 > ( t );
            // Units: dimensionless; Initial value: 0.05
            const value_type var_sodium_channel_m_gate__m = thrust::get< 2 > ( t );
            // Units: dimensionless; Initial value: 0.6
            const value_type var_sodium_channel_h_gate__h = thrust::get< 3 > ( t );
            // Units: dimensionless; Initial value: 0.325
            const value_type var_potassium_channel_n_gate__n = thrust::get< 4 > ( t );
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
            thrust::get< 5 > ( t )  = d_dt_membrane__V;
            thrust::get< 6 > ( t )  = d_dt_sodium_channel_m_gate__m;
            thrust::get< 7 > ( t )  = d_dt_sodium_channel_h_gate__h;
            thrust::get< 8 > ( t )  = d_dt_potassium_channel_n_gate__n;
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
                t_iterator,
                boost::begin(ys) + 0 * m_N,
                boost::begin(ys) + 1 * m_N,
                boost::begin(ys) + 2 * m_N,
                boost::begin(ys) + 3 * m_N,
                boost::begin(dydts) + 0 * m_N,
                boost::begin(dydts) + 1 * m_N,
                boost::begin(dydts) + 2 * m_N,
                boost::begin(dydts) + 3 * m_N)),
            thrust::make_zip_iterator(thrust::make_tuple(
                t_iterator + m_N,
                boost::begin(ys) + 1 * m_N,
                boost::begin(ys) + 2 * m_N,
                boost::begin(ys) + 3 * m_N,
                boost::begin(ys) + 4 * m_N,
                boost::begin(dydts) + 1 * m_N,
                boost::begin(dydts) + 2 * m_N,
                boost::begin(dydts) + 3 * m_N,
                boost::begin(dydts) + 4 * m_N)),
            ode_functor());
    }
};

int main()
{
    const size_t N = 100;
    cout << "t" << ",";
    cout << "membrane_voltage0"  << ",";
    cout << "membrane_voltage1"  << ",";
    cout << "membrane_voltage2"  << ",";
    cout << "membrane_voltage3"  << ",";
    cout << "membrane_voltage4"  << ",";
    cout << "membrane_voltage5"  << ",";
    cout << "membrane_voltage6"  << ",";
    cout << "membrane_voltage7"  << ",";
    cout << "membrane_voltage8"  << ",";
    cout << "membrane_voltage9"  << ",";
    cout << "membrane_voltage10"  << ",";
    cout << "membrane_voltage11"  << ",";
    cout << "membrane_voltage12"  << ",";
    cout << "membrane_voltage13"  << ",";
    cout << "membrane_voltage14"  << ",";
    cout << "membrane_voltage15"  << ",";
    cout << "membrane_voltage16"  << ",";
    cout << "membrane_voltage17"  << ",";
    cout << "membrane_voltage18"  << ",";
    cout << "membrane_voltage19"  << ",";
    cout << "membrane_voltage20"  << ",";
    cout << "membrane_voltage21"  << ",";
    cout << "membrane_voltage22"  << ",";
    cout << "membrane_voltage23"  << ",";
    cout << "membrane_voltage24"  << ",";
    cout << "membrane_voltage25"  << ",";
    cout << "membrane_voltage26"  << ",";
    cout << "membrane_voltage27"  << ",";
    cout << "membrane_voltage28"  << ",";
    cout << "membrane_voltage29"  << ",";
    cout << "membrane_voltage30"  << ",";
    cout << "membrane_voltage31"  << ",";
    cout << "membrane_voltage32"  << ",";
    cout << "membrane_voltage33"  << ",";
    cout << "membrane_voltage34"  << ",";
    cout << "membrane_voltage35"  << ",";
    cout << "membrane_voltage36"  << ",";
    cout << "membrane_voltage37"  << ",";
    cout << "membrane_voltage38"  << ",";
    cout << "membrane_voltage39"  << ",";
    cout << "membrane_voltage40"  << ",";
    cout << "membrane_voltage41"  << ",";
    cout << "membrane_voltage42"  << ",";
    cout << "membrane_voltage43"  << ",";
    cout << "membrane_voltage44"  << ",";
    cout << "membrane_voltage45"  << ",";
    cout << "membrane_voltage46"  << ",";
    cout << "membrane_voltage47"  << ",";
    cout << "membrane_voltage48"  << ",";
    cout << "membrane_voltage49"  << ",";
    cout << "membrane_voltage50"  << ",";
    cout << "membrane_voltage51"  << ",";
    cout << "membrane_voltage52"  << ",";
    cout << "membrane_voltage53"  << ",";
    cout << "membrane_voltage54"  << ",";
    cout << "membrane_voltage55"  << ",";
    cout << "membrane_voltage56"  << ",";
    cout << "membrane_voltage57"  << ",";
    cout << "membrane_voltage58"  << ",";
    cout << "membrane_voltage59"  << ",";
    cout << "membrane_voltage60"  << ",";
    cout << "membrane_voltage61"  << ",";
    cout << "membrane_voltage62"  << ",";
    cout << "membrane_voltage63"  << ",";
    cout << "membrane_voltage64"  << ",";
    cout << "membrane_voltage65"  << ",";
    cout << "membrane_voltage66"  << ",";
    cout << "membrane_voltage67"  << ",";
    cout << "membrane_voltage68"  << ",";
    cout << "membrane_voltage69"  << ",";
    cout << "membrane_voltage70"  << ",";
    cout << "membrane_voltage71"  << ",";
    cout << "membrane_voltage72"  << ",";
    cout << "membrane_voltage73"  << ",";
    cout << "membrane_voltage74"  << ",";
    cout << "membrane_voltage75"  << ",";
    cout << "membrane_voltage76"  << ",";
    cout << "membrane_voltage77"  << ",";
    cout << "membrane_voltage78"  << ",";
    cout << "membrane_voltage79"  << ",";
    cout << "membrane_voltage80"  << ",";
    cout << "membrane_voltage81"  << ",";
    cout << "membrane_voltage82"  << ",";
    cout << "membrane_voltage83"  << ",";
    cout << "membrane_voltage84"  << ",";
    cout << "membrane_voltage85"  << ",";
    cout << "membrane_voltage86"  << ",";
    cout << "membrane_voltage87"  << ",";
    cout << "membrane_voltage88"  << ",";
    cout << "membrane_voltage89"  << ",";
    cout << "membrane_voltage90"  << ",";
    cout << "membrane_voltage91"  << ",";
    cout << "membrane_voltage92"  << ",";
    cout << "membrane_voltage93"  << ",";
    cout << "membrane_voltage94"  << ",";
    cout << "membrane_voltage95"  << ",";
    cout << "membrane_voltage96"  << ",";
    cout << "membrane_voltage97"  << ",";
    cout << "membrane_voltage98"  << ",";
    cout << "membrane_voltage99"  << ",";
    cout << "m0"  << ",";
    cout << "m1"  << ",";
    cout << "m2"  << ",";
    cout << "m3"  << ",";
    cout << "m4"  << ",";
    cout << "m5"  << ",";
    cout << "m6"  << ",";
    cout << "m7"  << ",";
    cout << "m8"  << ",";
    cout << "m9"  << ",";
    cout << "m10"  << ",";
    cout << "m11"  << ",";
    cout << "m12"  << ",";
    cout << "m13"  << ",";
    cout << "m14"  << ",";
    cout << "m15"  << ",";
    cout << "m16"  << ",";
    cout << "m17"  << ",";
    cout << "m18"  << ",";
    cout << "m19"  << ",";
    cout << "m20"  << ",";
    cout << "m21"  << ",";
    cout << "m22"  << ",";
    cout << "m23"  << ",";
    cout << "m24"  << ",";
    cout << "m25"  << ",";
    cout << "m26"  << ",";
    cout << "m27"  << ",";
    cout << "m28"  << ",";
    cout << "m29"  << ",";
    cout << "m30"  << ",";
    cout << "m31"  << ",";
    cout << "m32"  << ",";
    cout << "m33"  << ",";
    cout << "m34"  << ",";
    cout << "m35"  << ",";
    cout << "m36"  << ",";
    cout << "m37"  << ",";
    cout << "m38"  << ",";
    cout << "m39"  << ",";
    cout << "m40"  << ",";
    cout << "m41"  << ",";
    cout << "m42"  << ",";
    cout << "m43"  << ",";
    cout << "m44"  << ",";
    cout << "m45"  << ",";
    cout << "m46"  << ",";
    cout << "m47"  << ",";
    cout << "m48"  << ",";
    cout << "m49"  << ",";
    cout << "m50"  << ",";
    cout << "m51"  << ",";
    cout << "m52"  << ",";
    cout << "m53"  << ",";
    cout << "m54"  << ",";
    cout << "m55"  << ",";
    cout << "m56"  << ",";
    cout << "m57"  << ",";
    cout << "m58"  << ",";
    cout << "m59"  << ",";
    cout << "m60"  << ",";
    cout << "m61"  << ",";
    cout << "m62"  << ",";
    cout << "m63"  << ",";
    cout << "m64"  << ",";
    cout << "m65"  << ",";
    cout << "m66"  << ",";
    cout << "m67"  << ",";
    cout << "m68"  << ",";
    cout << "m69"  << ",";
    cout << "m70"  << ",";
    cout << "m71"  << ",";
    cout << "m72"  << ",";
    cout << "m73"  << ",";
    cout << "m74"  << ",";
    cout << "m75"  << ",";
    cout << "m76"  << ",";
    cout << "m77"  << ",";
    cout << "m78"  << ",";
    cout << "m79"  << ",";
    cout << "m80"  << ",";
    cout << "m81"  << ",";
    cout << "m82"  << ",";
    cout << "m83"  << ",";
    cout << "m84"  << ",";
    cout << "m85"  << ",";
    cout << "m86"  << ",";
    cout << "m87"  << ",";
    cout << "m88"  << ",";
    cout << "m89"  << ",";
    cout << "m90"  << ",";
    cout << "m91"  << ",";
    cout << "m92"  << ",";
    cout << "m93"  << ",";
    cout << "m94"  << ",";
    cout << "m95"  << ",";
    cout << "m96"  << ",";
    cout << "m97"  << ",";
    cout << "m98"  << ",";
    cout << "m99"  << ",";
    cout << "h0"  << ",";
    cout << "h1"  << ",";
    cout << "h2"  << ",";
    cout << "h3"  << ",";
    cout << "h4"  << ",";
    cout << "h5"  << ",";
    cout << "h6"  << ",";
    cout << "h7"  << ",";
    cout << "h8"  << ",";
    cout << "h9"  << ",";
    cout << "h10"  << ",";
    cout << "h11"  << ",";
    cout << "h12"  << ",";
    cout << "h13"  << ",";
    cout << "h14"  << ",";
    cout << "h15"  << ",";
    cout << "h16"  << ",";
    cout << "h17"  << ",";
    cout << "h18"  << ",";
    cout << "h19"  << ",";
    cout << "h20"  << ",";
    cout << "h21"  << ",";
    cout << "h22"  << ",";
    cout << "h23"  << ",";
    cout << "h24"  << ",";
    cout << "h25"  << ",";
    cout << "h26"  << ",";
    cout << "h27"  << ",";
    cout << "h28"  << ",";
    cout << "h29"  << ",";
    cout << "h30"  << ",";
    cout << "h31"  << ",";
    cout << "h32"  << ",";
    cout << "h33"  << ",";
    cout << "h34"  << ",";
    cout << "h35"  << ",";
    cout << "h36"  << ",";
    cout << "h37"  << ",";
    cout << "h38"  << ",";
    cout << "h39"  << ",";
    cout << "h40"  << ",";
    cout << "h41"  << ",";
    cout << "h42"  << ",";
    cout << "h43"  << ",";
    cout << "h44"  << ",";
    cout << "h45"  << ",";
    cout << "h46"  << ",";
    cout << "h47"  << ",";
    cout << "h48"  << ",";
    cout << "h49"  << ",";
    cout << "h50"  << ",";
    cout << "h51"  << ",";
    cout << "h52"  << ",";
    cout << "h53"  << ",";
    cout << "h54"  << ",";
    cout << "h55"  << ",";
    cout << "h56"  << ",";
    cout << "h57"  << ",";
    cout << "h58"  << ",";
    cout << "h59"  << ",";
    cout << "h60"  << ",";
    cout << "h61"  << ",";
    cout << "h62"  << ",";
    cout << "h63"  << ",";
    cout << "h64"  << ",";
    cout << "h65"  << ",";
    cout << "h66"  << ",";
    cout << "h67"  << ",";
    cout << "h68"  << ",";
    cout << "h69"  << ",";
    cout << "h70"  << ",";
    cout << "h71"  << ",";
    cout << "h72"  << ",";
    cout << "h73"  << ",";
    cout << "h74"  << ",";
    cout << "h75"  << ",";
    cout << "h76"  << ",";
    cout << "h77"  << ",";
    cout << "h78"  << ",";
    cout << "h79"  << ",";
    cout << "h80"  << ",";
    cout << "h81"  << ",";
    cout << "h82"  << ",";
    cout << "h83"  << ",";
    cout << "h84"  << ",";
    cout << "h85"  << ",";
    cout << "h86"  << ",";
    cout << "h87"  << ",";
    cout << "h88"  << ",";
    cout << "h89"  << ",";
    cout << "h90"  << ",";
    cout << "h91"  << ",";
    cout << "h92"  << ",";
    cout << "h93"  << ",";
    cout << "h94"  << ",";
    cout << "h95"  << ",";
    cout << "h96"  << ",";
    cout << "h97"  << ",";
    cout << "h98"  << ",";
    cout << "h99"  << ",";
    cout << "n0"  << ",";
    cout << "n1"  << ",";
    cout << "n2"  << ",";
    cout << "n3"  << ",";
    cout << "n4"  << ",";
    cout << "n5"  << ",";
    cout << "n6"  << ",";
    cout << "n7"  << ",";
    cout << "n8"  << ",";
    cout << "n9"  << ",";
    cout << "n10"  << ",";
    cout << "n11"  << ",";
    cout << "n12"  << ",";
    cout << "n13"  << ",";
    cout << "n14"  << ",";
    cout << "n15"  << ",";
    cout << "n16"  << ",";
    cout << "n17"  << ",";
    cout << "n18"  << ",";
    cout << "n19"  << ",";
    cout << "n20"  << ",";
    cout << "n21"  << ",";
    cout << "n22"  << ",";
    cout << "n23"  << ",";
    cout << "n24"  << ",";
    cout << "n25"  << ",";
    cout << "n26"  << ",";
    cout << "n27"  << ",";
    cout << "n28"  << ",";
    cout << "n29"  << ",";
    cout << "n30"  << ",";
    cout << "n31"  << ",";
    cout << "n32"  << ",";
    cout << "n33"  << ",";
    cout << "n34"  << ",";
    cout << "n35"  << ",";
    cout << "n36"  << ",";
    cout << "n37"  << ",";
    cout << "n38"  << ",";
    cout << "n39"  << ",";
    cout << "n40"  << ",";
    cout << "n41"  << ",";
    cout << "n42"  << ",";
    cout << "n43"  << ",";
    cout << "n44"  << ",";
    cout << "n45"  << ",";
    cout << "n46"  << ",";
    cout << "n47"  << ",";
    cout << "n48"  << ",";
    cout << "n49"  << ",";
    cout << "n50"  << ",";
    cout << "n51"  << ",";
    cout << "n52"  << ",";
    cout << "n53"  << ",";
    cout << "n54"  << ",";
    cout << "n55"  << ",";
    cout << "n56"  << ",";
    cout << "n57"  << ",";
    cout << "n58"  << ",";
    cout << "n59"  << ",";
    cout << "n60"  << ",";
    cout << "n61"  << ",";
    cout << "n62"  << ",";
    cout << "n63"  << ",";
    cout << "n64"  << ",";
    cout << "n65"  << ",";
    cout << "n66"  << ",";
    cout << "n67"  << ",";
    cout << "n68"  << ",";
    cout << "n69"  << ",";
    cout << "n70"  << ",";
    cout << "n71"  << ",";
    cout << "n72"  << ",";
    cout << "n73"  << ",";
    cout << "n74"  << ",";
    cout << "n75"  << ",";
    cout << "n76"  << ",";
    cout << "n77"  << ",";
    cout << "n78"  << ",";
    cout << "n79"  << ",";
    cout << "n80"  << ",";
    cout << "n81"  << ",";
    cout << "n82"  << ",";
    cout << "n83"  << ",";
    cout << "n84"  << ",";
    cout << "n85"  << ",";
    cout << "n86"  << ",";
    cout << "n87"  << ",";
    cout << "n88"  << ",";
    cout << "n89"  << ",";
    cout << "n90"  << ",";
    cout << "n91"  << ",";
    cout << "n92"  << ",";
    cout << "n93"  << ",";
    cout << "n94"  << ",";
    cout << "n95"  << ",";
    cout << "n96"  << ",";
    cout << "n97"  << ",";
    cout << "n98"  << ",";
    cout << "n99"  << ",";
    cout << endl;
    state_type ys(400);
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
    
    integrate_const(runge_kutta4< state_type >(), 
                    system,
                    ys, 0.0, 50.0, 0.01,
                    observer);
}



#endif
