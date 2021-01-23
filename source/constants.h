//
// Created by ilaria on 2021-01-23.
//
#pragma once

/*Number of components*/
inline constexpr int DIM = 3;
inline constexpr int NC = 2;
extern int Lx, Ly, Lz, N;

#define C_TWO_PI (6.2831853071795864769252867665590058L)
#define C_PI (3.1415926535897932384626433832795029L)
/*MaxP corresponds to the number of q states of the theta field. Division of phase values in <-pi, pi>. Must have MaxP uneven to have a symmetric interval around zero.*/
//static constexpr int MaxP=129;
#define MaxP (129)

#undef  SQR
#define SQR(x) ((x)*(x))

namespace settings {

#ifdef NDEBUG
    inline constexpr bool debug = false;
#else
    inline constexpr bool debug = true;
#endif

}