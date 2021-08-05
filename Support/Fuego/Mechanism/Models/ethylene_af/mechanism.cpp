#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"



/*save atomic weights into array */
void atomicWeight(amrex::Real *  awt)
{
    awt[0] = 14.006700; /*N */
    awt[1] = 1.007970; /*H */
    awt[2] = 15.999400; /*O */
    awt[3] = 12.011150; /*C */

    return;
}



/*get atomic weight for all elements */
void CKAWT( amrex::Real *  awt)
{
    atomicWeight(awt);
}



/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * ncf)
{
    int id; /*loop counter */
    int kd = 4; 
    /*Zero ncf */
    for (id = 0; id < kd * 29; ++ id) {
         ncf[id] = 0; 
    }

    /*N2 */
    ncf[ 0 * kd + 0 ] = 2; /*N */

    /*H */
    ncf[ 1 * kd + 1 ] = 1; /*H */

    /*H2 */
    ncf[ 2 * kd + 1 ] = 2; /*H */

    /*O */
    ncf[ 3 * kd + 2 ] = 1; /*O */

    /*OH */
    ncf[ 4 * kd + 2 ] = 1; /*O */
    ncf[ 4 * kd + 1 ] = 1; /*H */

    /*O2 */
    ncf[ 5 * kd + 2 ] = 2; /*O */

    /*H2O2 */
    ncf[ 6 * kd + 1 ] = 2; /*H */
    ncf[ 6 * kd + 2 ] = 2; /*O */

    /*H2O */
    ncf[ 7 * kd + 1 ] = 2; /*H */
    ncf[ 7 * kd + 2 ] = 1; /*O */

    /*HO2 */
    ncf[ 8 * kd + 1 ] = 1; /*H */
    ncf[ 8 * kd + 2 ] = 2; /*O */

    /*CO */
    ncf[ 9 * kd + 3 ] = 1; /*C */
    ncf[ 9 * kd + 2 ] = 1; /*O */

    /*CH3 */
    ncf[ 10 * kd + 3 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 3; /*H */

    /*CH2O */
    ncf[ 11 * kd + 1 ] = 2; /*H */
    ncf[ 11 * kd + 3 ] = 1; /*C */
    ncf[ 11 * kd + 2 ] = 1; /*O */

    /*CO2 */
    ncf[ 12 * kd + 3 ] = 1; /*C */
    ncf[ 12 * kd + 2 ] = 2; /*O */

    /*CH4 */
    ncf[ 13 * kd + 3 ] = 1; /*C */
    ncf[ 13 * kd + 1 ] = 4; /*H */

    /*C2H2 */
    ncf[ 14 * kd + 3 ] = 2; /*C */
    ncf[ 14 * kd + 1 ] = 2; /*H */

    /*C2H4 */
    ncf[ 15 * kd + 3 ] = 2; /*C */
    ncf[ 15 * kd + 1 ] = 4; /*H */

    /*CH2CO */
    ncf[ 16 * kd + 3 ] = 2; /*C */
    ncf[ 16 * kd + 1 ] = 2; /*H */
    ncf[ 16 * kd + 2 ] = 1; /*O */

    /*C2H6 */
    ncf[ 17 * kd + 3 ] = 2; /*C */
    ncf[ 17 * kd + 1 ] = 6; /*H */

    /*C */
    ncf[ 18 * kd + 3 ] = 1; /*C */

    /*CH */
    ncf[ 19 * kd + 3 ] = 1; /*C */
    ncf[ 19 * kd + 1 ] = 1; /*H */

    /*HCO */
    ncf[ 20 * kd + 1 ] = 1; /*H */
    ncf[ 20 * kd + 3 ] = 1; /*C */
    ncf[ 20 * kd + 2 ] = 1; /*O */

    /*TXCH2 */
    ncf[ 21 * kd + 3 ] = 1; /*C */
    ncf[ 21 * kd + 1 ] = 2; /*H */

    /*SXCH2 */
    ncf[ 22 * kd + 3 ] = 1; /*C */
    ncf[ 22 * kd + 1 ] = 2; /*H */

    /*C2H3 */
    ncf[ 23 * kd + 3 ] = 2; /*C */
    ncf[ 23 * kd + 1 ] = 3; /*H */

    /*C2H5 */
    ncf[ 24 * kd + 3 ] = 2; /*C */
    ncf[ 24 * kd + 1 ] = 5; /*H */

    /*HCCO */
    ncf[ 25 * kd + 1 ] = 1; /*H */
    ncf[ 25 * kd + 3 ] = 2; /*C */
    ncf[ 25 * kd + 2 ] = 1; /*O */

    /*CH3CHO */
    ncf[ 26 * kd + 1 ] = 4; /*H */
    ncf[ 26 * kd + 3 ] = 2; /*C */
    ncf[ 26 * kd + 2 ] = 1; /*O */

    /*CH2CHO */
    ncf[ 27 * kd + 1 ] = 3; /*H */
    ncf[ 27 * kd + 3 ] = 2; /*C */
    ncf[ 27 * kd + 2 ] = 1; /*O */

    /*C2H5O */
    ncf[ 28 * kd + 3 ] = 2; /*C */
    ncf[ 28 * kd + 1 ] = 5; /*H */
    ncf[ 28 * kd + 2 ] = 1; /*O */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(4);
    ename[0] = "N";
    ename[1] = "H";
    ename[2] = "O";
    ename[3] = "C";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(29);
    kname[0] = "N2";
    kname[1] = "H";
    kname[2] = "H2";
    kname[3] = "O";
    kname[4] = "OH";
    kname[5] = "O2";
    kname[6] = "H2O2";
    kname[7] = "H2O";
    kname[8] = "HO2";
    kname[9] = "CO";
    kname[10] = "CH3";
    kname[11] = "CH2O";
    kname[12] = "CO2";
    kname[13] = "CH4";
    kname[14] = "C2H2";
    kname[15] = "C2H4";
    kname[16] = "CH2CO";
    kname[17] = "C2H6";
    kname[18] = "C";
    kname[19] = "CH";
    kname[20] = "HCO";
    kname[21] = "TXCH2";
    kname[22] = "SXCH2";
    kname[23] = "C2H3";
    kname[24] = "C2H5";
    kname[25] = "HCCO";
    kname[26] = "CH3CHO";
    kname[27] = "CH2CHO";
    kname[28] = "C2H5O";
}

/*compute the sparsity pattern of the chemistry Jacobian */
#ifdef COMPILE_JACOBIAN
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<29; l++) {
                c_d[l] = 1.0/ 29.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<30; k++) {
        for (int l=0; l<30; l++) {
            if(J_h[ 30 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;

    return;
}
#endif



/*compute the sparsity pattern of the system Jacobian */
#ifdef COMPILE_JACOBIAN
void SPARSITY_INFO_SYST( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<29; k++) {
                c_d[k] = 1.0/ 29.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<30; k++) {
        for (int l=0; l<30; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 30 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;

    return;
}
#endif



/*compute the sparsity pattern of the simplified (for preconditioning) system Jacobian */
#ifdef COMPILE_JACOBIAN
void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, int * consP)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<29; k++) {
                c_d[k] = 1.0/ 29.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<30; k++) {
        for (int l=0; l<30; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 30 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    nJdata[0] = nJdata_tmp;

    return;
}
#endif


/*compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0 */
#ifdef COMPILE_JACOBIAN
void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, int * consP, int NCELLS)
{
    int offset_row;
    int offset_col;

    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<29; k++) {
                c_d[k] = 1.0/ 29.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        offset_row = nc * 30;
        offset_col = nc * 30;
        for (int k=0; k<30; k++) {
            for (int l=0; l<30; l++) {
                if(J_h[30*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l + offset_row; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
            colPtrs[offset_col + (k + 1)] = nJdata_tmp;
        }
    }

    return;
}
#endif

/*compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0 */
#ifdef COMPILE_JACOBIAN
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base)
{
    int offset;
    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<29; k++) {
                c_d[k] = 1.0/ 29.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    if (base == 1) {
        rowPtrs[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            offset = nc * 30;
            for (int l=0; l<30; l++) {
                for (int k=0; k<30; k++) {
                    if(J_h[30*k + l] != 0.0) {
                        colVals[nJdata_tmp-1] = k+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
                rowPtrs[offset + (l + 1)] = nJdata_tmp;
            }
        }
    } else {
        rowPtrs[0] = 0;
        int nJdata_tmp = 0;
        for (int nc=0; nc<NCELLS; nc++) {
            offset = nc * 30;
            for (int l=0; l<30; l++) {
                for (int k=0; k<30; k++) {
                    if(J_h[30*k + l] != 0.0) {
                        colVals[nJdata_tmp] = k + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
                rowPtrs[offset + (l + 1)] = nJdata_tmp;
            }
        }
    }

    return;
}
#endif

/*compute the sparsity pattern of the system Jacobian */
/*CSR format BASE is user choice */
#ifdef COMPILE_JACOBIAN
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, int * consP, int NCELLS, int base)
{
    int offset;
    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<29; k++) {
                c_d[k] = 1.0/ 29.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif
    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            offset = nc * 30;
            for (int l=0; l<30; l++) {
                for (int k=0; k<30; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[30*k + l] != 0.0) {
                            colVals[nJdata_tmp-1] = k+1 + offset; 
                            nJdata_tmp = nJdata_tmp + 1; 
                        }
                    }
                }
                rowPtr[offset + (l + 1)] = nJdata_tmp;
            }
        }
    } else {
        rowPtr[0] = 0;
        int nJdata_tmp = 0;
        for (int nc=0; nc<NCELLS; nc++) {
            offset = nc * 30;
            for (int l=0; l<30; l++) {
                for (int k=0; k<30; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[30*k + l] != 0.0) {
                            colVals[nJdata_tmp] = k + offset; 
                            nJdata_tmp = nJdata_tmp + 1; 
                        }
                    }
                }
                rowPtr[offset + (l + 1)] = nJdata_tmp;
            }
        }
    }

    return;
}
#endif

/*compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU */
/*BASE 0 */
#ifdef COMPILE_JACOBIAN
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, int * consP)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<29; k++) {
                c_d[k] = 1.0/ 29.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<30; k++) {
        for (int l=0; l<30; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 30*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[30*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 30*k + l;
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
        }
        colPtrs[k+1] = nJdata_tmp;
    }

    return;
}
#endif

/*compute the sparsity pattern of the simplified (for precond) system Jacobian */
/*CSR format BASE is under choice */
#ifdef COMPILE_JACOBIAN
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<29; k++) {
                c_d[k] = 1.0/ 29.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int l=0; l<30; l++) {
            for (int k=0; k<30; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[30*k + l] != 0.0) {
                        colVals[nJdata_tmp-1] = k+1; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
            }
            rowPtr[l+1] = nJdata_tmp;
        }
    } else {
        rowPtr[0] = 0;
        int nJdata_tmp = 0;
        for (int l=0; l<30; l++) {
            for (int k=0; k<30; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[30*k + l] != 0.0) {
                        colVals[nJdata_tmp] = k; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
            }
            rowPtr[l+1] = nJdata_tmp;
        }
    }

    return;
}
#endif

#endif
