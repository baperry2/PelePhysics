#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

// Base class to do something pointless
class MultiplyByNumber {
public: 
  AMREX_FORCE_INLINE
  AMREX_GPU_HOST_DEVICE
  virtual amrex::Real multiply (amrex::Real in) = 0;
};

// Derived class, final function
class MultiplyByTwo : public MultiplyByNumber {
public:
  AMREX_FORCE_INLINE
  AMREX_GPU_HOST_DEVICE
  amrex::Real multiply (amrex::Real in) final {
    return in *	2.0;
  };
};

// Derived class, override function
class MultiplyByThree : public MultiplyByNumber {
public:
  AMREX_FORCE_INLINE
  AMREX_GPU_HOST_DEVICE
  amrex::Real multiply (amrex::Real in) override {
    return in *	3.0;
  };
};

// Wrapper around class, templated for convenience
template <typename MultiplyType>
class MultiplyWrapper {
public:
  AMREX_FORCE_INLINE
  AMREX_GPU_HOST_DEVICE
  MultiplyType* get_mult() {
    return(&mult);
  }
private:
  MultiplyType mult;
};

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);

  // initialize some dummy data
  int npts = 32;
  amrex::Box domain(amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
		    amrex::IntVect(AMREX_D_DECL(npts - 1, npts - 1, npts - 1)));
  amrex::BoxArray ba(domain);
  amrex::DistributionMapping dm(ba);
  int nvar{1};
  int nghost{0};
  amrex::MultiFab mf(ba, dm, nvar, nghost);
  mf.setVal(10.0);
  auto const& mf_arrs = mf.arrays();

  amrex::Print() << "Max val is " << mf.max(0) << std::endl;
  amrex::Print() << "Directly call MultiplyByThree" << std::endl;
  // Loop over data, use virtual function without polymorphism
  amrex::ParallelFor(mf, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
    MultiplyByThree mult; 
    mf_arrs[nbx](i,j,k) = mult.multiply(mf_arrs[nbx](i,j,k));  
  });
  amrex::Print() << "Max val is " << mf.max(0) << std::endl;  

  amrex::Print() << "Call MultiplyByTwo indirectly (function is marked `final`)" << std::endl;
  // Loop over data, use final virtual function in wrapped class
  amrex::ParallelFor(mf, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
    MultiplyWrapper<MultiplyByTwo> multwrap;
    mf_arrs[nbx](i,j,k) = multwrap.get_mult()->multiply(mf_arrs[nbx](i,j,k));
  });
  amrex::Print() << "Max val is " << mf.max(0) << std::endl;

  amrex::Print() << "Call MultiplyByThree indirectly (function is marked `override`)..." << std::endl;
  amrex::Print() << "  code will hang, unless compiled with `-fsycl-enable-function-pointers`, in which case it will segfault" << std::endl;
  // Loop over data, use overriden virtual function in wrapped class
  amrex::ParallelFor(mf, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
    MultiplyWrapper<MultiplyByThree> multwrap;
    mf_arrs[nbx](i,j,k) = multwrap.get_mult()->multiply(mf_arrs[nbx](i,j,k));
  });
  amrex::Print() << "Max val is " << mf.max(0) << std::endl;      
    
  amrex::Finalize();
  return 0;
}
