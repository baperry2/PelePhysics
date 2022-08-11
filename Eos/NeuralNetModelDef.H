#ifndef _NEURAL_NET_MD_H_
#define _NEURAL_NET_MD_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <memory>

#include "NeuralNetLayerDef.H"

namespace pele {
namespace physics {

AMREX_FORCE_INLINE
void read_fmt(const std::string& fmt, std::vector<char>& types,
    std::vector<std::vector<int>>& sizes)
{
  constexpr int READING_DTYPE = 0;
  constexpr int READING_PAREN = 1;
  constexpr int READING_SIZE  = 2;
  constexpr int READING_DELIM = 3;
  constexpr int TRAP_STATE = 4;
  
  int state = READING_DTYPE;
  
  std::istringstream is(fmt);
  int char_idx = -1;
  
  while(is.peek() != EOF && is.peek() != '\0')
  {
    switch(state)
    {
      case READING_DTYPE:
        char c0; is >> c0;
        types.push_back(c0);
        state = READING_PAREN;
        break;
      case READING_PAREN:
        char c1; is >> c1;
        if(c1 != '(')
        {
          state = TRAP_STATE;
        }
        else
        {
          state = READING_SIZE;  
        }
        break;
      case READING_SIZE:
        int i; is >> i;
        if(sizes.size() < types.size())
        {
          sizes.push_back(std::vector<int>());
        }
        sizes[sizes.size()-1].push_back(i);
        state = READING_DELIM;
        break;
      case READING_DELIM:
        char c3; is >> c3;
        if(c3 == ',')
        {
          state = READING_SIZE;
        }
        else if(c3 == ')')
        {
          state = READING_DTYPE;
        }
        else
        {
          state = TRAP_STATE;
        }
        break;
      default:
        std::string msg = "Error parsing neural net layer format. At character '";
        msg.append(std::string(1, fmt[char_idx]));
        msg.append("' and index " + std::to_string(char_idx) + ".");
        amrex::Abort(msg);
    }
    
    char_idx++;
  }
}

struct InitParams
{
  int ssize;
  int isize;
  int fsize;
  int inpsize;
  int nnets;
  bool cmlm_net = false;
};

class NNModelBase
{
public:
  
  AMREX_GPU_HOST_DEVICE
  virtual ~NNModelBase() {}
  
  virtual void init_from_stream(std::ifstream& file, const InitParams& ipar) = 0;
  
  AMREX_GPU_HOST_DEVICE
  virtual void init_from_buffers(amrex::Real* nnrd, int* nnid, const InitParams& ipar) = 0;
  
  // Evaluate net for given inputs and return pointer to output buffer
  AMREX_GPU_HOST_DEVICE
  virtual const amrex::Real* operator()(const amrex::Real* inputs) = 0;
  
  // Evaluate net and return pointer to the net output at the given index
  // Does not guarantee that any other outputs were calculated
  AMREX_GPU_HOST_DEVICE
  virtual const amrex::Real* operator()(const amrex::Real* inputs, int idx) = 0;
  
  // Evaluate net and return pointer to the net output at istart
  // Does not guarantee that outputs outside of the [istart, iend] range were calculated
  AMREX_GPU_HOST_DEVICE
  virtual const amrex::Real* operator()(const amrex::Real* inputs, int istart, int iend) = 0;
  
  // Evaluate net and store the results in outputs
  AMREX_GPU_HOST_DEVICE
  virtual void operator()(const amrex::Real* inputs, amrex::Real* outputs) = 0;
  
  virtual void print() = 0;
  
  AMREX_GPU_HOST_DEVICE
  virtual int nin() = 0;
  
  AMREX_GPU_HOST_DEVICE
  virtual int nout() = 0;
  
  AMREX_GPU_HOST_DEVICE
  virtual int nfparams() = 0;
  
  AMREX_GPU_HOST_DEVICE
  virtual int niparams(bool with_prefix=true) = 0;

  virtual void pack(amrex::Real* nnrd, int* nnid, bool store_prefix=true) = 0;
};

class NNModel: public NNModelBase
{
public:
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  NNModel() {}
  
  AMREX_GPU_HOST_DEVICE
  ~NNModel()
  {
    for(int i = 0; i < num_layers; i++)
    {
      delete layers[i];
    }
    delete[] layers;
    delete[] lsizes;
    delete[] outvec0;
    delete[] outvec1;
  }
  
  virtual void init_from_stream(std::ifstream& file, const InitParams& ipar)
  {
    cmlm_net = ipar.cmlm_net;
    inpsize = ipar.inpsize;
    
    int ssize, fsize, isize;
    ssize = ipar.ssize;
    fsize = ipar.fsize;
    isize = ipar.isize;
    
    file.read(reinterpret_cast<char*>(&num_layers), sizeof(int));
    
    layers = new NNLayer*[num_layers];
    
    for(int i = 0; i < num_layers; i++)
    {
      char name_arr[ssize];
      char fmt_arr[ssize];
      file.read(&name_arr[0], ssize*sizeof(char));
      file.read(&fmt_arr[0],  ssize*sizeof(char));
      std::string name = amrex::trim(std::string(&name_arr[0], ssize), std::string(1, '\0'));
      std::string fmt  = amrex::trim(std::string(&fmt_arr[0], ssize), std::string(1, '\0'));
      
      std::vector<char> types;
      std::vector<std::vector<int>> sizes;
      read_fmt(fmt, types, sizes);
      
      if(name == "Linear")
      {
        layers[i] = new LinearNNLayer(types, sizes, file, fsize, isize);
      }
      else if(name == "LeakyReLU")
      {
        layers[i] = new LeakyReluNNLayer(types, sizes, file, fsize, isize);
      }
      else if(name == "BatchNorm1d")
      {
        layers[i] = new BatchNorm1dNNLayer(types, sizes, file, fsize, isize);
      }
    }
    
    lsizes = new int[num_layers];
    lsizes[0] = (layers[0]->nout == -1) ? inpsize : layers[0]->nout;
    for(int i = 1; i < num_layers; i++)
    {
      lsizes[i] = (layers[i]->nout == -1) ? lsizes[i-1] : layers[i]->nout;
    }
    
    int maxsize = 0;
    for(int i = 0; i < num_layers; i++)
    {
      layers[i]->nout = lsizes[i];
      maxsize = (lsizes[i] > maxsize) ? lsizes[i] : maxsize;
    }
    
    outvec0 = new amrex::Real[maxsize];
    outvec1 = new amrex::Real[maxsize];
    
    // Do not need to store Real parameters for net, just the layers
    nreals = 0;
    // We want to store number of nets, `cmlm_net`, `inpsize`, `nreals`, `nints`
    // Edit `niparams` and `pack` functions if this changes
    nints = 5;
    
    for(int i = 0; i < num_layers; i++)
    {
      nreals += layers[i]->nfparams();
      nints  += layers[i]->niparams();
    }
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void init_from_buffers(amrex::Real* nnrd, int* nnid, const InitParams& ipar)
  {    
    cmlm_net = ipar.cmlm_net;
    inpsize = ipar.inpsize;
    
    nreals = nnid[0];
    nints = nnid[1];
    nnid += 2;
    
    if(cmlm_net)
    {
      num_layers = 7;
      layers = new NNLayer*[num_layers];
      lsizes = new int[num_layers];
      
      layers[0] = new LinearNNLayer(nnrd, nnid);
      nnrd += layers[0]->nfparams();
      nnid += layers[0]->niparams();
      
      layers[1] = new LeakyReluNNLayer(nnrd, nnid);
      nnrd += layers[1]->nfparams();
      nnid += layers[1]->niparams();
      
      layers[2] = new BatchNorm1dNNLayer(nnrd, nnid);
      nnrd += layers[2]->nfparams();
      nnid += layers[2]->niparams();
      
      layers[3] = new LinearNNLayer(nnrd, nnid);
      nnrd += layers[3]->nfparams();
      nnid += layers[3]->niparams();
      
      layers[4] = new LeakyReluNNLayer(nnrd, nnid);
      nnrd += layers[4]->nfparams();
      nnid += layers[4]->niparams();
      
      layers[5] = new BatchNorm1dNNLayer(nnrd, nnid);
      nnrd += layers[5]->nfparams();
      nnid += layers[5]->niparams();
      
      layers[6] = new LinearNNLayer(nnrd, nnid);
      
      lsizes[0] = (layers[0]->nout == -1) ? inpsize : layers[0]->nout;
      for(int i = 1; i < num_layers; i++)
      {
        lsizes[i] = (layers[i]->nout == -1) ? lsizes[i-1] : layers[i]->nout;
      }
      
      int maxsize = 0;
      for(int i = 0; i < num_layers; i++)
      {
        layers[i]->nout = lsizes[i];
        maxsize = (lsizes[i] > maxsize) ? lsizes[i] : maxsize;
      }
      
      outvec0 = new amrex::Real[maxsize];
      outvec1 = new amrex::Real[maxsize];
    }
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual const amrex::Real* operator()(const amrex::Real* inputs)
  {
    amrex::Real* invec = outvec0;
    amrex::Real* outvec = outvec1;
    amrex::Real* temp;
    int nin = inpsize;
  
    (*layers[0])(inputs, outvec, nin);
    nin = lsizes[0];
    
    for(int i = 1; i < num_layers; i++)
    {
      // Swap input and output buffers
      temp = invec;
      invec = outvec;
      outvec = temp;
      // Evaluate layer
      (*layers[i])(invec, outvec, nin);
      nin = lsizes[i];
    }
    
    return outvec;
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual const amrex::Real* operator()(const amrex::Real* inputs, int idx)
  {
    return (*this)(inputs) + idx;
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual const amrex::Real* operator()(const amrex::Real* inputs, int istart, int iend)
  {
    return (*this)(inputs) + istart;
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void operator()(const amrex::Real* inputs, amrex::Real* outputs)
  {
    int ninp = inpsize;
  
    if(num_layers == 1)
    {
      (*layers[0])(inputs, outputs, ninp);
    }
    else
    {
      amrex::Real* invec = outvec0;
      amrex::Real* outvec = outvec1;
      amrex::Real* temp;
      
      (*layers[0])(inputs, outvec, ninp);
      ninp = lsizes[0];
      
      for(int i = 1; i < num_layers-1; i++)
      {
        // Swap input and output buffers
        temp = invec;
        invec = outvec;
        outvec = temp;
        // Evaluate layer
        (*layers[i])(invec, outvec, ninp);
        ninp = lsizes[i];
      }
      
      invec = outvec;
      // Evaluate layer
      (*layers[num_layers-1])(invec, outputs, ninp);
    }
  }
  
  virtual void print()
  {
    amrex::Print() << "Number of layers: " << num_layers << std::endl;
    amrex::Print() << "Input size: " << inpsize << std::endl;
    amrex::Print() << std::endl;
    amrex::Print() << "Outgoing Sizes:" << std::endl;
    for(int i = 0; i < num_layers; i++)
    {
      amrex::Print() << layers[i]->nout << std::endl;
    }
    amrex::Print() << std::endl;
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual int nin()
  {
    return inpsize;
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual int nout()
  {
    return lsizes[num_layers-1];
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual int nfparams()
  {
    return nreals;
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual int niparams(bool with_prefix=true)
  {
    if(with_prefix)
    {
      return nints;
    }
    // Prefix length is also used in constructor, `pack` method
    return nints-3;
  }
  
  virtual void pack(amrex::Real* nnrd, int* nnid, bool store_prefix=true)
  {
    if(cmlm_net)
    {
      if(store_prefix)
      {
        // Store number of models, cmlm_net flag and input size
        // Edit constructor and `niparams` if you change the prefix length
        nnid[0] = 1;
        nnid[1] = static_cast<int>(cmlm_net);
        nnid[2] = inpsize;
        nnid += 3;
      }
      nnid[0] = nreals;
      nnid[1] = nints;
      nnid += 2;
      
      // Pack the layers into the arrays
      for(int i = 0; i < num_layers; i++)
      {
        layers[i]->pack(nnrd, nnid);
        nnrd += layers[i]->nfparams();
        nnid += layers[i]->niparams();
      }
    }
  }

private:
  
  NNLayer** layers;
  int* lsizes;
  int inpsize;
  int num_layers;
  amrex::Real* outvec0;
  amrex::Real* outvec1;
  bool cmlm_net;
  int nreals;
  int nints;
};

class NNArray: public NNModelBase
{
public:
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  NNArray() {}
  
  AMREX_GPU_HOST_DEVICE
  ~NNArray()
  {
    for(int i = 0; i < num_models; i++)
    {
      delete models[i];
    }
    delete[] models;
    delete[] out_idx;
    delete[] idx_map;
    delete[] outvec;
  }
  
  virtual void init_from_stream(std::ifstream& file, const InitParams& ipar)
  {
    cmlm_net = ipar.cmlm_net;
    inpsize = ipar.inpsize;
    num_models = ipar.nnets;
    
    int outsize = 0;
    models = new NNModel*[num_models];
    out_idx = new int[num_models+1];
    
    for(int i = 0; i < num_models; i++)
    {
      models[i] = new NNModel();
      models[i]->init_from_stream(file, ipar);
      out_idx[i] = outsize;
      outsize += models[i]->nout();
    }
    out_idx[num_models] = outsize;
    
    idx_map = new int[outsize];
    for(int i = 0; i < num_models; i++)
    {
      for(int j = out_idx[i]; j < out_idx[i+1]; j++)
      {
        idx_map[j] = i;
      }
    }
    
    outvec = new amrex::Real[outsize];
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void init_from_buffers(amrex::Real* nnrd, int* nnid, const InitParams& ipar)
  {
    num_models = ipar.nnets;  
    cmlm_net = ipar.cmlm_net;
    inpsize = ipar.inpsize;
    
    int outsize = 0;
    models = new NNModel*[num_models];
    out_idx = new int[num_models+1];
    
    for(int i = 0; i < num_models; i++)
    {
      models[i] = new NNModel();
      models[i]->init_from_buffers(nnrd, nnid, ipar);
      
      nnrd += models[i]->nfparams();
      nnid += models[i]->niparams(false);
      
      out_idx[i] = outsize;
      outsize += models[i]->nout();
    }
    out_idx[num_models] = outsize;
    
    idx_map = new int[outsize];
    for(int i = 0; i < num_models; i++)
    {
      for(int j = out_idx[i]; j < out_idx[i+1]; j++)
      {
        idx_map[j] = i;
      }
    }
    
    outvec = new amrex::Real[outsize];
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual const amrex::Real* operator()(const amrex::Real* inputs)
  {
    for(int i = 0; i < num_models; i++)
    {
      (*models[i])(inputs, outvec + out_idx[i]);
    }
    
    return outvec;
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual const amrex::Real* operator()(const amrex::Real* inputs, int idx)
  {
    int loc_idx = idx - out_idx[idx_map[idx]];
    const amrex::Real* out = (*models[idx_map[idx]])(inputs, loc_idx);
    outvec[idx] = *out;
    return outvec + idx;
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual const amrex::Real* operator()(const amrex::Real* inputs, int istart, int iend)
  {
    int mod_istart = idx_map[istart];
    int mod_iend = idx_map[iend];
    
    for(int i = mod_istart; i <= mod_iend; i++)
    {
      (*models[i])(inputs, outvec + out_idx[i]);
    }
    
    return outvec + istart;
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void operator()(const amrex::Real* inputs, amrex::Real* outputs)
  {
    for(int i = 0; i < num_models; i++)
    {
      (*models[i])(inputs, outputs + out_idx[i]);
    }
  }
  
  virtual void print()
  {
    for(int i = 0; i < num_models; i++)
    {
      models[i]->print();
      amrex::Print() << std::endl;
    }
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual int nin()
  {
    return inpsize;
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual int nout()
  {
    return out_idx[num_models-1];
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual int nfparams()
  {
    int nreals = 0;
    for(int i = 0; i < num_models; i++)
    {
      nreals += models[i]->nfparams();
    }
    return nreals;
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual int niparams(bool with_prefix=true)
  {
    int nints = 0;
    for(int i = 0; i < num_models; i++)
    {
      nints += models[i]->niparams(false);
    }
    
    if(with_prefix)
    {
      return nints+3;
    }
    return nints;
  }
  
  virtual void pack(amrex::Real* nnrd, int* nnid, bool store_prefix=true)
  {
    if(store_prefix)
    {
      // Store number of nets, cmlm_net flag and input size
      nnid[0] = num_models;
      nnid[1] = static_cast<int>(cmlm_net);
      nnid[2] = inpsize;
      nnid += 3;  
    }
    
    for(int i = 0; i < num_models; i++)
    {
      // The `false`s are because we have already stored the prefix data
      models[i]->pack(nnrd, nnid, false);
      nnrd += models[i]->nfparams();
      nnid += models[i]->niparams(false);
    }
  }

private:
  
  NNModel** models;
  int inpsize;
  int num_models;
  int* out_idx;
  int* idx_map;
  amrex::Real* outvec;
  bool cmlm_net;
};

// Factory functions for creating and initializing new NNModelBase instance on the heap
namespace NNModelFactory
{
  AMREX_FORCE_INLINE
  NNModelBase* from_ifstream(std::ifstream& file, bool cmlm_net=false)
  {
    InitParams ipar;
    file.read(reinterpret_cast<char*>(&ipar.ssize), sizeof(int));
    file.read(reinterpret_cast<char*>(&ipar.fsize), sizeof(int));
    file.read(reinterpret_cast<char*>(&ipar.isize), sizeof(int));
    file.read(reinterpret_cast<char*>(&ipar.inpsize), sizeof(int));
    file.read(reinterpret_cast<char*>(&ipar.nnets), sizeof(int));
    ipar.cmlm_net = cmlm_net;
    
    AMREX_ALWAYS_ASSERT(ipar.nnets > 0);
    
    NNModelBase* mptr;
    if(ipar.nnets == 1)
    {
      mptr = new NNModel();
    }
    else
    {
      mptr = new NNArray();
    }
    mptr->init_from_stream(file, ipar);
    
    return mptr;
  }
  
  AMREX_FORCE_INLINE
  NNModelBase* from_filename(std::string& fname, bool cmlm_net=false)
  {
    std::ifstream file(fname, std::ios::binary | std::ios::in);
    
    if(!file.is_open())
    {
      throw std::runtime_error("Unable to open input file " + fname + " .");
    }
    
    file.seekg(0, std::ios::beg);
    
    NNModelBase* model = from_ifstream(file, cmlm_net);
    file.close();
    return model;
  }
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  NNModelBase* from_buffers(amrex::Real* nnrd, int* nnid)
  {
    InitParams ipar;
    ipar.nnets = nnid[0];
    ipar.cmlm_net = nnid[1];
    ipar.inpsize = nnid[2];
    
    NNModelBase* mptr;
    if(ipar.nnets == 1)
    {
      mptr = new NNModel();
    }
    else
    {
      mptr = new NNArray();
    }
    mptr->init_from_buffers(nnrd, nnid+3, ipar);
    
    return mptr;
  }
} // namespace NNModelFactory
} // namespace physics
} // namespace pele

#endif
