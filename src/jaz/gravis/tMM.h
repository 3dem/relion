#ifndef __LIBGRAVIS_T_MM_H__
#define __LIBGRAVIS_T_MM_H__

#include "tVarMatrix.h"
#include <iostream>
#include <vector>

namespace gravis
{

  template <class T>
  class tConstMM;

  template <class T>
  class tMM
  {
    public:
      typedef tMatrix<T, 3, 1> Vector;

      std::string title;

      tVarVector< Vector > nu;
      tVarMatrix< Vector > D;

      const size_t& m() const
      {
        return D.h;
      }; // Number of vertices
      const size_t& k() const
      {
        return D.w;
      }; // Number of paramters

      tMM(const std::string& title="Morphable Model") : title(title), nu(title+"::nu"), D(title+"::D") { };

      inline
      void evaluate(tVectorView< Vector > &v, const tConstVectorView< T > &a) const
      {
        tConstMM<T> MM(*this);
        MM.evaluate(v, a);
      }

      inline
      void evaluate(tVarVector< Vector > &v, const tConstVectorView< T > &a) const
      {
        tConstMM<T> MM(*this);
        MM.evaluate(v, a);
      }

      inline
      void resize(size_t h, size_t w)
      {
        D.resize(h,w);
        nu.resize(h);
      };

      inline
      void clear()
      {
        matrix::clear(D);
        matrix::clear(nu);
      };

      // Check if the file to load has the right datatype
      bool load_is_compatible(const std::string& fn)
      {
        char mmid0[33] = "GRAVIS_MORPHABLE_MODEL          ";
        char mmid1[33] = "                                ";
        std::ifstream stream(fn.c_str(), std::ifstream::binary);
        uint8_t uint32_size;
        uint8_t T_size;
        uint32_t m,k;
        uint16_t endianness;
        stream.read(mmid1, 32);
        stream.read((char*)&endianness, 2);
        stream.read((char*)&uint32_size, 1);
        stream.read((char*)&T_size, 1);
        stream.read((char*)&m, sizeof(m));
        stream.read((char*)&k, sizeof(k));
        GRAVIS_CHECK( 0 == strncmp( mmid0, mmid1, 31 ),   "Not a gravis morphable model file" );
        GRAVIS_CHECK( endianness  == 0x0001,         "Wrong endianness");
        if (uint32_size != 4)
        {
          std::cerr << "Uint 32 size is " << uint32_size << std::endl;
        }
        GRAVIS_CHECK( uint32_size == 4,              "Wrong uint32_size size");
        return( T_size      == sizeof(T) );
      }

      void load(const std::string& fn)
      {
        char mmid0[33] = "GRAVIS_MORPHABLE_MODEL          ";
        char mmid1[33] = "                                ";
        std::ifstream stream(fn.c_str(), std::ifstream::binary);
        uint8_t uint32_size;
        uint8_t T_size;
        uint32_t m,k;
        uint16_t endianness;
        stream.read(mmid1, 32);
        stream.read((char*)&endianness, 2);
        stream.read((char*)&uint32_size, 1);
        stream.read((char*)&T_size, 1);
        stream.read((char*)&m, sizeof(m));
        stream.read((char*)&k, sizeof(k));
        GRAVIS_CHECK( 0 == strncmp( mmid0, mmid1, 31 ),   "Not a gravis morphable model file" );
        GRAVIS_CHECK( endianness  == 0x0001,         "Wrong endianness");
        if (uint32_size != 4)
        {
          std::cerr << "Uint 32 size is " << uint32_size << std::endl;
        }
        GRAVIS_CHECK( uint32_size == 4,              "Wrong uint32_size size");
        GRAVIS_CHECK( T_size      == sizeof(T),      "Wrong type in model file");
        resize(m, k);
        clear();
        stream.read((char*)D.data, sizeof(Vector)*D.size());
        stream.read((char*)nu.data, sizeof(Vector)*nu.size());
        char mmid2[33] = "                                ";
        stream.read(mmid2, 32);
        GRAVIS_CHECK( 0 == strncmp( mmid0, mmid2, 31 ),   "File did not end with the end marker" );
      }

      void save(const std::string& fn)
      {
        tConstMM<T> cm(*this);
        cm.save(fn);
      }

      // Create a new interpolated model from barycentric coordinates into the old model
      inline
      void interpolate(tMM<T> &out, const tConstMatrixView<size_t> &idx, const tConstMatrixView<T> &weight) const
      {
        tConstMM<T> cm(*this);
        cm.interpolate(out, idx, weight);
      }

      // Create a new model from chosen lines of the old model
      inline
      void submodel(tMM<T> &out,  const tConstVectorView<bool> &chosen) const
      {
        tConstMM<T> cm(*this);
        cm.submodel(out, chosen);
      }

      // Create a new model from chosen lines of the old model
      inline
      void submodel(tMM<T> &out,  const tConstVectorView<size_t> &chosen) const
      {
        tConstMM<T> cm(*this);
        cm.submodel(out, chosen);
      }

      // Create a new model from chosen lines of the old model
      inline
      void submodel(tMM<T> &out,  const std::vector<size_t> &chosen) const
      {
        tConstMM<T> cm(*this);
        cm.submodel(out, chosen);
      }

      // Create a new model from chosen lines of the old model
      inline
      void submodel(tMM<T> &out,  const std::vector<int> &chosen) const
      {
        tConstMM<T> cm(*this);
        cm.submodel(out, chosen);
      }

      // Create a new interpolated model from barycentric coordinates into the old model
      inline
      void interpolate(tMM<T> &out, const tConstMatrixView<size_t> &idx, const tConstMatrixView<T> &weight, const size_t& n_coeff) const
      {
        tConstMM<T> cm(*this);
        cm.interpolate(out, idx, weight, n_coeff);
      }

      // Create a new interpolated model from barycentric coordinates into the old model
      inline
      void interpolate(tMM<T> &out, const tConstMatrixView<int> &idx, const tConstMatrixView<T> &weight, const int& n_coeff) const
      {
        tConstMM<T> cm(*this);
        cm.interpolate(out, idx, weight, n_coeff);
      }

      // Create a new model from chosen lines of the old model
      inline
      void submodel(tMM<T> &out,  const tConstVectorView<bool> &chosen, const size_t& n_coeff) const
      {
        tConstMM<T> cm(*this);
        cm.submodel(out, chosen, n_coeff);
      }

      // Create a new model from chosen lines of the old model
      inline
      void submodel(tMM<T> &out,  const tConstVectorView<size_t> &chosen, const size_t& n_coeff) const
      {
        tConstMM<T> cm(*this);
        cm.submodel(out, chosen, n_coeff);
      }

      // Create a new model from chosen lines of the old model
      inline
      void submodel(tMM<T> &out,  const std::vector<size_t> &chosen, const size_t& n_coeff) const
      {
        tConstMM<T> cm(*this);
        cm.submodel(out, chosen, n_coeff);
      }

      // Create a new model from chosen lines of the old model
      inline
      void submodel(tMM<T> &out,  const std::vector<int> &chosen, const int& n_coeff) const
      {
        tConstMM<T> cm(*this);
        cm.submodel(out, chosen, n_coeff);
      }
  };

  template <class T>
  class tConstMM
  {
    public:
      typedef tMatrix<T, 3, 1> Vector;

      const tConstVectorView< Vector > nu;
      const tConstMatrixView< Vector > D;

      const size_t& m() const
      {
        return D.h;
      }; // Number of vertices
      const size_t& k() const
      {
        return D.w;
      }; // Number of paramters

      tConstMM(const tConstVectorView<Vector> &nu, const tConstMatrixView<Vector> &D) : nu(nu), D(D)
      {
        GRAVIS_CHECK( nu.h == D.h, "Morphable model is inconsistent" );
      };
      tConstMM(const tMM<T>      &o) : nu(o.nu), D(o.D) {};
      tConstMM(const tConstMM<T> &o) : nu(o.nu), D(o.D) {};
#ifdef MATLAB
      tConstMM(const tmxConstMatrixView<T,1> &_nu, const tmxConstMatrixView<T,2> &_D) :
        nu((Vector*)_nu.data, _nu.dims[0]/3) ,
        D((Vector*)_D.data, _D.dims[0]/3, _D.dims[1])
      {
        GRAVIS_CHECK( nu.h == D.h, "Morphable model is inconsistent" );
      };
#endif

      // If a is too short assumes zeros for the unset coefficents
      // If a is too long assumes zeros for the missing principal components
      inline
      void evaluate(tVectorView< Vector > &v, const tConstVectorView< T > &a) const
      {
        size_t K=std::min(k(), a.size());
        //            GRAVIS_CHECK( a.size()  == k(), "a and D are incompatible");
        GRAVIS_CHECK( v.size()  == m(), "v and nu are incompatible");
        GRAVIS_CHECK( nu.size() == m(), "k and nu are incompatible");
        // Apply the morphable model
#if 0
        v = nu;
        for (size_t j=0; j<k(); ++j)
          for (size_t i=0; i<m(); ++i)                 // Multiplication of data
            for (size_t d=0; d<3; ++d)  // TODO: Use BLAS
            {
              v[i][d] += a[j]*D(i,j)[d];
            }
#endif
        tVectorView< T > vv((T*)v.data, 3*v.size());
        tConstVectorView< T > vnu((T*)nu.data, 3*nu.size());
        tConstMatrixView< T > vD((T*)D.data, 3*D.h, K);
        tConstVectorView< T > va(a.data, K);
        matrix::addmult(vv, vnu, vD, va); // USING BLAS
      }

      inline
      void evaluate(tVarVector< Vector > &v, const tConstVectorView< T > &a) const
      {
        v.resize(m());
        tVectorView< tMatrix<T, 3, 1 > > vv(v);
        evaluate(vv, a);
      }

      void save(const std::string& fn)
      {
        char mmid[33] = "GRAVIS_MORPHABLE_MODEL          ";
        std::ofstream stream(fn.c_str(), std::ofstream::binary);
        uint8_t uint32_size = sizeof(uint32_t);
        if (uint32_size != 4)
        {
          std::cerr << "Uint 32 size is " << uint32_size << std::endl;
        }
        uint8_t T_size      = sizeof(T);
        uint32_t m_ = m(), k_ = k();
        uint16_t endianness = 0x0001;
        stream.write(mmid, 32);
        stream.write((char*)&endianness, 2);
        stream.write((char*)&uint32_size, 1);
        stream.write((char*)&T_size, 1);
        stream.write((char*)&(m_), sizeof(m_));
        stream.write((char*)&(k_), sizeof(k_));
        stream.write((char*)D.data, sizeof(Vector)*D.size());
        stream.write((char*)nu.data, sizeof(Vector)*nu.size());
        stream.write(mmid, 32);
      }

      // Create a new interpolated model from barycentric coordinates into the old model
      inline
      void interpolate(tMM<T> &out, const tConstMatrixView<size_t> &idx, const tConstMatrixView<T> &weight) const
      {
        const tConstMM<T> &model = *this;

        GRAVIS_CHECK( idx.w == weight.w && idx.h == weight.h, "idx and weight should be kxn and kxn");

        out.resize(idx.w, model.k());

        const size_t& n = idx.w;
        const size_t& t = idx.h;

        const size_t& K = model.k();

        // Initialize to zero
        out.clear();

        // Write out.nu
        for (size_t i=0; i<n; ++i)  // Run over output points
          for (size_t l=0; l<t; ++l)  // Run over rows to interpolate
            out.nu[i] += weight(l,i)*model.nu[idx(l,i)];
        // Write out.D
        for (size_t j=0; j<K; ++j)  // Run over output columns
          for (size_t i=0; i<n; ++i)  // Run over output points
            for (size_t l=0; l<t; ++l)  // Run over rows to interpolate
              out.D(i,j) += weight(l,i)*model.D(idx(l,i), j);
      }

      // Create a new model from chosen lines of the old model
      inline
      void submodel(tMM<T> &out,  const tConstVectorView<bool> &chosen) const
      {
        const tConstMM<T> &model = *this;


        GRAVIS_CHECK( chosen.h == model.m(), "Chosen and model are incompatible");

        size_t n = 0;
        for (size_t i=0; i<chosen.h; ++i)
          if (chosen[i])
            n += 1;

        out.resize(n, model.k());

        size_t I = 0;
        for (size_t i=0; i<model.m(); ++i)
          if (chosen[i])
          {
            out.nu[I] = model.nu[i];
            for (size_t j=0; j<model.k(); ++j)
              out.D(I,j) = model.D(i,j);
            ++I;
          }
      }

      // Create a new model from chosen lines of the old model
      inline
      void submodel(tMM<T> &out,  const tConstVectorView<size_t> &chosen) const
      {
        const tConstMM<T> &model = *this;

        size_t n = chosen.h;

        out.resize(n, model.k());

        for (size_t I=0; I<n; ++I)
        {
          size_t i = chosen[I];
          GRAVIS_CHECK(i < model.nu.h, "Accessing invalid position.");
          // std::cout << "Copy " << model.nu[i] << " from model.nu[" << i << "] to out.nu[" << I << "]" << std::endl;
          out.nu[I] = model.nu[i];
          for (size_t j=0; j<model.k(); ++j)
            out.D(I,j) = model.D(i,j);
        }
      }
      // Create a new model from chosen lines of the old model
      inline
      void submodel(tMM<T> &out,  const tConstVectorView<int> &chosen) const
      {
        const tConstMM<T> &model = *this;

        int n = chosen.h;

        out.resize(n, model.k());

        for (int I=0; I<n; ++I)
        {
          int i = chosen[I];
          GRAVIS_CHECK(i < model.nu.h, "Accessing invalid position.");
          // std::cout << "Copy " << model.nu[i] << " from model.nu[" << i << "] to out.nu[" << I << "]" << std::endl;
          out.nu[I] = model.nu[i];
          for (int j=0; j<model.k(); ++j)
            out.D(I,j) = model.D(i,j);
        }
      }


      // Create a new model from chosen lines of the old model
      inline
      void submodel(tMM<T> &out,  const std::vector<size_t> &chosen) const
      {
        tConstVectorView<size_t> vchosen(&chosen[0], chosen.size());
        submodel(out, vchosen);
      }
      // Create a new model from chosen lines of the old model
      inline
      void submodel(tMM<T> &out,  const std::vector<int> &chosen) const
      {
        tConstVectorView<int> vchosen(&chosen[0], chosen.size());
        submodel(out, vchosen);
      }

      // Create a new interpolated model from barycentric coordinates into the old model
      template <class Int>
      inline
      void interpolate(tMM<T> &out, const tConstMatrixView<Int> &idx, const tConstMatrixView<T> &weight, const Int& n_coeff) const
      {
        const tConstMM<T> &model = *this;

        GRAVIS_CHECK( idx.w == weight.w && idx.h == weight.h, "idx and weight should be kxn and kxn");

        out.resize(idx.w, n_coeff);

        const Int& n = idx.w;
        const Int& t = idx.h;

        const Int K = std::min<int>(n_coeff, model.k());

        // Initialize to zero
        out.clear();

        // Write out.nu
        for (Int i=0; i<n; ++i)  // Run over output points
          for (Int l=0; l<t; ++l)  // Run over rows to interpolate
            out.nu[i] += weight(l,i)*model.nu[idx(l,i)];
        // Write out.D
        for (Int j=0; j<K; ++j)  // Run over output columns
          for (Int i=0; i<n; ++i)  // Run over output points
            for (Int l=0; l<t; ++l)  // Run over rows to interpolate
              out.D(i,j) += weight(l,i)*model.D(idx(l,i), j);
      }

      // Create a new model from chosen lines of the old model
      inline
      void submodel(tMM<T> &out,  const tConstVectorView<bool> &chosen, const size_t& n_coeff) const
      {
        const tConstMM<T> &model = *this;


        GRAVIS_CHECK( chosen.h == model.m(), "Chosen and model are incompatible");

        size_t n = 0;
        for (size_t i=0; i<chosen.h; ++i)
          if (chosen[i])
            n += 1;

        out.resize(n, n_coeff);

        const size_t copy_coeff = std::min(n_coeff, model.k());

        size_t I = 0;
        for (size_t i=0; i<model.m(); ++i)
          if (chosen[i])
          {
            out.nu[I] = model.nu[i];
            for (size_t j=0; j<copy_coeff; ++j)
              out.D(I,j) = model.D(i,j);
            ++I;
          }
      }

      // Create a new model from chosen lines of the old model
      template <class Int>
      inline
      void submodel(tMM<T> &out,  const tConstVectorView<Int> &chosen, const Int& n_coeff) const
      {
        const tConstMM<T> &model = *this;

        Int n = chosen.h;

        out.resize(n, n_coeff);

        const Int copy_coeff = std::min<int>(n_coeff, model.k());

        for (Int I=0; I<n; ++I)
        {
          Int i = chosen[I];
          GRAVIS_CHECK(i < model.nu.h, "Accessing invalid position.");
          // std::cout << "Copy " << model.nu[i] << " from model.nu[" << i << "] to out.nu[" << I << "]" << std::endl;
          out.nu[I] = model.nu[i];
          for (Int j=0; j<copy_coeff; ++j)
            out.D(I,j) = model.D(i,j);
        }
      }

      // Create a new model from chosen lines of the old model
      template <class Int>
      inline
      void submodel(tMM<T> &out,  const std::vector<Int> &chosen, const Int& n_coeff) const
      {
        tConstVectorView<Int> vchosen(&chosen[0], chosen.size());
        submodel(out, vchosen, n_coeff);
      }
  };

}

#endif
