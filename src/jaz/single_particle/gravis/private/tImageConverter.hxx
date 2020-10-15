/******************************************************************************
**        Title: tImageConverter.hxx
**  Description: Traits and functions used to convert different pixel types.
**               Required by reader/writer classes. INTERNAL USE ONLY !
**
**       Author: Jean-Sebastien Pierrard, 2005
**               Computer Science Department, University Basel (CH)
**
******************************************************************************/
namespace gravis
{
  namespace priv
  {

    template <class T>
    struct IPT_Traits { };

    template <>
    struct IPT_Traits<unsigned char>
    {
      typedef unsigned char Value_t;
      typedef int           Promote_t;
      typedef float         RealPromote_t;

      static inline Value_t fullValue ()
      {
        return 255;
      }
      static inline Value_t zeroValue ()
      {
        return   0;
      }
    };


    template <>
    struct IPT_Traits<float>
    {
      typedef float Value_t;
      typedef float Promote_t;
      typedef float RealPromote_t;

      static inline Value_t fullValue ()
      {
        return 1.0f;
      }
      static inline Value_t zeroValue ()
      {
        return 0.0f;
      }
    };


    template <>
    struct IPT_Traits<double>
    {
      typedef double Value_t;
      typedef double Promote_t;
      typedef double RealPromote_t;

      static inline Value_t fullValue ()
      {
        return 1.0;
      }
      static inline Value_t zeroValue ()
      {
        return 0.0;
      }
    };



    // Integral Pixel Type Traits, used to convert pixel components
    template <class fromT, class toT>
    struct IPC_Traits { };

#define _DeclareIPConversion(FROM_T, TO_T, CONVERSION) \
template<> struct IPC_Traits<FROM_T, TO_T> { \
	static inline TO_T convert (FROM_T value) { \
		return (TO_T)(CONVERSION); \
	} \
}

    _DeclareIPConversion(unsigned char, unsigned char, value);
    _DeclareIPConversion(unsigned char, float, (1.0f/255.0f)*float(value));
    _DeclareIPConversion(unsigned char, double, (1.0/255.0)*double(value));
    _DeclareIPConversion(double, unsigned char, 255.0f*value);
    _DeclareIPConversion(double, float, value);
    _DeclareIPConversion(float, unsigned char, 255.0f*value);
    _DeclareIPConversion(float, double, value);
    _DeclareIPConversion(float, float, value);
    _DeclareIPConversion(double, double, value);



    template <class T>
    struct pixelTypeConverter_1
    {
      static inline void f (const T& src, T& dst)
      {
        dst = src;
      }
    };


    template <class S, class D>
    struct pixelTypeConverter_2
    {
      static inline void f (const S& src, D& dst)
      {
        dst = IPC_Traits<S, D>::convert(src);
      };
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tRGBA<S>, tRGBA<D> >
    {
      static inline void f (const tRGBA<S>& src, tRGBA<D>& dst)
      {
        dst = tRGBA<D>(
                IPC_Traits<S, D>::convert(src.r),
                IPC_Traits<S, D>::convert(src.g),
                IPC_Traits<S, D>::convert(src.b),
                IPC_Traits<S, D>::convert(src.a)
              );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tRGB_A<S>, tRGB_A<D> >
    {
      static inline void f (const tRGB_A<S>& src, tRGB_A<D>& dst)
      {
        dst.set(
          IPC_Traits<S, D>::convert(src.r),
          IPC_Traits<S, D>::convert(src.g),
          IPC_Traits<S, D>::convert(src.b),
          IPC_Traits<S, D>::convert(src.a)
        );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tRGBA<S>, tRGB_A<D> >
    {
      static inline void f (const tRGBA<S>& src, tRGB_A<D>& dst)
      {
        dst.set(
          IPC_Traits<S, D>::convert(src.r),
          IPC_Traits<S, D>::convert(src.g),
          IPC_Traits<S, D>::convert(src.b),
          IPC_Traits<S, D>::convert(src.a)
        );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tGray_A<S>, tRGBA<D> >
    {
      static inline void f (const tGray_A<S>& src, tRGBA<D>& dst)
      {
        dst.set(
          IPC_Traits<S, D>::convert(src.g),
          IPC_Traits<S, D>::convert(src.g),
          IPC_Traits<S, D>::convert(src.g),
          IPC_Traits<S, D>::convert(src.a)
        );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tGray_A<S>, tRGB<D> >
    {
      static inline void f (const tGray_A<S>& src, tRGB<D>& dst)
      {
        dst.set(
          IPC_Traits<S, D>::convert(src.g),
          IPC_Traits<S, D>::convert(src.g),
          IPC_Traits<S, D>::convert(src.g)
        );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tRGB_A<S>, tRGBA<D> >
    {
      static inline void f (const tRGB_A<S>& src, tRGBA<D>& dst)
      {
        dst.set(
          IPC_Traits<S, D>::convert(src.r),
          IPC_Traits<S, D>::convert(src.g),
          IPC_Traits<S, D>::convert(src.b),
          IPC_Traits<S, D>::convert(src.a)
        );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tRGB<S>, tRGB<D> >
    {
      static inline void f (const tRGB<S>& src, tRGB<D>& dst)
      {
        dst = tRGB<D>(
                IPC_Traits<S, D>::convert(src.r),
                IPC_Traits<S, D>::convert(src.g),
                IPC_Traits<S, D>::convert(src.b)
              );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tRGB<S>, tRGBA<D> >
    {
      static inline void f (const tRGB<S>& src, tRGBA<D>& dst)
      {
        dst = tRGBA<D>(
                IPC_Traits<S, D>::convert(src.r),
                IPC_Traits<S, D>::convert(src.g),
                IPC_Traits<S, D>::convert(src.b),
                IPT_Traits<D>::fullValue()
              );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tRGB<S>, tRGB_A<D> >
    {
      static inline void f (const tRGB<S>& src, tRGB_A<D>& dst)
      {
        dst.set(
          IPC_Traits<S, D>::convert(src.r),
          IPC_Traits<S, D>::convert(src.g),
          IPC_Traits<S, D>::convert(src.b),
          IPT_Traits<D>::fullValue()
        );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tRGBA<S>, tRGB<D> >
    {
      static inline void f (const tRGBA<S>& src, tRGB<D>& dst)
      {
        dst = tRGB<D>(
                IPC_Traits<S, D>::convert(src.r),
                IPC_Traits<S, D>::convert(src.g),
                IPC_Traits<S, D>::convert(src.b)
              );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tRGB_A<S>, tRGB<D> >
    {
      static inline void f (const tRGB_A<S>& src, tRGB<D>& dst)
      {
        dst.set(
          IPC_Traits<S, D>::convert(src.r),
          IPC_Traits<S, D>::convert(src.g),
          IPC_Traits<S, D>::convert(src.b)
        );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<S, tRGBA<D> >
    {
      static inline void f (const S& src, tRGBA<D>& dst)
      {
        dst = tRGBA<D>(
                IPC_Traits<S, D>::convert(src),
                IPT_Traits<D>::fullValue()
              );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<S, tRGB_A<D> >
    {
      static inline void f (const S& src, tRGB_A<D>& dst)
      {
        dst.set(
          IPC_Traits<S, D>::convert(src),
          IPT_Traits<D>::fullValue()
        );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<S, tRGB<D> >
    {
      static inline void f (const S& src, tRGB<D>& dst)
      {
        dst = tRGB<D>(
                IPC_Traits<S, D>::convert(src)
              );
      }
    };


    // Convert from tRGB of type S to grayvalue of type D
    template <class S, class D>
    struct pixelTypeConverter_2<tRGB<S>, D>
    {
      static inline void f (const tRGB<S>& src, D& dst)
      {
        dst = IPC_Traits<S, D>::convert(src.grayValue());
      }
    };

    // Convert from tRGBA of type S to grayvalue of type D
    template <class S, class D>
    struct pixelTypeConverter_2<tRGBA<S>, D>
    {
      static inline void f (const tRGBA<S>& src, D& dst)
      {
        dst = IPC_Traits<S, D>::convert(src.grayValue());
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tRGB_A<S>, D>
    {
      static inline void f (const tRGB_A<S>& src, D& dst)
      {
        dst = IPC_Traits<S, D>::convert(src.grayValue());
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tGray_A<S>, D>
    {
      static inline void f (const tGray_A<S>& src, D& dst)
      {
        dst = IPC_Traits<S, D>::convert(src.grayValue());
      }
    };

    // tBGR

    template <class S, class D>
    struct pixelTypeConverter_2<tBGR<S>, tRGB<D> >
    {
      static inline void f (const tBGR<S>& src, tRGB<D>& dst)
      {
        dst = tRGB<D>(
                IPC_Traits<S, D>::convert(src.b),
                IPC_Traits<S, D>::convert(src.g),
                IPC_Traits<S, D>::convert(src.r)
              );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tRGB<S>, tBGR<D> >
    {
      static inline void f (const tRGB<S>& src, tBGR<D>& dst)
      {
        dst = tBGR<D>(
                IPC_Traits<S, D>::convert(src.b),
                IPC_Traits<S, D>::convert(src.g),
                IPC_Traits<S, D>::convert(src.r)
              );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<S, tBGR<D> >
    {
      static inline void f (const S& src, tBGR<D>& dst)
      {
        dst = tBGR<D>(
                IPC_Traits<S, D>::convert(src)
              );
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tRGBA<S>, tBGR<D> >
    {
      static inline void f (const tRGBA<S>& src, tBGR<D>& dst)
      {
        dst = tBGR<D>(
                IPC_Traits<S, D>::convert(src.b),
                IPC_Traits<S, D>::convert(src.g),
                IPC_Traits<S, D>::convert(src.r)
              );
      }
    };
    /*
    template <class S, class D>
    struct pixelTypeConverter_2<tRGB_A<S>, tBGR<D> > {
            static inline void f (const tRGB_A<S>& src, tBGR<D>& dst) {
                    dst.set(
                            IPC_Traits<S, D>::convert(src.b),
                            IPC_Traits<S, D>::convert(src.g),
                            IPC_Traits<S, D>::convert(src.r)
                    );
            }
    };
    */
    template <class S, class D>
    struct pixelTypeConverter_2<tBGR<S>, D >
    {
      static inline void f (const tBGR<S>& src, D& dst)
      {
        dst = IPC_Traits<S, D>::convert(src.b);
      }
    };

    template <class S, class D>
    struct pixelTypeConverter_2<tBGR<S>, tRGBA<D> >
    {
      static inline void f (const tBGR<S>& src, tRGBA<D>& dst)
      {
        dst = tRGBA<D>(
                IPC_Traits<S, D>::convert(src.b),
                IPC_Traits<S, D>::convert(src.g),
                IPC_Traits<S, D>::convert(src.r),
                IPT_Traits<D>::fullValue()
              );
      }
    };
    /*
    template <class S, class D>
    struct pixelTypeConverter_2<tBGR<S>, tRGB_A<D> > {
            static inline void f (const tBGR<S>& src, tRGB_A<D>& dst) {
                    dst.set(
                            IPC_Traits<S, D>::convert(src.b),
                            IPC_Traits<S, D>::convert(src.g),
                            IPC_Traits<S, D>::convert(src.r),
                            IPT_Traits<D>::fullValue()
                    );
            }
    };
    */
    template <class S, class D> inline
    void pixelTypeConverter (const S& src, D& dst)
    {
      pixelTypeConverter_2<S, D>::f(src, dst);
    }

    template <class T> inline
    void pixelTypeConverter	(const T& src, T& dst)
    {
      pixelTypeConverter_1<T>::f(src, dst);
    }



    template <class STYPE, class DTYPE>
    class tImageConverter
    {
      public:
        DTYPE operator() (const STYPE& src) const
        {

          DTYPE dst(src.cols(), src.rows());
          typename STYPE::iterator src_it = src.begin();
          typename DTYPE::iterator dst_it = dst.begin();

          for (; dst_it != dst.end(); ++src_it, ++dst_it)
          {
            pixelTypeConverter(*src_it, *dst_it);
          }
          return dst;
        }

        void convert(const STYPE& src, DTYPE& dst) const
        {

          dst.resize(src.cols(), src.rows());
          typename STYPE::iterator src_it = src.begin();
          typename DTYPE::iterator dst_it = dst.begin();

          for (; dst_it!=dst.end(); ++src_it, ++dst_it)
          {
            pixelTypeConverter(*src_it, *dst_it);
          }
        }
    };


  } /* Close Namespace "priv" */
} /* Close Namespace "gravis" */
