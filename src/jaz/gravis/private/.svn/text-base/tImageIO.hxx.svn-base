/******************************************************************************
**        Title: tImageIO.hxx
**  Description: Traits for tImage class.
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
    struct tImage_Traits { };

    template<class T>
    struct tImage_Traits<tImage<T> >
    {
      typedef float Pixel_t;
      typedef float CComp_t;
      static inline int nofComponents ()
      {
        return 1;
      }
    };

    template <class T>
    struct tImage_Traits<tImage<tRGB<T> > >
    {
      typedef tRGB<T> Pixel_t;
      typedef T       CComp_t;
      static inline int nofComponents ()
      {
        return 3;
      }
    };

    template <class T>
    struct tImage_Traits<tImage<tRGBA<T> > >
    {
      typedef tRGBA<T> Pixel_t;
      typedef T        CComp_t;
      static inline int nofComponents ()
      {
        return 4;
      }
    };

    template <class T>
    struct tImage_Traits<tImage<tRGB_A<T> > >
    {
      typedef tRGB_A<T> Pixel_t;
      typedef T         CComp_t;
      static inline int nofComponents ()
      {
        return 4;
      }
    };
    /*
    template <class T>
    struct tImage_Traits<tImage<tGray_A<T> > > {
    	typedef tGray_A<T> Pixel_t;
    	typedef T          CComp_t;
    	static inline int nofComponents () { return 2; }
    };
    */
  } /* Close namespace "priv" */
} /* Close namespace "gravis" */
